# read_break/parser.py

import yaml
import jinja2
from typing import Dict, Any, Optional, Callable
import os
from read_break.logic import wobble_match
from read_break.registry import HAMMING_FUNCS
from read_break.io import FastqReader, FastqWriter
import re

from functools import lru_cache
import ast
from jinja2 import Environment, StrictUndefined, meta

ENV = Environment(undefined=StrictUndefined)   # module-level – create once

@lru_cache(maxsize=None)        # key = literal template string
def _compile(template: str):
    return ENV.from_string(template)


class ReadParser:
    """
    Core engine for executing a declarative read-parsing pipeline.

    Each pipeline step describes an operation (e.g., 'match', 'extract', 'hamming_test'),
    which operates on a FASTQ read (read 1 or read 2), and stores intermediate results
    in a shared context for later steps to use.
    """
    

    def __init__(self, pipeline_cfg: Dict[str, Any], globals_cfg: Optional[Dict[str, Any]] = None, 
                 globals_namespace: str = "params", base_dir: Optional[str] = None):
        # Store the list of pipeline steps, each a dict with operation and parameters
        self.steps = pipeline_cfg["pipeline"]

        # Jinja environment to render dynamic expressions from context (e.g., {{ s1_start + 15 }})
        #self.env = jinja2.Environment(undefined=jinja2.StrictUndefined)
        ## declared the Jinja environment at the module level # getting rid of self.env

        # Optional global constants (e.g. LT_LEN, sequences) passed in as 'cfg'
        
        # Namespace for accessing globals in templates (default: "cfg")
        self.globals_namespace = globals_namespace
        self.globals = globals_cfg.copy() if globals_cfg else {}
        # Base directory for resolving relative paths (e.g., whitelist files)
        self.base_dir = base_dir
        self._resolve_globals()
        
        # Load whitelist files during initialization
        self._load_whitelists()

        ## Load regex patterns
        self._load_regex_patterns()
        
        # Initialize parse log
        self.reset_parse_log()

        ## Freeze constants that reference only globals
        self._freeze_constants()

    def _freeze_constants(self) -> None:
        """
        Render any step-field that references **only** global params
        (i.e. does not depend on per-read context) exactly **once**.
        """
        for step in self.steps:
            for k, v in list(step.items()):
                if isinstance(v, str) and "{{" in v:
                    # Which variables does this template use?
                    ast_ = ENV.parse(v)
                    undeclared = meta.find_undeclared_variables(ast_)

                    # If all vars are from the globals namespace, it's constant
                    if undeclared <= {self.globals_namespace}:
                        step[k] = self._render(v, {self.globals_namespace: self.globals})

    def _load_regex_patterns(self) -> None:
        """Load regex patterns from the globals configuration."""
        self.regex_patterns = {}
        regex_configs = self.globals.get("regex_patterns", {})
        for pattern_name, pattern_dict in regex_configs.items():
            pattern_type = self._get_param(pattern_dict, "type", self.globals, default="full")
            sequence = self._get_param(pattern_dict, "sequence", self.globals)
            try:
                if pattern_type == "full":
                    compiled_pattern = re.compile(sequence)
                elif pattern_type == "full_or_tail":
                    min_tail = self._get_param(pattern_dict, "min_tail", self.globals, default=4, convert_type=int)
                    # Avoid issues if min_tail > len(sequence)
                    if min_tail > len(sequence):
                        raise ValueError(
                            f"min_tail ({min_tail}) is longer than sequence length ({len(sequence)}) for '{pattern_name}'")
                    alternates = [sequence[:i] for i in range(min_tail, len(sequence) + 1)]
                    partial = "(%s)$" % "|".join(alternates)
                    pattern_str = f"{sequence}|{partial}"
                    compiled_pattern = re.compile(pattern_str)
                else:
                    print(f"Unknown regex type '{pattern_type}' for '{pattern_name}'")
                    compiled_pattern = None
            except Exception as e:
                print(f"Error compiling regex pattern '{pattern_name}': {e}")
                compiled_pattern = None
            self.regex_patterns[pattern_name] = compiled_pattern

    def reset_parse_log(self) -> None:
        """Resets the parse log statistics."""
        self.parse_log = {
            "total_reads": 0,
            "successful_reads": 0,
            "failed_reads": 0,
            "failures_by_step": {step.get("id", f"step_{i}"): 0 for i, step in enumerate(self.steps)}
        }

    def _log_failure(self, step_id: str, reason: str) -> None:
        """Record a failure in the parse log."""
        # This method is deprecated and no longer used
        pass

    def get_parse_log(self) -> Dict[str, Any]:
        """Returns the current parse log with statistics."""
        log = dict(self.parse_log)
        if log["total_reads"] > 0:
            log["success_rate"] = round((log["successful_reads"] / log["total_reads"]) * 100, 2)
        return log

    def _resolve_globals(self) -> None:
        """
        Evaluate any Jinja templates that appear inside the globals dict
        (e.g. VT2_START: "{{ params.R2_S2 | length }}").
        Runs until a pass makes no further changes, handling simple
        dependencies without trying to detect cycles.
        """
        changed = True
        while changed:
            changed = False
            for key, val in list(self.globals.items()):
                rendered = self._render(val, {self.globals_namespace: self.globals})
                if rendered != val:
                    self.globals[key] = rendered
                    changed = True

    def _render(self, template_like: Any, ctx: Dict[str, Any]) -> Any:
        """
        Fast, safe rendering:
        • non-strings      → return as-is
        • plain strings    → return as-is
        • Jinja templates  → compile once (LRU), render every call
        • result           → literal-eval if it looks like a Python literal,
                            otherwise return the rendered string verbatim.
        """
        if not (isinstance(template_like, str) and "{{" in template_like):
            return template_like

        rendered = _compile(template_like).render(**ctx)
        try:
            return ast.literal_eval(rendered)  # safe numbers / lists / dicts / bools
        except (ValueError, SyntaxError):
            return rendered      # leave non-literal strings untouched


    def _get_param(self, step: Dict[str, Any], param_name: str, context: Dict[str, Any], 
                   default: Any = None, convert_type: Optional[Callable] = None) -> Any:
        """
        Centralized parameter retrieval with template rendering and type conversion.
        
        Args:
            step: The pipeline step dictionary
            param_name: Name of the parameter to retrieve
            context: Current execution context
            default: Default value if parameter is not found
            convert_type: Optional function to convert the result (e.g., int, str)
            
        Returns:
            The rendered and converted parameter value
        """
        # Get the parameter, or use default if not found
        if param_name not in step and default is None:
            raise KeyError(f"Required parameter '{param_name}' not found in step '{step.get('id', 'unknown')}'")
        
        value = step.get(param_name, default)
        
        # Render the value if it's a template
        rendered_value = self._render(value, {**context, self.globals_namespace: self.globals})
        
        # Convert the type if requested
        if convert_type is not None and rendered_value is not None:
            try:
                return convert_type(rendered_value)
            except (ValueError, TypeError) as e:
                raise ValueError(f"Failed to convert '{param_name}' to {convert_type.__name__}: {e}")
        
        return rendered_value

    def __str__(self) -> str:
        lines = ["ReadParser Pipeline:"]
        for i, step in enumerate(self.steps):
            step_id = step.get("id", f"step_{i}")
            op = step["op"]
            read_id = step.get("read", None)
            must_pass = step.get("must_pass", True)
            if op == 'test' or op == 'compute':
                lines.append(f"[{i:3} ]  {step_id}: {op} ({step['expression']}, must pass: {must_pass})")
            else:
                lines.append(f"[{i:3} ]  {step_id}: {op} (read {read_id}, must pass: {must_pass})")
        return "\n".join(lines)

    def to_yaml(self) -> str:
        return yaml.dump({"pipeline": self.steps}, sort_keys=False)

    def parse(self, name: str, seq1: str, qual1: str, seq2: str, qual2: str) -> Optional[Dict[str, Any]]:
        """
        Execute the parsing pipeline on a single read pair.
        Returns a context dictionary of extracted values if successful,
        or None if any step marked 'must_pass' fails.
        """
        # Select which read to operate on: read 1 or read 2
        reads = {1: seq1, 2: seq2}
        context: Dict[str, Any] = {"read_id": name,
                                   "len_seq1": len(seq1),
                                   "len_seq2": len(seq2)}

        # Update total reads count
        self.parse_log["total_reads"] += 1

        for step in self.steps:
            step_id = step.get("id", "unknown")
            if "read" not in step:
                read_id = 0
                seq = ""
            else:
                read_id = step["read"]
                seq = reads[read_id]
            op = step["op"]
            must_pass = step.get("must_pass", True)


            # Process based on operation type
            try:
                if op == "match":
                    success = self._process_match(step, seq, context)
                elif op == "extract":
                    success = self._process_extract(step, seq, context)
                elif op == "hamming_test":
                    success = self._process_hamming_test(step, seq, context)
                elif op == "test":
                    success = self._process_test(step, seq, context)
                elif op == "regex_search":
                    success = self._process_regex_search(step, seq, context)
                elif op == "compute":
                    success = self._process_compute(step, seq, context)
                else:
                    raise NotImplementedError(f"Unknown operation: {op}")

                if not success:
                    ## if this step fails, optional or not, record the failure
                    self.parse_log['failures_by_step'][step_id] += 1
                    ## but if it must pass, mark the read as failed and return
                    if must_pass:
                        # increment failed_reads once per read
                        self.parse_log["failed_reads"] += 1
                        return {
                            "read_id": name,
                            "status": "fail",
                            "failed_step": step_id,
                            "message": f"{op.capitalize()} operation failed"
                        }
                
            except Exception as e:
                if must_pass:
                    self.parse_log["failed_reads"] += 1
                    self.parse_log["failures_by_step"][step_id] += 1
                    return {
                        "read_id": name,
                        "status": "fail",
                        "failed_step": step_id,
                        "message": f"{op.capitalize()} operation error: {e}"
                    }
                # If step is optional, continue to next step

        # Add final required metadata
        context["status"] = "ok"
        self.parse_log["successful_reads"] += 1
        return context

    def _process_compute(self, step, seq, context):
        """
        Required keys:
          expression : Jinja template to evaluate
          store_as   : context key for the result
        Optional:
          pass_if    : Jinja boolean; step succeeds iff True
        """
        value = self._get_param(step, "expression", context)
        context[step["store_as"]] = value

        # Evaluate optional pass_if condition
        if "pass_if" in step:
            condition = self._get_param(step, "pass_if", context, convert_type=bool)
            return bool(condition)
        return True

    def _process_regex_search(self, step, seq, context):
        """
        Run a named, precompiled regex against seq.
        - Returns True if matched, False if not matched.
        """
        pattern_name = step["pattern"]
        rx = self.regex_patterns.get(pattern_name)
        if rx is None:
            raise KeyError(f"No regex named '{pattern_name}' in self.regex_patterns")
        m = rx.search(seq)
        if m is not None:
            context[step["store_pos_as"]] = m.start()
            ## only store the match if requested
            if "store_match_as" in step:
                context[step["store_match_as"]] = m.group()
            return True
        else:
            context[step["store_pos_as"]] = step.get("default", None)
            return False

    def _process_match(self, step: Dict[str, Any], seq: str, context: Dict[str, Any]) -> bool:
        """Process a match operation and return success status."""
        ref = self._get_param(step, "ref", context)
        max_wobble = self._get_param(step, "max_wobble", context, convert_type=int)
        max_hamming = self._get_param(step, "max_mismatch", context, convert_type=int)
        base_offset = self._get_param(step, "base_offset", context, default=0, convert_type=int)
        hamming_fn_name = self._get_param(step, "hamming_fn", context)
        hamming_fn = HAMMING_FUNCS[hamming_fn_name]

        offset = wobble_match(
            seq,
            ref,
            max_wobble=max_wobble,
            max_hamming=max_hamming,
            base_offset=base_offset,
            hamming_func=hamming_fn
        )

        if offset == -1:
            if step.get("must_pass", True):
                return False
            context[step["store_pos_as"]] = None
        else:
            context[step["store_pos_as"]] = offset
        return True

    def _process_extract(self, step: Dict[str, Any], seq: str, context: Dict[str, Any]) -> bool:
        """Process an extract operation and return success status."""
        start = self._get_param(step, "start", context, convert_type=int)
        length = self._get_param(step, "length", context, convert_type=int)
        frag = seq[start:start + length]
        context[step["store_seq_as"]] = frag
        is_ok = True
        # optional whitelist check
        if "whitelist" in step:
            wl_key = step["whitelist"]
            whitelist = self.whitelist_sets.get(wl_key, set())
            is_ok = frag in whitelist
            ## store the booklean result (default name: <id>_ok if not provided)
            match_name = step.get("store_match_as", f"{step['id']}_ok")
            context[match_name] = is_ok
        return is_ok

    def _process_hamming_test(self, step: Dict[str, Any], seq: str, context: Dict[str, Any]) -> bool:
        """Process a hamming_test operation and return success status."""
        ref = self._get_param(step, "ref", context)
        start = self._get_param(step, "start", context, convert_type=int)
        length = self._get_param(step, "length", context, convert_type=int)
        max_mismatch = self._get_param(step, "max_mismatch", context, convert_type=int)
        hamming_fn_name = self._get_param(step, "hamming_fn", context)
        hamming_fn = HAMMING_FUNCS[hamming_fn_name]

        target_seq = seq[start:start + length]
        distance = hamming_fn(ref, target_seq)        
        result = distance <= max_mismatch
        
        # Store boolean result under custom name or step ID
        context[step.get("store_result_as", step["id"])] = result
        
        return result or not step.get("must_pass", True)



    def _process_test(
        self, step: Dict[str, Any], seq: str, context: Dict[str, Any]
    ) -> bool:
        """
        Evaluate a Boolean expression that can reference any variables
        already present in `context` or in constants.

        Returns True if the test passes, False otherwise.
        """
        result = self._get_param(step, "expression", context)
        # Store the Boolean outcome (default key = step id)
        store_name = step.get("store_result_as", step["id"])
        context[store_name] = result
        return result  # True = pass, False = fail


    def _load_whitelists(self) -> None:
        """Load barcode whitelist files during initialization."""
        self.whitelist_sets = {}
        whitelist_configs = self.globals.get("barcode_whitelists", {})
        
        for bc_name, path in whitelist_configs.items():
            try:
                # Resolve the path relative to base_dir if provided
                if self.base_dir:
                    full_path = os.path.join(self.base_dir, path)
                else:
                    full_path = path

                with open(full_path, 'r') as f:
                    self.whitelist_sets[bc_name] = set(line.strip() for line in f if line.strip())
                print(f"Loaded {len(self.whitelist_sets[bc_name])} barcodes for {bc_name}")
            except Exception as e:
                print(f"Error loading whitelist for {bc_name}: {e}")
                self.whitelist_sets[bc_name] = set()





def read_clip_and_write(
    reader: FastqReader,
    parser: ReadParser,
    writer: FastqWriter,
    *,
    default_start_r1: int = 0,
    default_end_r1:   int = -1,
    default_start_r2: int = 0,
    default_end_r2:   int = -1,
    default_read_tag: str = ""
) -> None:
    """
    Iterate through `reader`; parse each pair; trim, tag, and write.

    Parser context keys looked for:
        start_r1, end_r1, start_r2, end_r2, read_tag
    Any missing key falls back to the corresponding default_* argument.

    An end position of -1 (or None) means "slice to end of read".
    """
    for read_id, seq1, qual1, seq2, qual2 in reader:
        ctx = parser.parse(read_id, seq1, qual1, seq2, qual2)
        if not ctx or ctx.get("status") != "ok":
            continue

        s1 = ctx.get("start_r1", default_start_r1)
        e1 = ctx.get("end_r1",   default_end_r1)
        s2 = ctx.get("start_r2", default_start_r2)
        e2 = ctx.get("end_r2",   default_end_r2)
        tag = ctx.get("read_tag", default_read_tag)

        # -1 or None means to the end of the read
        e1 = len(seq1) if (e1 is None or e1 == -1) else e1
        e2 = len(seq2) if (e2 is None or e2 == -1) else e2

        trimmed_seq1  = seq1 [s1:e1]
        trimmed_qual1 = qual1[s1:e1]
        trimmed_seq2  = seq2 [s2:e2]
        trimmed_qual2 = qual2[s2:e2]

        new_id = f"{read_id}/1_{tag}" if tag else f"{read_id}/1"
        writer.write(
            (new_id, trimmed_seq1, trimmed_qual1, trimmed_seq2, trimmed_qual2)
        )