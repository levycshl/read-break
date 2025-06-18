# Read-Break

Declarative **read-parsing and clipping** for paired-end FASTQ data.  
You describe *what* to extract or test in a YAML spec; `ReadParser` executes the steps, logs
success/failure statistics, and (optionally) writes trimmed reads via `read_clip_and_write`.

* No global state: pure functions where possible  
* Streaming data: constant memory, gzip-in / gzip-out  
* Embed within larger pipelines (e.g. [`bag-pipe`](https://github.com/levycshl/bag-pipe))
* Enables functions like `match`, `extract`, `hamming_test`, `compute`, `regex_search`, and more.

## Install

```bash
pip install read-break
# or: pip install git+https://github.com/levycshl/read-break.git
````

Python ≥ 3.11 required.

---

## Quick start

```bash
# A minimal run that clips reads per examples/test_pipeline.yaml
python -m read_break.cli \
    --config examples/test_pipeline.yaml \
    --r1   sample.R1.fastq.gz \
    --r2   sample.R2.fastq.gz \
    --out  clipped/
```

*See `parsers/` for more specs, including the full **bagpipe\_2N.yaml** that exercises almost every operation.*

---

## Pipeline specification

The pipeline is defined in a YAML file, with two main sections: `params` and `pipeline`.
The parameters section defines global variables, while the pipeline section describes the steps to execute.
Any statement can be templated using Jinja syntax, allowing dynamic values based on the input reads.

Each step also contains a must_pass flag, which determines whether the read should fail and stop processing if the step does not pass.

For example:
```yaml
params:
  # Adapters 
  adapter_r1: "CCAAACACACCCAA"
  adapter_r2: "GAGCGGACTCTGCG"

  # Upstream motif (Read 1)
  up1_seq: "CATG"
  up1_start: 63
  up1_end: 67
  up1_len: "{{ params.up1_end - params.up1_start }}" # Calculated length
  up1_max_hamming: 1...

pipeline:
    # 1. Check upstream motif on Read 1
  - id: check_upstream_motif_r1
    description: Verify upstream motif in R1 within Hamming distance
    read: 1
    op: hamming_test # Uses hamming_test operation
    ref: "{{ params.up1_seq }}"
    start: "{{ params.up1_start }}"
    length: "{{ params.up1_len }}"
    max_mismatch: "{{ params.up1_max_hamming }}"
    hamming_fn: "hamming" # Standard Hamming distance
    store_result_as: motif_ok # Store boolean result
    must_pass: true # Fail read if motif doesn't match
```


Every step is a mapping with at least:

| Field  | Type | Required | Description                            |
| ------ | ---- | -------- | -------------------------------------- |
| `id`   | str  | ✓        | Unique label for logs                  |
| `op`   | str  | ✓        | One of the operation names below       |
| `read` | int  | ✗        | `1` or `2` (omit for non-sequence ops) |

### `op: match`

Locate the first instance of `ref` starting at base_offset and 
searching up to `max_wobble` bp and allowing `max_mismatch` mismatches.

| Field          | Type       | Required | Description              |
| -------------- | ---------- | -------- | ------------------------ |
| `ref`          | str / tmpl | ✓        | Reference sequence       |
| `hamming_fn`   | str        | ✓        | Key from `HAMMING_FUNCS` |
| `max_wobble`   | int        | ✓        | Shift window             |
| `max_mismatch` | int        | ✓        | Allowed mismatches       |
| `base_offset`  | int        | –        | Start index (default 0)  |
| `store_pos_as` | str        | ✓        | Save offset here         |

### `op: extract`

Slice a fragment and optionally whitelist-check it.

| Field             | Type      | Required | Description                         |
| ----------------- | --------- | -------- | ----------------------------------- |
| `start`, `length` | str / int | ✓        | Can be templated                    |
| `store_seq_as`    | str       | ✓        | Variable name                       |
| `whitelist`       | str       | –        | Key in `globals.barcode_whitelists` |
| `store_match_as`  | str       | –        | Where to store the boolean          |

### `op: hamming_test`

Compare a read segment to `ref` using any registered Hamming-like metric.

| Field                                                  | Type | Required | Description             |
| ------------------------------------------------------ | ---- | -------- | ----------------------- |
| `ref`, `start`, `length`, `hamming_fn`, `max_mismatch` | –    | ✓        |                         |
| `store_result_as`                                      | str  | –        | Variable for True/False |

### `op: regex_search` 

Search `seq` with a pre-compiled pattern defined in `globals.regex_patterns`.

| Field            | Type | Required | Description             |
| ---------------- | ---- | -------- | ----------------------- |
| `pattern`        | str  | ✓        | Name of pattern         |
| `store_pos_as`   | str  | ✓        | Start index (or `None`) |
| `store_match_as` | str  | –        | Save the actual match   |

### `op: test` 

Evaluate an arbitrary Jinja expression that returns a Boolean.

| Field             | Type | Required | Description         |
| ----------------- | ---- | -------- | ------------------- |
| `expression`      | str  | ✓        | Jinja expression    |
| `store_result_as` | str  | –        | Variable for result |

### `op: compute` 

Render a Jinja expression and stash the result (number or string).

| Field        | Type | Required | Description                             |
| ------------ | ---- | -------- | --------------------------------------- |
| `expression` | str  | ✓        | Jinja expression                        |
| `store_as`   | str  | ✓        | Variable name                           |
| `pass_if`    | str  | –        | Boolean expression that decides success |

---

## Helper modules

| Module                | Highlights                                                                                                  |
| --------------------- | ----------------------------------------------------------------------------------------------------------- |
| `read_break.logic`    | `reverse_complement`, `seq_to_int`, plain `hamming` plus asymmetric variants (`hamming35`), `wobble_match`  |
| `read_break.registry` | `HAMMING_FUNCS` = `hamming`, `hammingTC` (T→C-blind), `hammingAG` (A→G-blind) – extendable via PRs          |
| `read_break.io`       | `FastqReader` / `FastqWriter` for streaming gzipped FASTQ; both are context-manager-friendly                |
| `read_break.parser`   | `ReadParser` (core engine) + `read_clip_and_write` convenience wrapper                                      |

---

## Logging & metrics

`ReadParser` keeps a `parse_log` you can inspect after a run:

```python
>>> parser.get_parse_log()
{'total_reads': 1_000_000,
 'successful_reads': 973_451,
 'failed_reads': 26_549,
 'success_rate': 97.35,
 'failures_by_step': {'s1_match': 11023, 'umi_wl': 5943, ...}}
```

---

## Roadmap / In-development

* **Failure piping** – a first-class API to stream failed-read records (with step ID & message) to file or Kafka for later QC.  The internal hooks exist; the external interface is being finalised.
* **CLI subcommand `validate`** – static checker for YAML specs.
* **More metrics** – Levenshtein-based matchers and indel-aware wobble.

---

## Citation

If this tool helps your research, please cite:

```
@software{read_break,
  author  = {Levy, Dan},
  title   = {Read-Break: Declarative parsing of paired-end reads},
  version = {0.1.0},
  url     = {https://github.com/levycshl/read-break},
  year    = {2025}
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
