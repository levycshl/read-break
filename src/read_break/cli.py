#!/usr/bin/env python3
"""
Command-line interface for read-break.

Provides the main CLI functionality to parse and clip paired-end FASTQ reads
based on a YAML configuration pipeline.
"""

import argparse
import sys
from pathlib import Path
import yaml

from read_break.parser import ReadParser, read_clip_and_write
from read_break.io import FastqReader, FastqWriter


def main(argv=None):
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Parse and clip paired-end FASTQ reads using declarative YAML pipeline"
    )
    parser.add_argument(
        "--config", 
        type=Path, 
        required=True, 
        help="YAML pipeline configuration file"
    )
    parser.add_argument(
        "--r1", 
        type=Path, 
        required=True, 
        help="Path to R1 FASTQ file (can be gzipped)"
    )
    parser.add_argument(
        "--r2", 
        type=Path, 
        required=True, 
        help="Path to R2 FASTQ file (can be gzipped)"
    )
    parser.add_argument(
        "--out", 
        type=Path, 
        required=True, 
        help="Output directory for clipped reads"
    )
    parser.add_argument(
        "--prefix", 
        default="clipped", 
        help="Prefix for output files (default: clipped)"
    )
    
    args = parser.parse_args(argv)
    
    # Validate inputs
    if not args.config.exists():
        sys.exit(f"Error: Configuration file not found: {args.config}")
    if not args.r1.exists():
        sys.exit(f"Error: R1 file not found: {args.r1}")
    if not args.r2.exists():
        sys.exit(f"Error: R2 file not found: {args.r2}")
    
    # Create output directory if it doesn't exist
    args.out.mkdir(parents=True, exist_ok=True)
    
    # Load configuration and create parser
    try:
        cfg = yaml.safe_load(args.config.read_text())
        read_parser = ReadParser(cfg, cfg.get("params", {}))
    except Exception as e:
        sys.exit(f"Error loading configuration: {e}")
    
    # Set up I/O
    try:
        reader = FastqReader(args.r1, args.r2)
        writer = FastqWriter(str(args.out), args.prefix)
    except Exception as e:
        sys.exit(f"Error setting up I/O: {e}")
    
    # Process reads
    print(f"Processing reads from {args.r1.name} and {args.r2.name}")
    print(f"Output directory: {args.out}")
    print(f"Configuration: {args.config.name}")
    
    try:
        with reader, writer:
            read_clip_and_write(reader, read_parser, writer)
        
        # Show summary statistics
        stats = read_parser.get_parse_log()
        print(f"\nProcessing complete!")
        print(f"Total reads processed: {stats.get('total_reads', 0):,}")
        print(f"Successful reads: {stats.get('successful_reads', 0):,}")
        print(f"Failed reads: {stats.get('failed_reads', 0):,}")
        if stats.get('total_reads', 0) > 0:
            success_rate = stats.get('successful_reads', 0) / stats.get('total_reads', 1) * 100
            print(f"Success rate: {success_rate:.1f}%")
            
    except Exception as e:
        sys.exit(f"Error during processing: {e}")


if __name__ == "__main__":
    main()