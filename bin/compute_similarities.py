#!/usr/bin/env python3
"""Compute pairwise consensus/reference similarities and write them to CSV."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence

from split_orfs import (
    run_similarity_comparisons,
    split_fastas,
    split_genbank_records,
    validate_args,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
	"""Parse command-line arguments."""
	parser = argparse.ArgumentParser(
		description=(
			"Compute pairwise similarities between consensus FASTA records and "
			"all entries from a multi-record GenBank file."
		)
	)
	parser.add_argument(
		"genbank",
		type=Path,
		help="Multi-record GenBank file containing candidate references.",
	)
	parser.add_argument(
		"fastas",
		type=Path,
		nargs="+",
		help="One or more FASTA files containing consensus genomes.",
	)
	parser.add_argument(
		"-o",
		"--output",
		type=Path,
		default=Path("similarities.csv"),
		help="Output CSV file (default: similarities.csv).",
	)
	parser.add_argument(
		"-j",
		"--cores",
		type=int,
		default=os.cpu_count() or 1,
		help="Maximum number of worker processes (default: all available cores).",
	)
	return parser.parse_args(argv)


def write_similarities_csv(
	genbank: Path,
	fastas: Sequence[Path],
	output: Path,
	cores: int,
) -> int:
	"""Compute all pairwise similarities and write them to CSV.

	Returns the number of comparison rows written.
	"""
	output.parent.mkdir(parents=True, exist_ok=True)

	with TemporaryDirectory(prefix="compute_similarities_") as temp_dir:
		temp_root = Path(temp_dir)
		genbank_paths = split_genbank_records(genbank, temp_root / "genbank_files")
		fasta_paths, _record_ids = split_fastas(fastas, temp_root / "consensi_dir")
		rows = run_similarity_comparisons(genbank_paths, fasta_paths, cores)

		count = 0
		with output.open("w", encoding="utf-8", newline="") as handle:
			writer = csv.writer(handle)
			writer.writerow(["consensus_id", "genbank_path", "identity"])
			for consensus_id, genbank_path, identity in rows:
				writer.writerow([consensus_id, str(genbank_path), f"{identity:.8f}"])
				count += 1

	return count


def main(argv: Sequence[str] | None = None) -> int:
	"""CLI entry point."""
	args = parse_args(argv)
	try:
		validate_args(args)
		written = write_similarities_csv(
			genbank=args.genbank,
			fastas=args.fastas,
			output=args.output,
			cores=args.cores,
		)
		print(f"Wrote {written} comparison rows to {args.output}.")
	except (FileNotFoundError, ValueError) as error:
		print(f"error: {error}", file=sys.stderr)
		return 2

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
