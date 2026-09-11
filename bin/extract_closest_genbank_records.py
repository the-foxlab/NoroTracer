#!/usr/bin/env python3
"""Select top-X closest GenBank references per consensus from a similarity CSV.

The script reads consensus IDs from FASTA files, filters a similarity CSV to those
consensus sequences, keeps the top-X references per consensus by similarity, and
writes unique selected GenBank records into one output GenBank file.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence, cast

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


LOGGER = logging.getLogger("extract_closes_genbank_records")


def configure_logging(verbose: bool) -> None:
	"""Configure logging to stdout."""
	level = logging.DEBUG if verbose else logging.INFO
	logging.basicConfig(
		level=level,
		stream=sys.stdout,
		format="%(levelname)s: %(message)s",
	)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
	"""Parse command-line arguments."""
	parser = argparse.ArgumentParser(
		description=(
			"Keep the top-X closest GenBank references per consensus from a "
			"similarity CSV and write selected records to one GenBank file."
		)
	)
	parser.add_argument(
		"fastas",
		type=Path,
		nargs="+",
		help="One or more FASTA files with consensus genomes.",
	)
	parser.add_argument(
		"similarity_csv",
		type=Path,
		help=(
			"Similarity CSV with columns equivalent to ref/cons/sim, or "
			"genbank_path/consensus_id/identity."
		),
	)
	parser.add_argument(
		"topx",
		type=int,
		nargs="?",
		help="Number of closest GenBank references to keep per consensus.",
	)
	parser.add_argument(
		"--topx-range",
		type=str,
		help=(
			"Range start,end for count-only reporting. Example: 1,10 reports "
			"unique reference counts for top1 through top10."
		),
	)
	parser.add_argument(
		"-o",
		"--output-genbank",
		type=Path,
		help="Output GenBank filename containing selected reference records.",
	)
	parser.add_argument(
		"--count-only",
		action="store_true",
		help=(
			"Only report the total number of unique GenBank entries selected by "
			"topx; do not write output file."
		),
	)
	parser.add_argument(
		"--genbank-source",
		type=Path,
		help=(
			"Optional multi-record GenBank source file used when reference values "
			"in the CSV are record IDs instead of file paths."
		),
	)
	parser.add_argument(
		"--verbose",
		action="store_true",
		help="Enable debug logging.",
	)
	return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
	"""Validate user input arguments."""
	if args.topx is None and args.topx_range is None:
		raise ValueError("Provide either topx or --topx-range")

	if args.topx is not None and args.topx < 1:
		raise ValueError("topx must be at least 1")

	if args.topx is not None and args.topx_range is not None:
		raise ValueError("Use either topx or --topx-range, not both")

	for fasta in args.fastas:
		if not fasta.is_file():
			raise FileNotFoundError(f"Input file does not exist: {fasta}")

	if not args.similarity_csv.is_file():
		raise FileNotFoundError(f"Input file does not exist: {args.similarity_csv}")

	if not args.count_only and args.output_genbank is None:
		raise ValueError("--output-genbank is required unless --count-only is used")

	if args.topx_range is not None:
		if not args.count_only:
			raise ValueError("--topx-range is only allowed with --count-only")
		if args.output_genbank is not None:
			raise ValueError("--topx-range cannot be used with --output-genbank")

	if args.genbank_source is not None and not args.genbank_source.is_file():
		raise FileNotFoundError(f"Input file does not exist: {args.genbank_source}")


def read_consensus_ids(fastas: Iterable[Path]) -> set[str]:
	"""Read all consensus IDs from FASTA input files."""
	consensus_ids: set[str] = set()
	for fasta in fastas:
		count = 0
		for record in SeqIO.parse(fasta, "fasta"):
			consensus_ids.add(record.id)
			count += 1
		if count == 0:
			raise ValueError(f"FASTA file contains no records: {fasta}")

	if not consensus_ids:
		raise ValueError("No consensus IDs found in FASTA inputs.")

	return consensus_ids


def load_similarity_table(path: Path) -> pd.DataFrame:
	"""Load and normalize similarity CSV columns to ref/cons/sim."""
	df = cast(pd.DataFrame, pd.read_csv(path))
	columns = set(df.columns)

	if {"ref", "cons", "sim"}.issubset(columns):
		normalized = cast(pd.DataFrame, df[["ref", "cons", "sim"]].copy())
	elif {"genbank_path", "consensus_id", "identity"}.issubset(columns):
		normalized = cast(
			pd.DataFrame,
			df[["genbank_path", "consensus_id", "identity"]].copy(),
		)
		normalized.columns = ["ref", "cons", "sim"]
	else:
		raise ValueError(
			"Unsupported CSV columns. Expected either ref/cons/sim or "
			"genbank_path/consensus_id/identity."
		)

	normalized["ref"] = normalized["ref"].astype(str)
	normalized["cons"] = normalized["cons"].astype(str)
	normalized["sim"] = pd.to_numeric(normalized["sim"], errors="raise")
	return normalized


def select_top_references(
	similarity_df: pd.DataFrame, consensus_ids: set[str], topx: int
) -> set[str]:
	"""Return unique selected references using top-X per consensus logic."""
	dd: dict[str, list[tuple[str, float]]] = defaultdict(list)
	filtered = cast(pd.DataFrame, similarity_df[similarity_df["cons"].isin(consensus_ids)])
	for ref, cons, sim in filtered[["ref", "cons", "sim"]].itertuples(index=False, name=None):
		dd[str(cons)].append((str(ref), float(sim)))

	selected_refs = {
		item[0]
		for values in dd.values()
		for item in sorted(values, key=lambda x: -x[1])[:topx]
	}
	return selected_refs


def parse_topx_range(raw_range: str) -> tuple[int, int]:
	"""Parse and validate a topX range string in the form start,end."""
	parts = [part.strip() for part in raw_range.split(",")]
	if len(parts) != 2:
		raise ValueError("--topx-range must be formatted as start,end (e.g. 1,10)")

	try:
		start = int(parts[0])
		end = int(parts[1])
	except ValueError as error:
		raise ValueError("--topx-range values must be integers") from error

	if start < 1 or end < 1:
		raise ValueError("--topx-range values must be at least 1")
	if start > end:
		raise ValueError("--topx-range start must be less than or equal to end")

	return start, end


def report_topx_range_counts(
	similarity_df: pd.DataFrame,
	consensus_ids: set[str],
	raw_range: str,
) -> None:
	"""Log unique selected reference counts for every topX in a range."""
	start, end = parse_topx_range(raw_range)
	for topx in range(start, end + 1):
		selected_refs = select_top_references(similarity_df, consensus_ids, topx)
		LOGGER.info("topx=%d unique_genbank_entries=%d", topx, len(selected_refs))


def _load_records_from_ref_paths(ref_paths: Iterable[str]) -> list[SeqRecord]:
	"""Load GenBank records from selected reference file paths."""
	records: list[SeqRecord] = []
	for ref_path in sorted(set(ref_paths)):
		path = Path(ref_path)
		if not path.is_file():
			raise FileNotFoundError(
				f"Selected reference path from CSV does not exist: {path}"
			)
		parsed = list(SeqIO.parse(path, "genbank"))
		if len(parsed) != 1:
			raise ValueError(
				f"Expected one GenBank record in {path}, found {len(parsed)}"
			)
		records.append(parsed[0])
	return records


def _load_records_from_source_ids(ref_ids: Iterable[str], source: Path) -> list[SeqRecord]:
	"""Load GenBank records by record ID from a multi-record source file."""
	records_by_id: dict[str, SeqRecord] = {}
	for record in SeqIO.parse(source, "genbank"):
		if record.id in records_by_id:
			raise ValueError(f"Duplicate GenBank record ID in source file: {record.id!r}")
		records_by_id[record.id] = record

	records: list[SeqRecord] = []
	missing: list[str] = []
	for ref_id in sorted(set(ref_ids)):
		record = records_by_id.get(ref_id)
		if record is None:
			missing.append(ref_id)
		else:
			records.append(record)

	if missing:
		preview = ", ".join(missing[:10])
		raise ValueError(
			"Some selected references were not found in --genbank-source: "
			f"{preview}"
		)

	return records


def write_selected_records(
	selected_refs: set[str], output_genbank: Path, genbank_source: Path | None
) -> int:
	"""Write selected GenBank records to output file and return count."""
	if all(Path(ref).is_file() for ref in selected_refs):
		records = _load_records_from_ref_paths(selected_refs)
	elif genbank_source is not None:
		records = _load_records_from_source_ids(selected_refs, genbank_source)
	else:
		raise ValueError(
			"Selected references are not all readable file paths. Provide "
			"--genbank-source to resolve references by record ID."
		)

	output_genbank.parent.mkdir(parents=True, exist_ok=True)
	written = SeqIO.write(records, output_genbank, "genbank")
	return int(written)


def main(argv: Sequence[str] | None = None) -> int:
	"""CLI entrypoint."""
	args = parse_args(argv)
	configure_logging(args.verbose)
	try:
		validate_args(args)
		consensus_ids = read_consensus_ids(args.fastas)
		similarity_df = load_similarity_table(args.similarity_csv)

		if args.topx_range is not None:
			report_topx_range_counts(similarity_df, consensus_ids, args.topx_range)
			return 0

		assert args.topx is not None
		selected_refs = select_top_references(similarity_df, consensus_ids, args.topx)

		LOGGER.info(
			"topx=%d selected %d unique GenBank entries across %d consensus sequences",
			args.topx,
			len(selected_refs),
			len(consensus_ids),
		)

		if args.count_only:
			return 0

		written = write_selected_records(
			selected_refs=selected_refs,
			output_genbank=args.output_genbank,
			genbank_source=args.genbank_source,
		)
		LOGGER.info("Wrote %d GenBank entries to %s", written, args.output_genbank)
	except (FileNotFoundError, ValueError, pd.errors.ParserError) as error:
		LOGGER.error("%s", error)
		return 2

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
