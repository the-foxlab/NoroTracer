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
from typing import Iterable, Iterator, Sequence, cast

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
		"fasta",
		type=Path,
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
		required=True,
		help=(
			"Multi-record GenBank source file used to resolve selected GenBank IDs."
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

	if not args.fasta.is_file():
			raise FileNotFoundError(f"Input file does not exist: {args.fasta}")

	if not args.similarity_csv.is_file():
		raise FileNotFoundError(f"Input file does not exist: {args.similarity_csv}")

	if not args.count_only and args.output_genbank is None:
		raise ValueError("--output-genbank is required unless --count-only is used")

	if args.topx_range is not None:
		if not args.count_only:
			raise ValueError("--topx-range is only allowed with --count-only")
		if args.output_genbank is not None:
			raise ValueError("--topx-range cannot be used with --output-genbank")

	if not args.genbank_source.is_file():
		raise FileNotFoundError(f"Input file does not exist: {args.genbank_source}")



def read_consensus_ids(fasta: Path) -> set[str]:
	"""Read all consensus IDs from FASTA input files."""
	consensus_ids: set[str] = set()

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
	if not set(df.columns) == {"consensus_id", "genbank_id", "identity"}:
		raise ValueError(f"Unexpected CSV columns: {df.columns}")

	normalized = cast(pd.DataFrame, df[["consensus_id", "genbank_id", "identity"]])
	normalized["consensus_id"] = normalized["consensus_id"].astype(str)
	normalized["genbank_id"] = normalized["genbank_id"].astype(str)
	normalized["identity"] = pd.to_numeric(normalized["identity"], errors="raise")
	return normalized


def select_top_references(
	similarity_df: pd.DataFrame, consensus_ids: set[str], topx: int
) -> set[str]:
	"""Return unique selected references using top-X per consensus logic."""
	dd: dict[str, list[tuple[str, float]]] = defaultdict(list)

	LOGGER.debug("Selecting top %d references per consensus", topx)
	LOGGER.debug("Before filtering, similarity_df has %d rows", len(similarity_df))

	filtered = cast(pd.DataFrame, similarity_df[similarity_df["consensus_id"].isin(consensus_ids)])

	LOGGER.debug("After filtering, similarity_df has %d rows", len(filtered))

	for ref, cons, sim in filtered[["genbank_id", "consensus_id", "identity"]].itertuples(index=False, name=None):
		dd[str(cons)].append((str(ref), float(sim)))

	selected_refs = {
		item[0]
		for values in dd.values()
		for item in sorted(values, key=lambda x: -x[1])[:topx]
	}
	LOGGER.debug("Selected %d unique references across %d consensus sequences", len(selected_refs), len(dd))
	
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


def _validate_source_ids(ref_ids: set[str], source: Path) -> None:
	"""Ensure all requested GenBank IDs exist in source before writing."""
	available_ids = {record.id for record in SeqIO.parse(source, "genbank")}
	missing = sorted(ref_ids - available_ids)
	if missing:
		preview = ", ".join(missing[:10])
		raise ValueError(
			"Some selected references were not found in --genbank-source: "
			f"{preview}"
		)


def _load_records_from_source_ids(ref_ids: Iterable[str], source: Path) -> Iterator[SeqRecord]:
	"""Yield GenBank records by record ID from a multi-record source file."""
	selected = set(ref_ids)
	for record in SeqIO.parse(source, "genbank"):
		if record.id in selected:
			yield record


def write_selected_records(
	selected_refs: set[str], output_genbank: Path, genbank_source: Path
) -> int:
	"""Write selected GenBank records to output file and return count."""
	_validate_source_ids(selected_refs, genbank_source)
	records = _load_records_from_source_ids(selected_refs, genbank_source)

	output_genbank.parent.mkdir(parents=True, exist_ok=True)
	written = SeqIO.write(records, output_genbank, "genbank")
	return int(written)


def main(argv: Sequence[str] | None = None) -> int:
	"""CLI entrypoint."""
	args = parse_args(argv)
	configure_logging(args.verbose)
	try:
		validate_args(args)
		consensus_ids = read_consensus_ids(args.fasta)
		similarity_df = load_similarity_table(args.similarity_csv)

		LOGGER.info(
			"Loaded %d similarity rows for %d consensus sequences",
			len(similarity_df),
			len(consensus_ids),
		)
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
