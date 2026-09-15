#!/usr/bin/env python3
"""Compute pairwise consensus/reference similarities and write them to CSV."""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from itertools import product
from multiprocessing import Pool
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterable, Iterator, Sequence

from Bio import SeqIO
from dark.fasta import FastaReads
from gb2seq.alignment import Gb2Alignment
from gb2seq.features import Features


LOGGER = logging.getLogger("compute_similarities")


def configure_logging(verbose: bool) -> None:
	"""Configure logging output."""
	level = logging.DEBUG if verbose else logging.INFO
	logging.basicConfig(level=level, stream=sys.stdout, format="%(levelname)s: %(message)s")


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
	parser.add_argument(
		"--verbose",
		action="store_true",
		help="Enable debug logging.",
	)
	return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
	"""Validate command-line arguments."""
	if args.cores < 1:
		raise ValueError("--cores must be at least 1")

	for path in (args.genbank, *args.fastas):
		if not path.is_file():
			raise FileNotFoundError(f"Input file does not exist: {path}")


def split_genbank_records(genbank: Path, directory: Path) -> list[Path]:
	"""Write each GenBank record to its own file and return all output paths."""
	directory.mkdir(parents=True, exist_ok=True)
	output_paths: list[Path] = []
	record_ids: set[str] = set()

	for record in SeqIO.parse(genbank, "genbank"):
		if record.id in record_ids:
			raise ValueError(f"Duplicate GenBank record ID: {record.id!r}")
		record_ids.add(record.id)
		output_path = directory / f"{record.id}.gb"
		SeqIO.write(record, output_path, "genbank")
		output_paths.append(output_path)

	if not output_paths:
		raise ValueError(f"No GenBank records found in file: {genbank}")

	return output_paths


def split_fastas(fastas: Iterable[Path], directory: Path) -> tuple[list[Path], set[str]]:
	"""Write each consensus record to its own FASTA file."""
	directory.mkdir(parents=True, exist_ok=True)
	output_paths: list[Path] = []
	record_ids: set[str] = set()

	for fasta in fastas:
		count = 0
		for record in SeqIO.parse(fasta, "fasta"):
			if record.id in record_ids:
				raise ValueError(
					f"Duplicate FASTA record ID {record.id!r}; record IDs must be unique."
				)
			record_ids.add(record.id)
			output_path = directory / f"{record.id}.fasta"
			SeqIO.write(record, output_path, "fasta")
			output_paths.append(output_path)
			count += 1
		if count == 0:
			raise ValueError(f"FASTA file contains no records: {fasta}")

	return output_paths, record_ids


def compare_to_reference(paths: tuple[Path, Path]) -> tuple[str, str, float]:
	"""Return consensus ID, reference record ID, and aligned nucleotide identity."""
	genbank_path, fasta_path = paths
	features = Features(genbank_path)
	reads = list(FastaReads(fasta_path))
	if len(reads) != 1:
		raise ValueError(f"Expected one FASTA record in {fasta_path}, found {len(reads)}")

	read = reads[0]
	alignment = Gb2Alignment(read, features, aligner="edlib")
	reference_sequence = alignment.referenceAligned.sequence
	genome_sequence = alignment.genomeAligned.sequence
	assert len(reference_sequence) == len(genome_sequence), "Aligned sequences must be the same length"
	identity = sum(
		reference_nt == genome_nt
		for reference_nt, genome_nt in zip(reference_sequence, genome_sequence)
	) / len(reference_sequence)
	return read.id.split(" ")[0], features.reference.id, identity


def run_similarity_comparisons(
	genbank_paths: Sequence[Path], fasta_paths: Sequence[Path], cores: int
) -> Iterator[tuple[str, str, float]]:
	"""Yield all pairwise similarity results from compare_to_reference."""
	number_of_jobs = len(genbank_paths) * len(fasta_paths)
	if number_of_jobs == 0:
		return

	worker_count = min(cores, number_of_jobs)
	chunksize = max(1, number_of_jobs // (worker_count * 4))
	comparisons = product(genbank_paths, fasta_paths)

	LOGGER.info(
		"Comparing %d consensus sequences with %d references (%d comparisons) using %d core(s).",
		len(fasta_paths),
		len(genbank_paths),
		number_of_jobs,
		worker_count,
	)
	with Pool(processes=worker_count) as pool:
		for processed, result in enumerate(
			pool.imap_unordered(compare_to_reference, comparisons, chunksize=chunksize),
			start=1,
		):
			if processed % 10_000 == 0 or processed == number_of_jobs:
				LOGGER.info(
					"Processed %d/%d (%.2f%%).",
					processed,
				number_of_jobs,
				processed / number_of_jobs * 100,
				)
			yield result


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
			writer.writerow(["consensus_id", "genbank_id", "identity"])
			for consensus_id, genbank_id, identity in rows:
				writer.writerow([consensus_id, genbank_id, f"{identity:.8f}"])
				count += 1

	return count


def main(argv: Sequence[str] | None = None) -> int:
	"""CLI entry point."""
	args = parse_args(argv)
	configure_logging(args.verbose)
	try:
		validate_args(args)
		written = write_similarities_csv(
			genbank=args.genbank,
			fastas=args.fastas,
			output=args.output,
			cores=args.cores,
		)
		LOGGER.info("Wrote %d comparison rows to %s.", written, args.output)
	except (FileNotFoundError, ValueError) as error:
		LOGGER.error("%s", error)
		return 2

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
