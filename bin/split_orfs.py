#!/usr/bin/env python3
"""Extract norovirus ORF1 and ORF2 using the closest GenBank reference.

For every record in the input FASTA files, the script finds the most similar
reference from a multi-record GenBank file and uses gb2seq to cut ORF1 and
ORF2. The reference comparisons can be distributed across multiple processes.
"""

from __future__ import annotations

import argparse
import os
import sys
from itertools import product
from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dark.fasta import FastaReads
from gb2seq.alignment import Gb2Alignment, ReferenceInsertionError
from gb2seq.features import Features


ORF1_NAMES = {
    "ORF-1",
    "ORF1",
    "ORF1 polyprotein",
    "non structural polyprotein",
    "non structural protein",
    "non-structrual polyprotein",
    "non-structural polyporotein precursor",
    "non-structural polyprotein",
    "non-structural polyprotein precursor",
    "non-structural protein",
    "non-structure protein",
    "non-strustural polyprotein",
    "nonstructrual polyprotein",
    "nonstructural polyprotein",
    "nonstructural protein",
}

ORF2_NAMES = {
    "ORF-2",
    "ORF2",
    "VP1",
    "VP1 capsid protein",
    "VP1 protein",
    "capid protein VP1",
    "capsid VP1",
    "capsid protein VP1",
    "major capsid protein",
    "major capsid protein VP1",
    "major capsid region",
    "major structural protein",
    "major viral capsid protein",
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Choose the closest GenBank reference for each consensus sequence "
            "and extract ORF1 and ORF2 with gb2seq."
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
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory (default: current directory).",
    )
    parser.add_argument(
        "-j",
        "--cores",
        type=int,
        default=os.cpu_count() or 1,
        help="Maximum number of worker processes (default: all available cores).",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.cores < 1:
        raise ValueError("--cores must be at least 1")
    for path in (args.genbank, *args.fastas):
        if not path.is_file():
            raise FileNotFoundError(f"Input file does not exist: {path}")


def split_genbank(genbank: Path, directory: Path) -> list[Path]:
    """Write references containing recognized ORF1 and ORF2 products."""
    directory.mkdir(parents=True, exist_ok=True)
    output_paths: list[Path] = []
    record_ids: set[str] = set()
    total = 0

    for record in SeqIO.parse(genbank, "genbank"):
        total += 1
        if record.id in record_ids:
            raise ValueError(f"Duplicate GenBank record ID: {record.id!r}")
        record_ids.add(record.id)
        feature_names = {
            product_name
            for feature in record.features
            for product_name in feature.qualifiers.get("product", [])
        }
        if feature_names.intersection(ORF1_NAMES) and feature_names.intersection(ORF2_NAMES):
            output_path = directory / f"{record.id}.gb"
            SeqIO.write(record, output_path, "genbank")
            output_paths.append(output_path)

    print(f"Wrote {len(output_paths)} of {total} GenBank records.")
    if not output_paths:
        raise ValueError(
            "No GenBank record contains recognized product names for both ORF1 and ORF2."
        )
    return output_paths


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


def compare_to_reference(paths: tuple[Path, Path]) -> tuple[str, Path, float]:
    """Return consensus ID, reference path, and aligned nucleotide identity."""
    genbank_path, fasta_path = paths
    features = Features(genbank_path)
    reads = list(FastaReads(fasta_path))
    if len(reads) != 1:
        raise ValueError(f"Expected one FASTA record in {fasta_path}, found {len(reads)}")

    read = reads[0]
    alignment = Gb2Alignment(read, features, aligner="edlib")
    reference_sequence = alignment.referenceAligned.sequence
    genome_sequence = alignment.genomeAligned.sequence
    identity = sum(
        reference_nt == genome_nt
        for reference_nt, genome_nt in zip(reference_sequence, genome_sequence)
    ) / len(reference_sequence)
    return read.id, genbank_path, identity


def run_similarity_comparisons(
    genbank_paths: Sequence[Path], fasta_paths: Sequence[Path], cores: int
) -> Iterator[tuple[str, Path, float]]:
    """Yield all pairwise similarity results from compare_to_reference."""
    number_of_jobs = len(genbank_paths) * len(fasta_paths)
    if number_of_jobs == 0:
        return

    worker_count = min(cores, number_of_jobs)
    chunksize = max(1, number_of_jobs // (worker_count * 4))
    comparisons = product(genbank_paths, fasta_paths)

    print(
        f"Comparing {len(fasta_paths)} consensus sequences with "
        f"{len(genbank_paths)} references ({number_of_jobs} comparisons) "
        f"using {worker_count} core(s)."
    )
    with Pool(processes=worker_count) as pool:
        for processed, result in enumerate(
            pool.imap_unordered(compare_to_reference, comparisons, chunksize=chunksize),
            start=1,
        ):
            if processed % 10_000 == 0 or processed == number_of_jobs:
                print(
                    f"Processed {processed}/{number_of_jobs} "
                    f"({processed / number_of_jobs:.2%})."
                )
            yield result


def find_closest_references(
    genbank_paths: Sequence[Path], fasta_paths: Sequence[Path], cores: int
) -> dict[str, tuple[Path, float]]:
    closest: dict[str, tuple[Path, float]] = {}

    for consensus_id, genbank_path, identity in run_similarity_comparisons(
        genbank_paths, fasta_paths, cores
    ):
        if consensus_id not in closest or identity > closest[consensus_id][1]:
            closest[consensus_id] = (genbank_path, identity)

    return closest


def matching_feature_name(features: Features, aliases: set[str], orf: str) -> str:
    matches = sorted(aliases & set(features.keys()))
    if not matches:
        raise ValueError(
            f"Closest reference {features.reference.id} has no recognized {orf} feature."
        )
    return matches[0]


def yield_orf_sequences(
    fasta_file: Path,
    orf: str,
    closest: dict[str, tuple[Path, float]],
) -> Iterator[SeqRecord]:
    aliases = ORF1_NAMES if orf == "ORF1" else ORF2_NAMES

    for read in FastaReads(fasta_file):
        genbank_path, _identity = closest[read.id]
        features = Features(genbank_path)
        feature_name = matching_feature_name(features, aliases, orf)
        alignment = Gb2Alignment(read, features)
        try:
            _, cut_orf = alignment.ntSequences(feature_name)
        except ReferenceInsertionError:
            print(
                f"Reference insertion in {features.reference.id} for {read.id}; "
                "retrying with reference gaps allowed.",
                file=sys.stderr,
            )
            _, cut_orf = alignment.ntSequences(
                feature_name, raiseOnReferenceGaps=False
            )

        yield SeqRecord(
            seq=Seq(cut_orf.sequence),
            id=f"{cut_orf.id.split()[0]}_{orf}",
            description=f"{orf} of {cut_orf.id}",
        )


def output_path(fasta: Path, output_dir: Path, orf: str) -> Path:
    return output_dir / f"{fasta.stem}_{orf.lower()}_gb2seqed.fasta"


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        genbank_paths = split_genbank(
            args.genbank, args.output_dir / "genbank_files"
        )
        consensus_paths, _record_ids = split_fastas(
            args.fastas, args.output_dir / "consensi_dir"
        )
        closest = find_closest_references(genbank_paths, consensus_paths, args.cores)

        for fasta in args.fastas:
            for orf in ("ORF1", "ORF2"):
                destination = output_path(fasta, args.output_dir, orf)
                count = SeqIO.write(
                    yield_orf_sequences(fasta, orf, closest), destination, "fasta"
                )
                print(f"Wrote {count} {orf} sequences to {destination}.")
    except (FileNotFoundError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())