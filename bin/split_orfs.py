#!/usr/bin/env python3
"""Extract norovirus ORF1 and ORF2 using a precomputed closest-reference table.

The script reads consensus sequences from FASTA input, loads closest GenBank
IDs from a similarity CSV, resolves those IDs in a multi-record GenBank source,
and uses gb2seq to cut ORF1 and ORF2.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterator, Sequence, cast

from tempfile import TemporaryDirectory

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dark.fasta import FastaReads
from gb2seq.alignment import Gb2Alignment, ReferenceInsertionError
from gb2seq.features import Features


LOGGER = logging.getLogger("split_orfs")


def configure_logging(verbose: bool) -> None:
    """Configure console logging at INFO or DEBUG level.

    Args:
        verbose: Enable DEBUG-level messages when true.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, stream=sys.stdout, format="%(levelname)s: %(message)s")


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
    """Parse command-line options for ORF extraction.

    Args:
        argv: Arguments to parse, or None to use the process command line.

    Returns:
        The parsed input paths and runtime options.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Extract ORF1 and ORF2 using precomputed closest GenBank IDs "
            "from a similarity table."
        )
    )
    parser.add_argument(
        "genbank_source",
        type=Path,
        help="Multi-record GenBank file containing reference genomes.",
    )
    parser.add_argument(
        "similarity_csv",
        type=Path,
        help=(
            "Similarity CSV with columns consensus_id, genbank_id, identity. "
            "If multiple rows per consensus exist, the highest identity row is used."
        ),
    )
    parser.add_argument(
        "fastas",
        type=Path,
        nargs="+",
        help="One or more FASTA files containing consensus genomes.",
    )
    parser.add_argument(
        "-p",
        "--output-prefix",
        type=Path,
        default=Path("split_orfs"),
        help=(
            "Prefix for output FASTA files. The script writes "
            "<prefix>_orf1.fasta and <prefix>_orf2.fasta."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    """Check that all required input paths point to existing files.

    Args:
        args: Parsed arguments containing the GenBank, similarity, and FASTA paths.

    Raises:
        FileNotFoundError: If any required input path is not a file.
    """
    for path in (args.genbank_source, args.similarity_csv, *args.fastas):
        if not path.is_file():
            raise FileNotFoundError(f"Input file does not exist: {path}")


def build_reference_paths(genbank_source: Path, directory: Path) -> dict[str, Path]:
    """Write each GenBank record to a separate file.

    Args:
        genbank_source: Multi-record GenBank input file.
        directory: Directory in which to write the individual records.

    Returns:
        A mapping from each GenBank record ID to its output path.

    Raises:
        ValueError: If the input contains duplicate IDs or no records.
    """
    directory.mkdir(parents=True, exist_ok=True)
    reference_paths: dict[str, Path] = {}
    record_ids: set[str] = set()

    for record in SeqIO.parse(genbank_source, "genbank"):
        if record.id in record_ids:
            raise ValueError(f"Duplicate GenBank record ID: {record.id!r}")
        record_ids.add(record.id)
        output_path = directory / f"{record.id}.gb"
        SeqIO.write(record, output_path, "genbank")
        reference_paths[record.id] = output_path

    if not reference_paths:
        raise ValueError(f"No GenBank records found in file: {genbank_source}")

    return reference_paths


def load_closest_reference_map(similarity_csv: Path) -> dict[str, str]:
    """Select the highest-identity GenBank reference for each consensus.

    Args:
        similarity_csv: CSV containing consensus IDs, GenBank IDs, and identities.

    Returns:
        A mapping from each consensus ID to its highest-identity GenBank ID.

    Raises:
        ValueError: If required columns are absent or the CSV has no rows.
    """
    table = cast(pd.DataFrame, pd.read_csv(similarity_csv))
    required_columns = {"consensus_id", "genbank_id", "identity"}
    if not required_columns.issubset(table.columns):
        raise ValueError(
            "Similarity CSV must contain consensus_id, genbank_id, and identity columns."
        )

    table["consensus_id"] = table["consensus_id"].astype(str)
    table["genbank_id"] = table["genbank_id"].astype(str)
    table["identity"] = pd.to_numeric(table["identity"], errors="raise")

    max_identity_rows = cast(
        pd.DataFrame,
        table.loc[
            table.groupby("consensus_id")["identity"].idxmax(),
            ["consensus_id", "genbank_id"],
        ],
    )
    closest = {
        consensus_id: genbank_id
        for consensus_id, genbank_id in max_identity_rows.itertuples(index=False, name=None)
    }
    if not closest:
        raise ValueError(f"Similarity CSV has no rows: {similarity_csv}")

    return closest


def matching_feature_name(features: Features, aliases: set[str], orf: str) -> str:
    """Find the first recognized feature name for an ORF.

    Args:
        features: Annotated features from the selected reference record.
        aliases: Accepted names for the requested ORF.
        orf: ORF label used in the error message when no name matches.

    Returns:
        The lexicographically first feature name matching an alias.

    Raises:
        ValueError: If none of the aliases occur in the reference features.
    """
    matches = sorted(aliases & set(features.keys()))
    if not matches:
        raise ValueError(
            f"Closest reference {features.reference.id} has no recognized {orf} feature."
        )
    return matches[0]


def yield_orf_sequences(
    fasta_file: Path,
    orf: str,
    closest: dict[str, str],
    reference_paths: dict[str, Path],
) -> Iterator[SeqRecord]:
    """Yield one aligned ORF sequence record for each consensus in a FASTA file.

    Args:
        fasta_file: FASTA file containing consensus genome sequences.
        orf: ORF label, either ``ORF1`` or ``ORF2``.
        closest: Mapping from consensus IDs to closest GenBank IDs.
        reference_paths: Mapping from GenBank IDs to individual record files.

    Yields:
        SeqRecord objects containing the extracted consensus ORF sequences.

    Raises:
        ValueError: If a consensus or its selected GenBank reference cannot be found.
    """
    aliases = ORF1_NAMES if orf == "ORF1" else ORF2_NAMES

    for read in FastaReads(fasta_file):
        # make sure that that you get the fasta-id and not the description
        consensus_id = read.id.split()[0]
        if consensus_id not in closest:
            raise ValueError(
                f"No closest GenBank ID found for consensus sequence: {consensus_id}"
            )

        reference_id = closest[consensus_id]
        if reference_id not in reference_paths:
            raise ValueError(
                f"GenBank ID {reference_id!r} not found in source GenBank file."
            )

        features = Features(reference_paths[reference_id])
        feature_name = matching_feature_name(features, aliases, orf)
        alignment = Gb2Alignment(read, features)
        try:
            _, cut_orf = alignment.ntSequences(feature_name)
        except ReferenceInsertionError:
            LOGGER.warning(
                f"Reference insertion in {features.reference.id} for {read.id}; "
                "retrying with reference gaps allowed.",
            )
            _, cut_orf = alignment.ntSequences(
                feature_name, raiseOnReferenceGaps=False
            )

        yield SeqRecord(
            seq=Seq(cut_orf.sequence),
            # ensure to get the fasta-id and not the description by splitting on whitespace
            id=f"{cut_orf.id.split()[0]}_{orf}",
            description=f"{orf} of {cut_orf.id}",
        )


def yield_orf_sequences_for_fastas(
    fastas: Sequence[Path],
    orf: str,
    closest: dict[str, str],
    reference_paths: dict[str, Path],
) -> Iterator[SeqRecord]:
    """Yield extracted ORF sequences from multiple FASTA files.

    Args:
        fastas: FASTA files containing consensus genome sequences.
        orf: ORF label to extract from each sequence.
        closest: Mapping from consensus IDs to closest GenBank IDs.
        reference_paths: Mapping from GenBank IDs to individual record files.

    Yields:
        SeqRecord objects from each input file, in file and record order.
    """
    for fasta in fastas:
        yield from yield_orf_sequences(fasta, orf, closest, reference_paths)


def output_path(prefix: Path, orf: str) -> Path:
    """Construct the output FASTA path for an ORF.

    Args:
        prefix: Base path shared by the ORF output files.
        orf: ORF label, whose lowercase form is used in the suffix.

    Returns:
        A path formatted as ``<prefix>_<orf-lowercase>.fasta``.
    """
    return Path(f"{prefix}_{orf.lower()}.fasta")


def main(argv: Sequence[str] | None = None) -> int:
    """Run ORF extraction and write separate ORF1 and ORF2 FASTA files.

    Args:
        argv: Arguments to parse, or None to use the process command line.

    Returns:
        Exit status 0 on success or 2 for invalid or missing input files.
    """
    args = parse_args(argv)
    configure_logging(args.verbose)
    try:
        validate_args(args)
        args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

        with TemporaryDirectory(prefix="split_orfs_") as temp_dir:
            temp_root = Path(temp_dir)
            reference_paths = build_reference_paths(
                args.genbank_source, temp_root / "genbank_files"
            )
            closest = load_closest_reference_map(args.similarity_csv)

            for orf in ("ORF1", "ORF2"):
                destination = output_path(args.output_prefix, orf)
                count = SeqIO.write(
                    yield_orf_sequences_for_fastas(
                        args.fastas, orf, closest, reference_paths
                    ),
                    destination,
                    "fasta",
                )
                LOGGER.info("Wrote %d %s sequences to %s.", count, orf, destination)
    except (FileNotFoundError, ValueError) as error:
        LOGGER.error("%s", error)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())