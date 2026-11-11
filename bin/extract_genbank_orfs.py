#!/usr/bin/env python3
"""Extract annotated ORF1 and ORF2 CDS sequences from GenBank records."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterator, Sequence

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from split_orfs import ORF1_NAMES, ORF2_NAMES, output_path


LOGGER = logging.getLogger("extract_genbank_orfs")
ORF_QUALIFIERS = ("gene", "product", "note")


def matching_cds_feature(
    record: SeqRecord,
    aliases: set[str],
    orf: str,
) -> SeqFeature | None:
    """Find the unique CDS feature annotated as the requested ORF.

    Args:
        record: GenBank record containing annotated features.
        aliases: Accepted annotation values for the requested ORF.
        orf: ORF label, used in ambiguity errors.

    Returns:
        The matching CDS feature, or None if the record has no match.

    Raises:
        ValueError: If multiple CDS features match the requested ORF.
    """
    normalized_aliases = {alias.casefold() for alias in aliases}
    matches = []

    for feature in record.features:
        if feature.type != "CDS":
            continue

        qualifier_values = {
            str(value).strip().casefold()
            for qualifier in ORF_QUALIFIERS
            for value in feature.qualifiers.get(qualifier, [])
        }
        if normalized_aliases & qualifier_values:
            matches.append(feature)

    if len(matches) > 1:
        raise ValueError(
            f"GenBank record {record.id!r} has multiple CDS features for {orf}."
        )

    return matches[0] if matches else None


def yield_genbank_orf_sequences(
    genbank_file: Path,
    orf: str,
    aliases: set[str],
) -> Iterator[SeqRecord]:
    """Yield the annotated nucleotide sequence for each record's ORF.

    Args:
        genbank_file: Multi-record GenBank input file.
        orf: ORF label to include in output record IDs.
        aliases: Accepted annotation values for the requested ORF.

    Yields:
        Extracted CDS records in the same order as the GenBank input.

    Raises:
        ValueError: If a GenBank record has multiple matching CDS features,
            or if the input contains no records.
    """
    found_record = False
    for record in SeqIO.parse(genbank_file, "genbank"):
        found_record = True
        feature = matching_cds_feature(record, aliases, orf)
        if feature is None:
            LOGGER.warning("No %s CDS found in GenBank record %s; skipping.", orf, record.id)
            continue

        yield SeqRecord(
            seq=feature.extract(record.seq),
            id=f"{record.id}_{orf}",
            description=f"{orf} from GenBank record {record.id}",
        )

    if not found_record:
        raise ValueError(f"No GenBank records found in file: {genbank_file}")


def extract_genbank_orfs(
    genbank_file: Path,
    output_prefix: Path,
) -> tuple[Path, Path]:
    """Write combined ORF1 and ORF2 nucleotide FASTA files.

    Args:
        genbank_file: Multi-record GenBank input file.
        output_prefix: Prefix for the two output FASTA paths.

    Returns:
        Paths to the ORF1 and ORF2 FASTA files, respectively.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If no GenBank records are present or a record has
            multiple matching CDS features.
    """
    if not genbank_file.is_file():
        raise FileNotFoundError(f"GenBank input file does not exist: {genbank_file}")

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_paths = []
    for orf, aliases in (("ORF1", ORF1_NAMES), ("ORF2", ORF2_NAMES)):
        destination = output_path(output_prefix, orf)
        with destination.open("w", encoding="utf-8") as handle:
            count = SeqIO.write(
                yield_genbank_orf_sequences(genbank_file, orf, aliases),
                handle,
                "fasta",
            )
        LOGGER.info("Wrote %d %s sequences to %s.", count, orf, destination)
        output_paths.append(destination)

    return output_paths[0], output_paths[1]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for GenBank ORF extraction.

    Args:
        argv: Arguments to parse, or None to use the process command line.

    Returns:
        The parsed GenBank input path and output prefix.
    """
    parser = argparse.ArgumentParser(
        description="Extract annotated ORF1 and ORF2 CDS sequences from GenBank records."
    )
    parser.add_argument(
        "genbank_file",
        type=Path,
        help="Multi-record GenBank file containing annotated CDS features.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("genbank_orfs"),
        help="Prefix for the ORF1 and ORF2 FASTA output files.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line GenBank ORF extraction.

    Args:
        argv: Arguments to parse, or None to use the process command line.

    Returns:
        Exit status 0 on success or 2 for invalid input or annotations.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args(argv)
    try:
        extract_genbank_orfs(args.genbank_file, args.output_prefix)
    except (FileNotFoundError, ValueError) as error:
        LOGGER.error("%s", error)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())