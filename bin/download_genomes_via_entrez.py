#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
from os import PathLike
from pathlib import Path
from typing import Callable, Iterable, Sequence, TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import download_from_ncbi as dn


LOGGER = logging.getLogger(__name__)

FilePath = str | PathLike[str]


def parse_args() -> argparse.Namespace:
    """Parse command line arguments for the downloader script."""
    parser = argparse.ArgumentParser(
        description="Download norovirus genomes and write genotype FASTA files"
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Email address for NCBI Entrez",
    )
    parser.add_argument(
        "--api-key",
        required=True,
        help="NCBI Entrez API key",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=7450,
        help="Minimum sequence length",
    )
    parser.add_argument(
        "--max-len",
        type=int,
        default=7650,
        help="Maximum sequence length",
    )
    parser.add_argument(
        "--g2-17-tax-id",
        default="552592",
        help="Taxonomy ID for GII.17",
    )
    parser.add_argument(
        "--g2-4-tax-id",
        default="489821",
        help="Taxonomy ID for GII.4",
    )
    parser.add_argument(
        "--g2-tax-id",
        default="122929",
        help="Taxonomy ID for GII",
    )
    parser.add_argument(
        "--blacklist",
        nargs="*",
        default=["KT589391.1"],
        help="Accession IDs to exclude",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory where output FASTA files are written",
    )
    parser.add_argument(
        "--g2-4-output",
        default="g2_4_full_genomes.fasta",
        help="Output FASTA filename for GII.4 records",
    )
    parser.add_argument(
        "--g2-17-output",
        default="g2_17_full_genomes.fasta",
        help="Output FASTA filename for GII.17 records",
    )
    parser.add_argument(
        "--combined-output",
        default="g2_4and17_full_genomes.fasta",
        help="Output FASTA filename for combined GII.4/GII.17 records",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    return parser.parse_args()


def configure_logging(verbose: bool) -> None:
    """Configure logging level and message format."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")


def is_g24(desc: str) -> bool:
    """Return True if description belongs to genotype GII.4 only."""
    return "GII.4" in desc and "GII.17" not in desc


def is_g217(desc: str) -> bool:
    """Return True if description belongs to genotype GII.17 only."""
    return "GII.17" in desc and "GII.4" not in desc


def filter_acc_list(acc_list: Sequence[str], blacklist: set[str]) -> list[str]:
    """Filter accession IDs by blacklist values."""
    return [acc for acc in acc_list if acc not in blacklist]


def write_genotype(
    fn_out: FilePath,
    gen_g2: Iterable[SeqRecord],
    gen_g2_gt: Iterable[SeqRecord],
    is_gt_func: Callable[[str], bool],
    fh_both: TextIO,
) -> tuple[set[str], list[SeqRecord]]:
    """Write genotype FASTA output for one target genotype.

    Args:
        fn_out: Output FASTA path for the current genotype subset.
        gen_g2: Base GII records to scan first.
        gen_g2_gt: Genotype-specific records used to add missing IDs.
        is_gt_func: Predicate that decides whether a record belongs to this genotype.
        fh_both: Open text handle for the combined output FASTA.

    Returns:
        A tuple of:
            - seen_ids: IDs already written for this genotype.
            - left_g2_records: Base GII records that did not match this genotype.
    """
    seen_ids: set[str] = set()
    left_g2_records: list[SeqRecord] = []

    with open(fn_out, "w", encoding="utf-8") as handle:
        for record in gen_g2:
            if is_gt_func(record.description):
                seen_ids.add(record.id)
                SeqIO.write([record], handle, "fasta")
                SeqIO.write([record], fh_both, "fasta")
            else:
                left_g2_records.append(record)

        for record in gen_g2_gt:
            if record.id not in seen_ids:
                SeqIO.write([record], handle, "fasta")
                SeqIO.write([record], fh_both, "fasta")

    return seen_ids, left_g2_records


def count_fasta_records(fn: FilePath) -> int:
    """Count records in a FASTA file."""
    with open(fn, "r", encoding="utf-8") as handle:
        return sum(1 for _ in SeqIO.parse(handle, "fasta"))


def check_for_unique_fasta_ids(fn: FilePath) -> bool:
    """Return True if all FASTA IDs are unique."""
    seen: set[str] = set()
    with open(fn, "r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in seen:
                return False
            seen.add(record.id)
    return True


def run(args: argparse.Namespace) -> None:
    """Run the download and FASTA writing workflow from CLI values."""
    email = args.email
    api_key = args.api_key
    min_len = args.min_len
    max_len = args.max_len

    g2_17_tax_id = args.g2_17_tax_id
    g2_4_tax_id = args.g2_4_tax_id
    g2_tax_id = args.g2_tax_id

    blacklist = set(args.blacklist)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    g24_file = output_dir / str(args.g2_4_output)
    g217_file = output_dir / str(args.g2_17_output)
    combined_file = output_dir / str(args.combined_output)

    LOGGER.info("Starting accession search")
    accs_g2_17 = dn.get_accessions_from_search(
        taxid=g2_17_tax_id,
        min_len=min_len,
        max_len=max_len,
        email=email,
        api_key=api_key,
    )
    LOGGER.info("Found %d GII.17 accessions", len(accs_g2_17))
    accs_g2_4 = dn.get_accessions_from_search(
        taxid=g2_4_tax_id,
        min_len=min_len,
        max_len=max_len,
        email=email,
        api_key=api_key,
    )
    LOGGER.info("Found %d GII.4 accessions", len(accs_g2_4))
    accs_g2 = dn.get_accessions_from_search(
        taxid=g2_tax_id,
        min_len=min_len,
        max_len=max_len,
        email=email,
        api_key=api_key,
    )
    LOGGER.info("Found %d GII accessions", len(accs_g2))

    LOGGER.info(
        "Raw accession counts: G2=%d, G2.4=%d, G2.17=%d",
        len(accs_g2),
        len(accs_g2_4),
        len(accs_g2_17),
    )

    accs_g2 = filter_acc_list(accs_g2, blacklist)
    accs_g2_4 = filter_acc_list(accs_g2_4, blacklist)
    accs_g2_17 = filter_acc_list(accs_g2_17, blacklist)

    LOGGER.info(
        "Filtered accession counts: G2=%d, G2.4=%d, G2.17=%d",
        len(accs_g2),
        len(accs_g2_4),
        len(accs_g2_17),
    )

    gen_g2_17 = dn.fetch_genbank_records(
        accessions=accs_g2_17,
        email=email,
        api_key=api_key,
    )
    gen_g2_4 = dn.fetch_genbank_records(
        accessions=accs_g2_4,
        email=email,
        api_key=api_key,
    )
    gen_g2 = dn.fetch_genbank_records(
        accessions=accs_g2,
        email=email,
        api_key=api_key,
    )

    with combined_file.open("w", encoding="utf-8") as combined_handle:
        _, lefties_g2 = write_genotype(g24_file, gen_g2, gen_g2_4, is_g24, combined_handle)
        _, _ = write_genotype(g217_file, lefties_g2, gen_g2_17, is_g217, combined_handle)

    stats = []
    for out_file in (g24_file, g217_file, combined_file):
        record_count = count_fasta_records(out_file)
        unique_ids = check_for_unique_fasta_ids(out_file)
        stats.append((out_file, record_count, unique_ids))

    for out_file, record_count, unique_ids in stats:
        LOGGER.info(
            "Output stats: file=%s records=%d unique_ids=%s",
            out_file,
            record_count,
            unique_ids,
        )

    assert all(unique_ids for _, _, unique_ids in stats), (
        "Duplicate FASTA IDs found in output files"
    )

    LOGGER.info("Wrote outputs: %s, %s, %s", g24_file, g217_file, combined_file)


def main() -> None:
    """Run the CLI application."""
    args = parse_args()
    configure_logging(args.verbose)
    run(args)


if __name__ == "__main__":
    main()
