#!/usr/bin/env python3
"""Tests for extracting annotated ORFs from GenBank records."""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO

from tests.conftest import RESOURCES_DIR, load_script_module

try:
    split_orfs = load_script_module("split_orfs", "split_orfs.py")
    extract_genbank_orfs = load_script_module(
        "extract_genbank_orfs", "extract_genbank_orfs.py"
    )
except ModuleNotFoundError:
    pytest.skip("optional runtime dependencies are not installed", allow_module_level=True)


def test_extract_genbank_orfs_writes_annotated_cds_sequences(tmp_path: Path) -> None:
    genbank_file = RESOURCES_DIR / "noro_gb_records.gb"
    records = list(SeqIO.parse(genbank_file, "genbank"))
    output_prefix = tmp_path / "references"

    orf1_path, orf2_path = extract_genbank_orfs.extract_genbank_orfs(
        genbank_file, output_prefix
    )

    assert orf1_path == split_orfs.output_path(output_prefix, "ORF1")
    assert orf2_path == split_orfs.output_path(output_prefix, "ORF2")

    for orf, aliases, output_file in (
        ("ORF1", split_orfs.ORF1_NAMES, orf1_path),
        ("ORF2", split_orfs.ORF2_NAMES, orf2_path),
    ):
        normalized_aliases = {alias.casefold() for alias in aliases}
        expected = []
        for record in records:
            for feature in record.features:
                qualifier_values = {
                    str(value).strip().casefold()
                    for qualifier in ("gene", "product", "note")
                    for value in feature.qualifiers.get(qualifier, [])
                }
                if feature.type == "CDS" and normalized_aliases & qualifier_values:
                    expected.append((f"{record.id}_{orf}", str(feature.extract(record.seq))))

        observed = [
            (record.id, str(record.seq))
            for record in SeqIO.parse(output_file, "fasta")
        ]
        assert observed == expected


def test_cli_accepts_genbank_input_and_output_prefix(tmp_path: Path) -> None:
    genbank_file = RESOURCES_DIR / "noro_gb_records.gb"
    output_prefix = tmp_path / "cli" / "references"

    exit_code = extract_genbank_orfs.main(
        [str(genbank_file), "--output-prefix", str(output_prefix)]
    )

    assert exit_code == 0
    assert split_orfs.output_path(output_prefix, "ORF1").is_file()
    assert split_orfs.output_path(output_prefix, "ORF2").is_file()