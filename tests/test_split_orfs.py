#!/usr/bin/env python3
"""Tests for the ORF extraction script."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.conftest import load_script_module

try:
	split_orfs = load_script_module("split_orfs", "split_orfs.py")
except ModuleNotFoundError:
	pytest.skip("optional runtime dependencies are not installed", allow_module_level=True)


def test_load_closest_reference_map_picks_highest_identity(similarity_csv: tuple[Path, object]) -> None:
	csv_path, df = similarity_csv
	closest = split_orfs.load_closest_reference_map(csv_path)

	expected = (
		df.loc[df.groupby("consensus_id")["identity"].idxmax()]
		.set_index("consensus_id")["genbank_id"]
		.to_dict()
	)
	assert closest == expected


def test_output_path_uses_expected_suffix() -> None:
	path = split_orfs.output_path(Path("results/sample"), "ORF1")
	assert path == Path("results/sample_orf1.fasta")


def test_validate_args_rejects_missing_genbank_input(tmp_path: Path) -> None:
	missing_path = tmp_path / "missing.gb"
	fasta_path = tmp_path / "consensus.fasta"
	fasta_path.write_text(">consensus_1\nATGC\n", encoding="utf-8")
	similarity_path = tmp_path / "similarities.csv"
	similarity_path.write_text("consensus_id,genbank_id,identity\nconsensus_1,ref_1,0.88\n", encoding="utf-8")

	args = type(
		"Args",
		(),
		{
			"genbank_source": missing_path,
			"similarity_csv": similarity_path,
			"fastas": [fasta_path],
			"output_prefix": tmp_path / "prefix",
			"verbose": False,
		},
	)()

	with pytest.raises(FileNotFoundError, match="Input file does not exist"):
		split_orfs.validate_args(args)
