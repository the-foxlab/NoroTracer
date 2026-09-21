#!/usr/bin/env python3
"""Tests for the closest reference extraction script."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from tests.conftest import load_script_module

try:
	extract = load_script_module("extract_closest_genbank_records", "extract_closest_genbank_records.py")
except ModuleNotFoundError:
	pytest.skip("optional runtime dependencies are not installed", allow_module_level=True)


def test_load_similarity_table_validates_expected_columns(tmp_path: Path) -> None:
	df = pd.DataFrame(
		[
			{"consensus_id": "consensus_1", "genbank_id": "ref_1", "identity": 0.95},
		],
	)
	csv_path = tmp_path / "similarities.csv"
	df.to_csv(csv_path, index=False)

	loaded = extract.load_similarity_table(csv_path)
	assert list(loaded.columns) == ["consensus_id", "genbank_id", "identity"]
	assert loaded.iloc[0]["genbank_id"] == "ref_1"


def test_select_top_references_keeps_highest_identity_per_consensus(similarity_csv: tuple[Path, pd.DataFrame]) -> None:
	_csv_path, df = similarity_csv
	consensus_ids = set(df["consensus_id"])
	selected = extract.select_top_references(df, consensus_ids, 1)

	expected = set(
		df.sort_values(["consensus_id", "identity"], ascending=[True, False])
		.drop_duplicates("consensus_id", keep="first")["genbank_id"]
	)

	assert selected == expected


def test_parse_topx_range_accepts_valid_range() -> None:
	assert extract.parse_topx_range("1,10") == (1, 10)


def test_parse_topx_range_rejects_invalid_range() -> None:
	with pytest.raises(ValueError, match="start must be less than or equal to end"):
		extract.parse_topx_range("10,1")


def test_validate_args_requires_output_when_not_count_only(tmp_path: Path) -> None:
	fasta_path = tmp_path / "consensus.fasta"
	fasta_path.write_text(">consensus_1\nATGC\n", encoding="utf-8")
	similarity_path = tmp_path / "similarities.csv"
	similarity_path.write_text("consensus_id,genbank_id,identity\nconsensus_1,ref_1,0.9\n", encoding="utf-8")
	genbank_path = tmp_path / "references.gb"
	genbank_path.write_text("", encoding="utf-8")

	args = type(
		"Args",
		(),
		{
			"fasta": fasta_path,
			"similarity_csv": similarity_path,
			"topx": 1,
			"topx_range": None,
			"output_genbank": None,
			"count_only": False,
			"genbank_source": genbank_path,
			"verbose": False,
		},
	)()

	with pytest.raises(ValueError, match="--output-genbank is required"):
		extract.validate_args(args)
