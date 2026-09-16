#!/usr/bin/env python3
"""Tests for the similarity computation CLI."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tests.conftest import load_script_module

try:
	compute_similarities = load_script_module("compute_similarities", "compute_similarities.py")
except ModuleNotFoundError:
	pytest.skip("optional runtime dependencies are not installed", allow_module_level=True)


def test_parse_args_accepts_required_inputs(tmp_path: Path) -> None:
	genbank_path = tmp_path / "refs.gb"
	fasta_path = tmp_path / "consensus.fasta"
	genbank_path.write_text("", encoding="utf-8")
	fasta_path.write_text(">consensus_1\nATGC\n", encoding="utf-8")

	args = compute_similarities.parse_args([
		str(genbank_path),
		str(fasta_path),
		"-o",
		"results.csv",
		"-j",
		"2",
	])

	assert args.genbank == genbank_path
	assert args.fastas == [fasta_path]
	assert args.output == Path("results.csv")
	assert args.cores == 2


def test_validate_args_rejects_invalid_core_count(tmp_path: Path) -> None:
	genbank_path = tmp_path / "refs.gb"
	fasta_path = tmp_path / "consensus.fasta"
	genbank_path.write_text("", encoding="utf-8")
	fasta_path.write_text(">consensus_1\nATGC\n", encoding="utf-8")

	args = argparse.Namespace(
		genbank=genbank_path,
		fastas=[fasta_path],
		output=tmp_path / "results.csv",
		cores=0,
		verbose=False,
	)

	with pytest.raises(ValueError, match="--cores must be at least 1"):
		compute_similarities.validate_args(args)


def test_split_fastas_writes_unique_records(consensus_fasta: tuple[Path, list]) -> None:
	fasta_path, records = consensus_fasta
	output_dir = fasta_path.parent / "split_fastas"

	paths, record_ids = compute_similarities.split_fastas([fasta_path], output_dir)

	assert len(paths) == len(records)
	assert len(record_ids) == len(records)
	assert {path.stem for path in paths} == {record.id for record in records}


def test_split_genbank_records_writes_unique_record_files(genbank_source: tuple[Path, list]) -> None:
	genbank_path, records = genbank_source
	output_dir = genbank_path.parent / "split_genbank"

	paths = compute_similarities.split_genbank_records(genbank_path, output_dir)

	assert len(paths) == len(records)
	assert {path.stem for path in paths} == {record.id for record in records}
	assert len(set(path.stem for path in paths)) == len(paths)


def test_compare_to_reference_exact_match_on_real_resource_data(tmp_path: Path) -> None:
	resource_dir = Path(__file__).resolve().parent / "resources"
	genbank_path = resource_dir / "noro_gb_records.gb"
	record = next(SeqIO.parse(genbank_path, "genbank"))
	genbank_single_record_path = tmp_path / f"{record.id}.gb"
	SeqIO.write(record, genbank_single_record_path, "genbank")

	fasta_path = tmp_path / "copied_consensus.fasta"
	SeqIO.write(
		SeqRecord(seq=record.seq, id="consensus_copy", description=""),
		fasta_path,
		"fasta",
	)

	result = compute_similarities.compare_to_reference((genbank_single_record_path, 
													 fasta_path))

	assert result[0] == "consensus_copy"
	assert result[1] == record.id
	assert result[2] == pytest.approx(1.0)


def test_compare_to_reference_matches_expected_identity_after_mutations(tmp_path: Path) -> None:
	resource_dir = Path(__file__).resolve().parent / "resources"
	genbank_path = resource_dir / "noro_gb_records.gb"
	record = next(SeqIO.parse(genbank_path, "genbank"))
	genbank_single_record_path = tmp_path / f"{record.id}.gb"
	SeqIO.write(record, genbank_single_record_path, "genbank")

	sequence = str(record.seq)
	assert len(sequence) >= 101, f"Sequence too short for mutation test: {len(sequence)} nt"
	mutated = list(sequence)
	positions = [0, 10, 100]
	for index in positions:
		mutated[index] = "A" if sequence[index] != "A" else "C"
	mutated_sequence = "".join(mutated)

	fasta_path = tmp_path / "mutated_consensus.fasta"
	SeqIO.write(
		SeqRecord(seq=Seq(mutated_sequence), id="consensus_mutated", description=""),
		fasta_path,
		"fasta",
	)

	result = compute_similarities.compare_to_reference((genbank_single_record_path, fasta_path))
	mismatch_count = sum(
		left != right
		for left, right in zip(sequence, mutated_sequence)
	)
	expected_identity = (len(sequence) - mismatch_count) / len(sequence)

	assert result[0] == "consensus_mutated"
	assert result[1] == record.id
	assert result[2] == pytest.approx(expected_identity)


def test_similarity_output_matches_ground_truth(tmp_path: Path) -> None:
	resource_dir = Path(__file__).resolve().parent / "resources"
	genbank_path = resource_dir / "noro_gb_records.gb"
	fasta_path = resource_dir / "consensus.fasta"
	expected_path = resource_dir / "similarities.csv"
	output_path = tmp_path / "computed_similarities.csv"

	written = compute_similarities.write_similarities_csv(
		genbank=genbank_path,
		fastas=[fasta_path],
		output=output_path,
		cores=1,
	)

	assert written > 0

	actual = pd.read_csv(output_path)
	expected = pd.read_csv(expected_path)
	actual_sorted = actual.sort_values(["consensus_id", "genbank_id", "identity"]).reset_index(drop=True)
	expected_sorted = expected.sort_values(["consensus_id", "genbank_id", "identity"]).reset_index(drop=True)

	pd.testing.assert_frame_equal(actual_sorted, expected_sorted)
