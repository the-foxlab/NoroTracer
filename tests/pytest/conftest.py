#!/usr/bin/env python3
"""Shared pytest fixtures for the bin scripts."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
RESOURCES_DIR = Path(__file__).resolve().parent / "resources"


def load_script_module(module_name: str, script_name: str):
	"""Load a script from the bin directory as a Python module."""
	module_path = BIN_DIR / script_name
	spec = importlib.util.spec_from_file_location(module_name, module_path)
	if spec is None or spec.loader is None:
		raise ImportError(f"Could not create a spec for {module_path}")

	module = importlib.util.module_from_spec(spec)
	sys.modules[module_name] = module
	spec.loader.exec_module(module)
	return module


@pytest.fixture
def consensus_fasta() -> tuple[Path, list[SeqRecord]]:
	"""Load the checked-in small FASTA fixture used by the tests."""
	fasta_path = RESOURCES_DIR / "consensus.fasta"
	if not fasta_path.exists():
		raise FileNotFoundError(f"Missing test fixture: {fasta_path}")

	with fasta_path.open("r", encoding="utf-8") as handle:
		records = list(SeqIO.parse(handle, "fasta"))
	return fasta_path, records


@pytest.fixture
def genbank_source() -> tuple[Path, list[SeqRecord]]:
	"""Load the checked-in small GenBank fixture used by the tests."""
	genbank_path = RESOURCES_DIR / "noro_gb_records.gb"
	if not genbank_path.exists():
		raise FileNotFoundError(f"Missing test fixture: {genbank_path}")

	with genbank_path.open("r", encoding="utf-8") as handle:
		records = list(SeqIO.parse(handle, "genbank"))
	return genbank_path, records


@pytest.fixture
def similarity_csv() -> tuple[Path, pd.DataFrame]:
	"""Load the checked-in similarity fixture used by the tests."""
	csv_path = RESOURCES_DIR / "similarities.csv"
	if not csv_path.exists():
		raise FileNotFoundError(f"Missing test fixture: {csv_path}")

	df = pd.read_csv(csv_path)
	return csv_path, df
