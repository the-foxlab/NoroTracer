# Norovirus phylogenetic analysis pipeline

This repository currently contains the implemented workflow for:

1. downloading a GenBank source file
2. computing pairwise consensus/reference similarities
3. selecting the closest GenBank records per consensus
4. splitting the selected reference records into ORF1 and ORF2 sequences for the closest matches

The current implementation does not yet include downstream alignment or tree-building steps. In the current `main.nf`, the workflow stops after `SplitClosestOrfs`.

## Current workflow status

The active Nextflow workflow in `main.nf` runs these stages:

- `DownloadGenbank`
- `ComputeSimilarities`
- `ExtractClosestGenbankRecords`
- `SplitClosestOrfs`

The scripts used by these stages are in `bin/`:

- `compute_similarities.py`
- `extract_closest_genbank_records.py`
- `split_orfs.py`

The downstream MAFFT/IQ-TREE steps described in earlier drafts are not yet implemented in the current pipeline.

## Requirements

- Nextflow
- Python 3
- Docker (for the gb2seq runtime image and containerized execution)
- optionally: conda, if you prefer to run the tests outside Docker

## Running the current workflow

```bash
nextflow run main.nf -params-file params.yaml
```

The current workflow is designed around the runtime parameters in `params.yaml` and uses the local project scripts in `bin/`.

## gb2seq runtime container

The Python steps that depend on `gb2seq` are executed in a custom runtime container built from the upstream `gb2seq` repository. The image is built from a fixed revision and tagged as `gb2seq:runtime`:

```bash
docker build --target runtime -t gb2seq:runtime -f - https://github.com/VirologyCharite/gb2seq.git\#9db97d50185970403bb3265c19ade641f0e9c260 < Dockerfile
```

This container is used in `main.nf` for the Python-based stages that rely on `gb2seq` and Biopython:

- `ComputeSimilarities`
- `ExtractClosestGenbankRecords`
- `SplitClosestOrfs`

A container image is expected to be available at a registry or a local Docker image store under the tag `gb2seq:runtime` when running the pipeline. Placeholder upload location: https://example.com/gb2seq-runtime-container

## Tests

The project includes a small pytest suite under `tests/`.

These tests are intended to validate the Python scripts and the resource-backed fixture data in `tests/resources/`.

Run the tests with the project environment, using the same environment that provides the `gb2seq` runtime as defined by the project Dockerfile:

```bash
conda run -n vevo_snake pytest -q tests
```

If you are using the containerized workflow or a containerized test environment, make sure the `gb2seq:runtime` image is built first and available locally before running the tests that exercise the `gb2seq`-dependent code paths.

## Repository files

- `main.nf`: current workflow definition
- `nextflow.config`: runtime settings and container configuration
- `bin/`: Python scripts for similarity computation, closest-record extraction, and ORF extraction
- `tests/`: pytest suite and small tracked test resources
- `params.yaml`: runtime parameters, including secret or local values
- `params.example.yaml`: example configuration template

## Notes

- This README reflects the current implemented state of the project.
- Any downstream MAFFT/IQ-TREE alignment and phylogenetic tree-generation logic is still planned and not yet part of the active workflow.
- The `gb2seq:runtime` image is required for the currently implemented gb2seq-based stages in the pipeline and validation tests.
