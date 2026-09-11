# Simple Norovirus Phylogenetics Pipeline (Nextflow)

This is a minimal and easy-to-read Nextflow pipeline based on your sketch:

1. download records and write fasta
2. split into ORF1 and ORF2
3. align ORF1, ORF2, and full genome (MAFFT)
4. build trees for each alignment (IQ-TREE)
5. create metadata

For Python-based steps, the pipeline calls a dummy `script.py`.

## Files

- `main.nf`: pipeline workflow and processes
- `nextflow.config`: local executor + Docker container settings
- `script.py`: dummy script for download/split/metadata

## Requirements

- Nextflow
- Python 3
- MAFFT
- IQ-TREE
- Optional: Docker (if you want containerized execution)

## Run

```bash
nextflow run main.nf -params-file params.yaml
```

Override default file names via params:

```bash
nextflow run main.nf -params-file params.yaml \
	--records_file input_records.fasta \
	--orf1_file part_orf1.fasta \
	--orf2_file part_orf2.fasta \
	--genome_file whole_genome.fasta \
	--metadata_file summary.tsv
```

To run with containers instead:

```bash
nextflow run main.nf -params-file params.yaml -profile docker
```

Outputs are written under `results/`:

- `results/01_download/records.fasta`
- `results/02_split/` with the configured ORF/genome filenames
- `results/03_align/*.aln.fasta`
- `results/04_trees/*.treefile`
- `results/05_metadata/` with the configured metadata filename

## Notes

- This is intentionally simple and can be extended later.
- Container image tags in `nextflow.config` can be changed if needed for your environment.
