process ComputeSimilarities {
    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'
    label 'process_medium'
    cpus 4

    input:
    tuple val(meta), path(consensus_fasta)
    path genbank_file

    output:
    tuple val(meta), path('similarities.csv'), emit: similarities

    script:
    """
    python3 ${projectDir}/bin/compute_similarities.py \
        ${genbank_file} \
        ${consensus_fasta} \
        -o similarities.csv \
        -j ${task.cpus}
    """
}