

process COMPUTE_SIMILARITIES {
    tag "$meta.id"
    label 'process_medium'

    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'

    input:
    tuple val(meta), path(consensus_fasta)
    path genbank_file

    output:
    tuple val(meta), path('similarities.csv'), emit: similarities
    tuple val("${task.process}"), val('compute_similarities.py'), eval("compute_similarities.py --version"), topic: versions, emit: versions_computesimilarities

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!workflow.containerEngine) {
        error "COMPUTE_SIMILARITIES requires a container runtime. Enable Docker, Apptainer, or Singularity."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    compute_similarities.py \
        ${genbank_file} \
        ${consensus_fasta} \
        -o similarities.csv \
        -j ${task.cpus}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch similarities.csv
    
    """
}
