

process SPLIT_CLOSEST_ORFS {
    tag "$meta.id"
    label 'process_low'

    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'

    input:
    tuple val(meta), path(consensus_fasta), path(similarity_csv)
    path genbank_file

    output:
    tuple val(meta), path('*_orf1.fasta', arity: '1'), emit: orf1s
    tuple val(meta), path('*_orf2.fasta', arity: '1'), emit: orf2s
    tuple val("${task.process}"), val('split_orfs.py'), eval("split_orfs.py --version"), topic: versions, emit: versions_splitclosestorf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    split_orfs.py \\
        ${genbank_file} \\
        ${similarity_csv} \\
        ${consensus_fasta} \\
        ${args} \\
        --output-prefix split_orfs
    """

    stub:
    """
    touch split_orfs_orf1.fasta
    touch split_orfs_orf2.fasta
    """
}
