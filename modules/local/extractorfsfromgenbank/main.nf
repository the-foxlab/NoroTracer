

process EXTRACT_ORFS_FROM_GENBANK {
    label 'process_low'

    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'

    input:
    path genbank_file

    output:
    path '*_orf1.fasta', arity: '1', emit: orf1s
    path '*_orf2.fasta', arity: '1', emit: orf2s
    tuple val("${task.process}"), val('extract_genbank_orfs.py'), eval("extract_genbank_orfs.py --version"), topic: versions, emit: versions_extractorfsfromgenbank

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    extract_genbank_orfs.py \\
        ${genbank_file} \\
        ${args} \\
        --output-prefix reduced_references
    """

    stub:
    """
    touch reduced_references_orf1.fasta
    touch reduced_references_orf2.fasta
    """
}
