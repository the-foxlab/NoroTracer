

process EXTRACT_CLOSEST_GENBANK_RECORDS {
    label 'process_low'

    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'

    input:
    path similarity_csv
    path consensus_fasta
    path genbank_file

    output:
    path 'reduced_references.gb', emit: reduced_genbank
    tuple val("${task.process}"), val('extract_closest_genbank_records.py'), eval("extract_closest_genbank_records.py --version"), topic: versions, emit: versions_extractclosestgenbankrecord

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    extract_closest_genbank_records.py \\
        ${consensus_fasta} \\
        ${similarity_csv} \\
        ${genbank_file} \\
        ${args} \\
        --output-genbank reduced_references.gb
    """

    stub:
    """
    touch reduced_references.gb
    """
}
