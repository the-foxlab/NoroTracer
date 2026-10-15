process ExtractClosestGenbankRecords {
    label 'process_low'
    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'

    input:
    path similarity_csv
    path consensus_fasta
    path genbank_file
    val topx

    output:
    path 'reduced_references.gb', emit: reduced_genbank

    script:
    """
    python3 ${projectDir}/bin/extract_closest_genbank_records.py \
        ${consensus_fasta} \
        ${similarity_csv} \
        ${topx} \
        --output-genbank reduced_references.gb \
        --genbank-source ${genbank_file}
    """
}