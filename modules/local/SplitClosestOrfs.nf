process SplitClosestOrfs {
    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'
    label 'process_low'

    input:
    tuple val(meta), path(consensus_fasta), path(similarity_csv)
    path genbank_file

    output:
    tuple val(meta), path('*_orf1.fasta', arity: '1..*'), emit: orf1s
    tuple val(meta), path('*_orf2.fasta', arity: '1..*'), emit: orf2s

    script:
    """
    python3 ${projectDir}/bin/split_orfs.py \
        ${genbank_file} \
        ${similarity_csv} \
        ${consensus_fasta} \
        --output-prefix split_orfs
    """
}