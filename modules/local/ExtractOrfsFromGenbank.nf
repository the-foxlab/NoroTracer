process ExtractOrfsFromGenbank {
    label 'process_low'
    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'

    input:
    path genbank_file
    

    output:
    path '*_orf1.fasta', arity: '1', emit: orf1s
    path '*_orf2.fasta', arity: '1', emit: orf2s

    script:
    """
    python3 ${projectDir}/bin/extract_genbank_orfs.py \
        ${genbank_file} \
        --output-prefix reduced_references \
    """
}