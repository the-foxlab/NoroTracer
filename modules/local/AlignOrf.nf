process AlignOrf {
    label 'process_medium'
    container 'community.wave.seqera.io/library/mafft:7.526--8484e078c0b635aa'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.aln.fasta'), emit: alignment

    script:
    """
    mafft --thread ${task.cpus} ${fasta} > ${fasta.simpleName}.aln.fasta
    """
}