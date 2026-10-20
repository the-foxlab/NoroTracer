process TreeOrf {
    label 'process_medium'
    container 'community.wave.seqera.io/library/iqtree:3.1.3--95869691c4fe61c2'

    input:
    tuple val(meta), path(aln)

    output:
    tuple val(meta), path('*.treefile'), emit: tree

    script:
    def args    = task.ext.args   ?: ''

    """
    iqtree -s ${aln} -nt ${task.cpus} -pre ${aln.simpleName} ${args}
    """
}