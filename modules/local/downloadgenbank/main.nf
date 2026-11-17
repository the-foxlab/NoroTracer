

process DOWNLOAD_GENBANK {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/curl:8.23.0--19f33ac87ab0ac0d'


    input:
    val download_url

    output:
    path "${prefix}.gb", emit: genbank_file
    tuple val("${task.process}"), val('curl'), eval('curl --version | sed "1!d; s/^curl //; s/ .*//"'), topic: versions, emit: versions_downloadgenbank

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "noroII"

    """
    curl -fL ${args} "${download_url}" -o "${prefix}.gb"
    """

    stub:
    def args = task.ext.args ?: ''    
    prefix = task.ext.prefix ?: "noroII"

    """
    touch "${prefix}.gb"
    
    """
}
