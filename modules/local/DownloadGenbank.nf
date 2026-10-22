process DownloadGenbank {
    label 'prcoess_single'
    container 'community.wave.seqera.io/library/curl:8.23.0--19f33ac87ab0ac0d'

    input:
    val download_url

    output:
    path "noroII.gb", emit: genbank_file

    script:
    """
    curl -4 -fL "${download_url}" -o "noroII.gb"
    """
}