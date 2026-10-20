process DownloadGenbank {
    label 'prcoess_single'
    container 'fedora:40'

    input:
    val download_url
    val output_filename

    output:
    path "${output_filename}", emit: genbank_file

    script:
    """
    curl -4 -fL "${download_url}" -o "${output_filename}"
    """
}