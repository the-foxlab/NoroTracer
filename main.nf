nextflow.enable.dsl=2


process DownloadGenbank {
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

process ComputeSimilarities {
    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'
    label 'python'
    cpus 4

    input:
    tuple val(meta), path(consensus_fasta)
    path genbank_file

    output:
    tuple val(meta),path('similarities.csv'), emit: similarities

    script:
    """
    python3 ${projectDir}/bin/compute_similarities.py \
        ${genbank_file} \
        ${consensus_fasta} \
        -o similarities.csv \
        -j ${task.cpus}
    """
}

process ExtractClosestGenbankRecords {
    label 'python'
    // container 'quay.io/biocontainers/biopython:1.70'
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

process SplitClosestOrfs {
    container 'ghcr.io/udogi/gb2seq-env:0.1.0-runtime'
    label 'python'

    input:
    path reduced_genbank
    path similarity_csv
    path consensus_fasta

    output:
    path '*_orf1.fasta', emit: orf1s
    path '*_orf2.fasta', emit: orf2s

    script:
    """
    python3 ${projectDir}/bin/split_orfs.py \
        ${reduced_genbank} \
        ${similarity_csv} \
        ${consensus_fasta} \
        --output-prefix split_orfs
    """
}

process AlignOrf1 {
    label 'mafft'

    input:
    path orf1_fasta

    output:
    path '*.aln.fasta'

    script:
    """
    mafft --thread ${task.cpus} ${orf1_fasta} > ${orf1_fasta.simpleName}.aln.fasta
    """
}

process AlignOrf2 {
    label 'mafft'

    input:
    path orf2_fasta

    output:
    path '*.aln.fasta'

    script:
    """
    mafft --thread ${task.cpus} ${orf2_fasta} > ${orf2_fasta.simpleName}.aln.fasta
    """
}

process AlignGenome {
    label 'mafft'

    input:
    path genome_fasta

    output:
    path '*.aln.fasta'

    script:
    """
    mafft --thread ${task.cpus} ${genome_fasta} > ${genome_fasta.simpleName}.aln.fasta
    """
}

process TreeOrf1 {
    label 'iqtree'

    input:
    path aln

    output:
    path '*.treefile'

    script:
    """
    iqtree2 -s ${aln} -nt ${task.cpus} -pre ${aln.simpleName}
    """
}

process TreeOrf2 {
    label 'iqtree'

    input:
    path aln

    output:
    path '*.treefile'

    script:
    """
    iqtree2 -s ${aln} -nt ${task.cpus} -pre ${aln.simpleName}
    """
}

process TreeGenome {
    label 'iqtree'

    input:
    path aln

    output:
    path '*.treefile'

    script:
    """
    iqtree2 -s ${aln} -nt ${task.cpus} -pre ${aln.simpleName}
    """
}

process CreateMetadata {
    label 'python'

    input:
    path treefiles

    output:
    path params.metadata_file

    script:
    """
    python3 ${projectDir}/script.py metadata --trees *.treefile --out ${params.metadata_file}
    """
}

workflow CREATE_ALIGNMENTS_AND_TREES {

    take:
    ch_consensus

    main:
    downloaded_genbank = DownloadGenbank(
        params.genbank_download_url,
        params.genbank_file,
    )
    
    genbank_file = downloaded_genbank.genbank_file

    ComputeSimilarities(ch_consensus, genbank_file)
    // reduced_genbank = ExtractClosestGenbankRecords(
    //     similarities.similarities_csv,
    //     consensus_fasta,
    //     genbank_file,
    //     params.similarity_topx,
    // )
    // split_orfs = SplitClosestOrfs(
    //     reduced_genbank.reduced_genbank,
    //     similarities.similarities_csv,
    //     consensus_fasta,
    // )

    emit:
    similarities_csv = ComputeSimilarities.out.similarities
    // reduced_genbank_file = reduced_genbank.reduced_genbank
    // split_orf1_output = split_orfs.orf1s
    // split_orf2_output = split_orfs.orf2s
}


workflow {

    main:

    // Temporary hack to create a channel of consensus FASTA files that match the IDs in the multi-FASTA file
     def multi_fasta = file(
        '/home/udo/shared/researchers/udo_gieraths/code/noro_phylogenetic_analysis/g2_17_consensus.fasta',
        checkIfExists: true
    )

    def root_directory = '/home/udo/shared/researchers/udo_gieraths/code/amplicon-nf/results'

    // Read the FASTA IDs into a regular Groovy set
    def rids = multi_fasta
        .splitFasta(record: [id: true])
        .collect { record -> record.id }
        .toSet()

    // Emit one [id, file] tuple per matching consensus file
    ch_consensus = channel
        .fromPath("${root_directory}/Noro-P*/*.fasta")
        .map { fasta ->
            def id = fasta.baseName.split(/\./)[0]
            tuple(id, fasta)
        }
        .filter { id, fasta -> id in rids }
        .map{id,fasta -> [id, fasta]}

    
    // That's later the point to wire the amplicon_nf consensus_fasta channel from the workflow AMPLICON_NF 
    CREATE_ALIGNMENTS_AND_TREES(ch_consensus)

    publish:
    similarities_csv = CREATE_ALIGNMENTS_AND_TREES.out.similarities_csv
    // reduced_genbank_file = CREATE_ALIGNMENTS_AND_TREES.out.reduced_genbank_file
    // split_orf1_output = CREATE_ALIGNMENTS_AND_TREES.out.split_orf1_output
    // split_orf2_output = CREATE_ALIGNMENTS_AND_TREES.out.split_orf2_output

}




output {
    similarities_csv {
        path {id, csv -> "similarity_outputs/${id}_similarities/"} 
        mode 'copy'
    }
    // reduced_genbank_file {
    //     path 'filtered_references'
    //     mode 'copy'
    // }
    // split_orf1_output {
    //     path 'split_orfs'
    //     mode 'copy'
    // }
    // split_orf2_output {
    //     path 'split_orfs'
    //     mode 'copy'
    // }
}