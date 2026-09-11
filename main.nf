nextflow.enable.dsl=2


process ComputeSimilarities {
    label 'python'

    input:
    path genbank_file
    path consensus_fasta

    output:
    path 'similarities.csv', emit: similarities_csv

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

    input:
    path similarity_csv
    path consensus_fasta
    path genbank_file
    val topx

    output:
    path 'reduced_references.gb', emit: reduced_genbank

    script:
    """
    python3 ${projectDir}/bin/extract_closes_genbank_records.py \
        ${consensus_fasta} \
        ${similarity_csv} \
        ${topx} \
        --output-genbank reduced_references.gb \
        --genbank-source ${genbank_file}
    """
}

process SplitClosestOrfs {
    label 'python'

    input:
    path reduced_genbank
    path consensus_fasta

    output:
    path 'outdir/*.fasta', emit: orfs

    script:
    """
    mkdir -p outdir
    python3 ${projectDir}/bin/split_orfs.py \
        ${reduced_genbank} \
        ${consensus_fasta} \
        -o outdir \
        -j ${task.cpus}
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

workflow {

    main:
    genbank_file = Channel.fromPath(params.genbank_file)
    consensus_fasta = Channel.fromPath(params.consensus_file)

    similarities = ComputeSimilarities(genbank_file, consensus_fasta)
    reduced_genbank = ExtractClosestGenbankRecords(
        similarities.similarities_csv,
        consensus_fasta,
        genbank_file,
        params.similarity_topx,
    )
    split_orfs = SplitClosestOrfs(reduced_genbank.reduced_genbank, consensus_fasta)

    publish:
    similarities_csv = similarities.similarities_csv
    reduced_genbank_file = reduced_genbank.reduced_genbank
    split_orfs_output = split_orfs.orfs
}

output {
    similarities_csv {
        path 'similarity_outputs'
        mode 'copy'
    }
    reduced_genbank_file {
        path 'filtered_references'
        mode 'copy'
    }
    split_orfs_output {
        path 'split_orfs'
        mode 'copy'
    }
}