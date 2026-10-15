nextflow.enable.dsl=2

include { DownloadGenbank } from './modules/DownloadGenbank'
include { ComputeSimilarities } from './modules/ComputeSimilarities'
include { ExtractClosestGenbankRecords } from './modules/ExtractClosestGenbankRecords'
include { SplitClosestOrfs } from './modules/SplitClosestOrfs'
include { AlignOrf } from './modules/AlignOrf'
include { TreeOrf } from './modules/TreeOrf'


workflow CREATE_ALIGNMENTS_AND_TREES {

    take:
    ch_consensus

    main:
    downloaded_genbank = DownloadGenbank(
        params.genbank_download_url,
        params.genbank_file,
    )
    
    genbank_file = downloaded_genbank.genbank_file

    ch_sim = ComputeSimilarities(ch_consensus, genbank_file)

    ch_combined_csv = ch_sim.similarities
    .map { meta, csv -> csv }
    .collectFile(
        name: 'all_similarities.csv',
        keepHeader: true,
        skip: 1
    )
    ch_multifasta = ch_consensus
    .map { meta, fasta -> fasta }
    .collectFile(
        name: 'combined.fasta',
        newLine: true
    )

    // all channels should just contain one file 
    reduced_genbank = ExtractClosestGenbankRecords(
        ch_combined_csv.first(),
        ch_multifasta.first(),
        genbank_file,
        params.similarity_topx,
    )
    
    SplitClosestOrfs(ch_consensus.join(ch_sim), reduced_genbank.reduced_genbank)

    ch_orf1_combined = SplitClosestOrfs.out.orf1s
    .flatMap { meta, files -> files }
    .collectFile(
        name: 'all_orf1.fasta',
        newLine: true
    )

    ch_orf2_combined = SplitClosestOrfs.out.orf2s
        .flatMap { meta, files -> files }
        .collectFile(
            name: 'all_orf2.fasta',
            newLine: true
        )

    ch_orfs = ch_orf1_combined
    .map { fasta -> tuple([id: 'orf1'], fasta) }
    .mix(
        ch_orf2_combined.map { fasta ->
            tuple([id: 'orf2'], fasta)
        }
    )

    orf_aln = AlignOrf(ch_orfs)
    ch_separated = orf_aln.alignment.branch { meta, fasta ->
        orf1: meta.id == 'orf1'
        orf2: meta.id == 'orf2'
        }

    tree = TreeOrf(orf_aln.alignment)

    emit:
    similarities_csv = ComputeSimilarities.out.similarities
    reduced_genbank_file = reduced_genbank.reduced_genbank
    split_orf1_output = SplitClosestOrfs.out.orf1s
    split_orf2_output = SplitClosestOrfs.out.orf2s
    orf1_aln_output = ch_separated.orf1
    orf2_aln_output = ch_separated.orf2
    tree_output = tree.tree
}


workflow {

    main:

    // Temporary hack to create a channel of consensus FASTA files that match the IDs in the multi-FASTA file
     def multi_fasta = file(
        params.used_fasta_ids,
        checkIfExists: true
    )

    def amplicon_nf_root_directory = params.amplicon_nf_result_dir

    // Read the FASTA IDs into a  set
    def rids = multi_fasta
        .splitFasta(record: [id: true])
        .collect { record -> record.id }
        .toSet()

    // Emit one [id, file] tuple per matching consensus file
    ch_consensus = channel
        .fromPath("${amplicon_nf_root_directory}/Noro-P*/*.fasta")
        .map { fasta ->
            def id = fasta.baseName.split(/\./)[0]
            tuple(id, fasta)
        }
        .filter { id, fasta -> id in rids }
        .map{id,fasta -> [[id: id], fasta]}

    
    // That's the point to wire together with the amplicon_nf consensus_fasta channel from the workflow AMPLICON_NF 
    CREATE_ALIGNMENTS_AND_TREES(ch_consensus)

    publish:
    similarities_csv = CREATE_ALIGNMENTS_AND_TREES.out.similarities_csv
    reduced_genbank_file = CREATE_ALIGNMENTS_AND_TREES.out.reduced_genbank_file
    split_orf1_output = CREATE_ALIGNMENTS_AND_TREES.out.split_orf1_output
    orf1_aln_output = CREATE_ALIGNMENTS_AND_TREES.out.orf1_aln_output
    orf2_aln_output = CREATE_ALIGNMENTS_AND_TREES.out.orf2_aln_output
    tree_output = CREATE_ALIGNMENTS_AND_TREES.out.tree_output
    split_orf2_output = CREATE_ALIGNMENTS_AND_TREES.out.split_orf2_output

}




output {
    similarities_csv {
        path {meta, csv -> "similarity_outputs/${meta.id}_similarities/"} 
        mode 'copy'
    }
    reduced_genbank_file {
        path 'filtered_references'
        mode 'copy'
    }
    split_orf1_output {
        path {meta, fasta -> "split_orfs/${meta.id}_orf1/"}
        mode 'copy'
    }
    split_orf2_output {
        path {meta, fasta -> "split_orfs/${meta.id}_orf2/"}
        mode 'copy'
    }
    orf1_aln_output {
        path {aln -> "alignments/orf1/"}
        mode 'copy'
    }
    orf2_aln_output {
        path {aln -> "alignments/orf2/"}
        mode 'copy'
    }
    tree_output {
        path {meta, tree -> "trees/${meta.id}/"}
        mode 'copy'
    }
}