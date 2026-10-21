/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { DownloadGenbank } from '../modules/local/DownloadGenbank'
include { ComputeSimilarities } from '../modules/local/ComputeSimilarities'
include { ExtractClosestGenbankRecords } from '../modules/local/ExtractClosestGenbankRecords'
include { SplitClosestOrfs } from '../modules/local/SplitClosestOrfs'
include { AlignOrf } from '../modules/local/AlignOrf'
include { TreeOrf } from '../modules/local/TreeOrf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NOROTRACER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    outdir

    main:

    def ch_versions = channel.empty()

     /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    START Temporary hack to create a channel of consensus FASTA files that match the IDs in the multi-FASTA file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    //  def multi_fasta = file(
    //     params.used_fasta_ids,
    //     checkIfExists: true
    // )

    // def amplicon_nf_root_directory = params.amplicon_nf_result_dir

    // // Read the FASTA IDs into a  set
    // def rids = multi_fasta
    //     .splitFasta(record: [id: true])
    //     .collect { record -> record.id }
    //     .toSet()

    // // Emit one [id, file] tuple per matching consensus file
    // ch_consensus = channel
    //     .fromPath("${amplicon_nf_root_directory}/Noro-P*/*.fasta")
    //     .map { fasta ->
    //         def id = fasta.baseName.split(/\./)[0]
    //         tuple(id, fasta)
    //     }
    //     .filter { id, fasta -> id in rids }
    //     .map{id,fasta -> [[id: id], fasta]}

     /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END Temporary hack to create a channel of consensus FASTA files that match the IDs in the multi-FASTA file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    downloaded_genbank = DownloadGenbank(
        params.genbank_download_url,
        params.genbank_file,
    )
    
    genbank_file = downloaded_genbank.genbank_file

    ch_sim = ComputeSimilarities(ch_samplesheet, genbank_file)

    ch_combined_csv = ch_sim.similarities
    .map { meta, csv -> csv.text.replace('\r\n', '\n') }
    .collectFile(
        name: 'all_similarities.csv',
        keepHeader: true,
        skip: 1
    )
    ch_multifasta = ch_samplesheet
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
    
    SplitClosestOrfs(ch_samplesheet.join(ch_sim), reduced_genbank.reduced_genbank)

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


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'norotracer_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    similarities_csv = ComputeSimilarities.out.similarities
    reduced_genbank_file = reduced_genbank.reduced_genbank
    split_orf1_output = SplitClosestOrfs.out.orf1s
    split_orf2_output = SplitClosestOrfs.out.orf2s
    orf1_aln_output = ch_separated.orf1
    orf2_aln_output = ch_separated.orf2
    tree_output = tree.tree
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
