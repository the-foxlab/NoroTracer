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
include { ExtractOrfsFromGenbank } from '../modules/local/ExtractOrfsFromGenbank'
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


    downloaded_genbank = DownloadGenbank(
        params.genbank_download_url,
    )
    
    ch_sim = ComputeSimilarities(ch_samplesheet, downloaded_genbank.genbank_file)

    ch_combined_csv = ch_sim.similarities
        .map { meta, csv -> csv.text.replace('\r\n', '\n') }
        .collectFile(name: 'all_similarities.csv', keepHeader: true, skip: 1)

    ch_multifasta = ch_samplesheet
        .map { meta, fasta -> fasta }
        .collectFile(name: 'combined.fasta', newLine: true)

    // all channels should just contain one file 
    reduced_genbank = ExtractClosestGenbankRecords(
        ch_combined_csv.first(),
        ch_multifasta.first(),
        downloaded_genbank.genbank_file,
        params.similarity_topx,
    )
    
    SplitClosestOrfs(ch_samplesheet.join(ch_sim), reduced_genbank.reduced_genbank)

    ch_orf1_combined = SplitClosestOrfs.out.orf1s
        .flatMap { meta, files -> files }
        .collectFile(name: 'all_orf1.fasta', newLine: true)

    ch_orf2_combined = SplitClosestOrfs.out.orf2s
        .flatMap { meta, files -> files }
        .collectFile(name: 'all_orf2.fasta',newLine: true)

    ch_orfs = ch_orf1_combined
        .map { fasta -> tuple([id: 'orf1'], fasta) }
        .concat(
            ch_orf2_combined.map { fasta -> tuple([id: 'orf2'], fasta)}
        )
    
    // get ORF1 and ORF2 from the redcued genbank file for alignment and tree building
    // this can already be a fasta file so i don't need the emboss step 
    ch_gb_orfs = ExtractOrfsFromGenbank(reduced_genbank.reduced_genbank)

    ch_orf1_with_ref = ch_gb_orfs.orf1s.concat(ch_orf1_combined)
        .collectFile(name: 'gb_and_cons_orf1.fasta', newLine: true)
        .map { fasta -> tuple([id: 'orf1'], fasta) }

    ch_orf2_with_ref = ch_gb_orfs.orf2s.concat(ch_orf2_combined)
        .collectFile(name: 'gb_and_cons_orf2.fasta', newLine: true)
        .map { fasta -> tuple([id: 'orf2'], fasta) }
        
    ch_orf_with_ref = ch_orf1_with_ref.mix(ch_orf2_with_ref)
    orf_aln = AlignOrf(ch_orf_with_ref)
    tree = TreeOrf(orf_aln.alignment)

    combined_aln_ch = ch_orf_with_ref.mix(orf_aln.alignment)
    combined_cons_orf = SplitClosestOrfs.out.orf1s
        .map{ meta, file -> [meta + [orf: 'orf1'], file] }
        .mix(
            SplitClosestOrfs.out.orf2s
                .map{ meta, file -> [meta + [orf: 'orf2'], file] }
        )
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
    all_similarities_csv = ch_combined_csv
    reduced_genbank_file = reduced_genbank.reduced_genbank
    split_orf_output = combined_cons_orf
    comb_aln_output = combined_aln_ch
    tree_output = tree.tree
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
