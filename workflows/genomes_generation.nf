include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet; paramsSummaryMap } from 'plugin/nf-validation'

def summary_params = paramsSummaryMap(workflow)

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}

validateParameters()


log.info paramsSummaryLog(workflow)

if (params.help) {
   log.info paramsHelp("nextflow run ebi-metagenomics/genomes-generation --help")
   exit 0
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/mgnify_logo.png")
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~
    DBs
~~~~~~~~~~~~~~~~~~
*/
eukcc_db           = file(params.eukcc_db, checkIfExists: true)
cat_db_folder      = file(params.cat_db_folder, checkIfExists: true)
cat_diamond_db     = file(params.cat_diamond_db, checkIfExists: true)
cat_taxonomy_db    = file(params.cat_taxonomy_db, checkIfExists: true)
gunc_db            = file(params.gunc_db, checkIfExists: true)
checkm2_db         = file(params.checkm2_db, checkIfExists: true)
busco_db           = file(params.busco_db, checkIfExists: true)
gtdbtk_db          = file(params.gtdbtk_db, checkIfExists: true)
rfam_rrna_models   = file(params.rfam_rrna_models, checkIfExists: true)

ref_genome         = file(params.ref_genome)
ref_genome_index   = file("${ref_genome.parent}/*.fa.*")

assembly_software  = file(params.assembly_software_file)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { PREPARE_TSV_FOR_UPLOADER    } from '../modules/local/genome_uploader/main'
include { FINALIZE_LOGGING            } from '../modules/local/utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PROCESS_INPUT        } from '../subworkflows/local/process_input_files'
include { DECONTAMINATION      } from '../subworkflows/local/decontamination'
include { ALIGN                } from '../subworkflows/local/alignment'
include { EUK_MAGS_GENERATION  } from '../subworkflows/local/euk_mags_generation'
include { PROK_MAGS_GENERATION } from '../subworkflows/local/prok_mags_generation'

include { QC_AND_MERGE_READS   } from '../subworkflows/local/qc_and_merge'
include { BINNING              } from '../subworkflows/local/mag_binning'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Run workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion summary
def multiqc_report    = []

workflow GGP {

    ch_versions    = Channel.empty()
    ch_log         = Channel.empty()

    euk_genomes    = Channel.empty()
    stats_euks     = Channel.empty()
    coverage_euks  = Channel.empty()
    prok_genomes   = Channel.empty()
    stats_proks    = Channel.empty()
    coverage_proks = Channel.empty()
    rna            = Channel.empty()
    taxonomy_euks  = Channel.empty()
    taxonomy_proks = Channel.empty()

    // ---- combine data for reads and contigs pre-processing ---- //
    groupReads = { meta, assembly, fq1, fq2 ->
        if (fq2 == []) {
            return tuple(meta, assembly, [fq1])
        }
        else {
            return tuple(meta, assembly, [fq1, fq2])
        }
    }
    assembly_and_runs = Channel.fromSamplesheet("samplesheet", header: true, sep: ',').map(groupReads) // [ meta, assembly_file, [raw_reads] ]

    // ---- pre-processing ---- //
    PROCESS_INPUT( assembly_and_runs ) // output: [ meta, assembly, [raw_reads] ]

    ch_versions = ch_versions.mix( PROCESS_INPUT.out.versions )

    tuple_assemblies = PROCESS_INPUT.out.assembly_and_reads.map{ meta, assembly, _ -> [meta, assembly] }

    // --- trimming reads ---- //
    QC_AND_MERGE_READS( PROCESS_INPUT.out.assembly_and_reads.map { meta, _, reads -> [meta, reads] } )

    ch_versions = ch_versions.mix( QC_AND_MERGE_READS.out.versions )

    // --- decontamination ---- //
    // We need a tuple as the alignment and decontamination module needs the input like that
    DECONTAMINATION( QC_AND_MERGE_READS.out.reads, ref_genome, ref_genome_index )

    ch_versions = ch_versions.mix( DECONTAMINATION.out.versions )

    // --- align reads to assemblies ---- //
    assembly_and_reads = tuple_assemblies.join( DECONTAMINATION.out.decontaminated_reads )

    ALIGN( assembly_and_reads ) // tuple (meta, fasta, [reads])

    ch_versions = ch_versions.mix( ALIGN.out.versions )

    // ---- binning ---- //
    BINNING( ALIGN.out.assembly_bam )
    ch_versions = ch_versions.mix( BINNING.out.versions )

    collectBinsFolder = { meta, bin_folder ->
        [ meta, bin_folder.listFiles().flatten() ]
    }
    concoct_collected_bins = BINNING.out.concoct_bins.map( collectBinsFolder )
    metabat_collected_bins = BINNING.out.metabat_bins.map( collectBinsFolder )
    maxbin_collected_bins = BINNING.out.maxbin_bins.map( collectBinsFolder )

    if ( !params.skip_euk ) {

        // ---- detect euk ---- //
        // input: tuple( meta, assembly_file, [raw_reads], concoct_folder, metabat_folder ), dbs...
        euk_input = tuple_assemblies.join(
            DECONTAMINATION.out.decontaminated_reads
        ).join(
            concoct_collected_bins
        ).join(
            metabat_collected_bins
        )

        EUK_MAGS_GENERATION( 
            euk_input,
            eukcc_db,
            busco_db,
            cat_db_folder,
            cat_diamond_db,
            cat_taxonomy_db
        )

        euk_genomes = euk_genomes.mix( EUK_MAGS_GENERATION.out.genomes )
        stats_euks = stats_euks.mix( EUK_MAGS_GENERATION.out.stats )
        coverage_euks = coverage_euks.mix( EUK_MAGS_GENERATION.out.coverage )
        taxonomy_euks = taxonomy_euks.mix( EUK_MAGS_GENERATION.out.taxonomy )
        ch_versions = ch_versions.mix( EUK_MAGS_GENERATION.out.versions )
        ch_log = ch_log.mix( EUK_MAGS_GENERATION.out.progress_log )
    }

    if ( !params.skip_prok ) {

        // input: tuple( meta, concoct, metabat, maxbin, depth_file)
        collected_binners_and_depth = concoct_collected_bins.join( maxbin_collected_bins ) \
            .join( metabat_collected_bins ) \
            .join( BINNING.out.metabat2depths )

        PROK_MAGS_GENERATION(
            collected_binners_and_depth,
            cat_db_folder,
            cat_diamond_db,
            cat_taxonomy_db,
            gunc_db,
            checkm2_db,
            gtdbtk_db,
            rfam_rrna_models
        )

        prok_genomes = prok_genomes.mix( PROK_MAGS_GENERATION.out.genomes )
        stats_proks = stats_proks.mix( PROK_MAGS_GENERATION.out.stats )
        coverage_proks = coverage_proks.mix( PROK_MAGS_GENERATION.out.coverage )
        rna = rna.mix( PROK_MAGS_GENERATION.out.rna )
        taxonomy_proks = taxonomy_proks.mix( PROK_MAGS_GENERATION.out.taxonomy )
        ch_versions = ch_versions.mix( PROK_MAGS_GENERATION.out.versions )
        ch_log = ch_log.mix( PROK_MAGS_GENERATION.out.progress_log )
    }

    if ( params.skip_euk && params.skip_prok ) {
        println "You skipped proks and euks. No results for MAGs. Exit."
        exit(1)
    }
    else {
        PREPARE_TSV_FOR_UPLOADER(
            euk_genomes.ifEmpty([]),
            prok_genomes.ifEmpty([]),
            assembly_software,
            stats_euks.ifEmpty([]),
            stats_proks.ifEmpty([]),
            coverage_euks.ifEmpty([]),
            coverage_proks.ifEmpty([]),
            rna.ifEmpty([]),
            taxonomy_euks.ifEmpty([]),
            taxonomy_proks.ifEmpty([])
        )
        ch_versions = ch_versions.mix( PREPARE_TSV_FOR_UPLOADER.out.versions )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    pipeline_logging = ch_log.collectFile(name: 'pipeline_logging.txt')
    FINALIZE_LOGGING(pipeline_logging, "structured_pipeline_logging_by_runs.txt")

    //
    // MODULE: MultiQC
    //


    workflow_summary    = WorkflowGenomesGeneration.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowGenomesGeneration.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix( CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect() )
    // TODO: MultiQC 1.17 doesn't support how we are running fastp
    //       should be fixed in the next release -> https://github.com/ewels/MultiQC/issues/2162
    // ch_multiqc_files = ch_multiqc_files.mix( QC_AND_MERGE_READS.out.mqc.map { map, json -> json }.collect().ifEmpty([]) )

    // MULTIQC(
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    // multiqc_report = MULTIQC.out.report.toList()

}
