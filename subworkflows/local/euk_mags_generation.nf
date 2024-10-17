/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BUSCO                                      } from '../../modules/local/busco/main'
include { EUKCC as EUKCC                             } from '../../modules/local/eukcc/main'
include { ALIGNMENT_LINKTABLE as LINKTABLE           } from '../../modules/local/align_linktable/main'
include { DREP                                       } from '../../modules/local/drep/main'
include { DREP as DREP_MAGS                          } from '../../modules/local/drep/main'
include { BUSCO_EUKCC_QC                             } from '../../modules/local/qc/main'
include { BAT                                        } from '../../modules/local/cat/bat/main'
include { BAT_TAXONOMY_WRITER                        } from '../../modules/local/bat_taxonomy_writer/main'
include { COVERAGE_RECYCLER as COVERAGE_RECYCLER_EUK } from '../../modules/local/coverage_recycler/main'
include { GZIP as GZIP_MAGS                          } from '../../modules/local/utils'
include { GZIP as GZIP_BINS                          } from '../../modules/local/utils'
include { CONCATENATE_QUALITY_FILES                  } from '../../modules/local/euk_utils'
include { MODIFY_QUALITY_FILE                        } from '../../modules/local/euk_utils'
include { FILTER_QUALITY                             } from '../../modules/local/euk_utils'


workflow EUK_MAGS_GENERATION {
    take:
    assemblies_reads_bins  // tuple( meta, assembly_file, [raw_reads], bins, depths, binner_name )
    eukcc_db
    busco_db
    cat_db_folder
    cat_diamond_db
    cat_taxonomy_db

    main:

    ch_versions = Channel.empty()
    ch_log      = Channel.empty()

    /* split the inputs */
    assemblies_reads_bins.multiMap { meta, assembly, reads, bin_folder, depths, binner_name ->
        assembly: [ meta, assembly ]
        reads: [ meta, reads ]
        bin_folder: [ meta, bin_folder ]
        depths: [ meta, depths ]
        binner_name: channel.value(binner_name)
    }.set {
        input
    }

    LINKTABLE( input.reads.join(input.assembly).join(input.bin_folder), input.binner_name) // output: tuple(meta, links.csv)

    EUKCC ( input.binner_name, LINKTABLE_CONCOCT.out.links_table.join ( input.bin_folder ), eukcc_db )

    ch_versions = ch_versions.mix( LINKTABLE.out.versions.first() )
    ch_versions = ch_versions.mix( EUKCC.out.versions.first() )


    

    // -- prepare quality file
    combine_quality = EUKCC.out.eukcc_csv.join( EUKCC_METABAT.out.eukcc_csv )

    // "genome,completeness,contamination" //
    functionCATCSV = { item ->
        def meta = item[0]
        def list_files = [item[1], item[2]]
        return tuple(meta, list_files)
    }

    CONCATENATE_QUALITY_FILES( combine_quality.map( functionCATCSV ), "quality_eukcc.csv" )

    quality = CONCATENATE_QUALITY_FILES.out.concatenated_result

    // -- qs50 -- //
    // [meta, concoct_bins, metabat_bins, merged_concoct, merged_metabat]
    // combine concoct, metabat bins with merged bins (if any)
    collect_data = quality.join( input.bins_concoct ) \
        .join( input.bins_metabat ) \
        .join( EUKCC_CONCOCT.out.eukcc_merged_bins ) \
        .join( EUKCC_METABAT.out.eukcc_merged_bins )

    FILTER_QUALITY( collect_data )

    // input: tuple (meta, genomes/*, quality_file)
    DREP( FILTER_QUALITY.out.qs50_filtered_genomes, params.euk_drep_args, "eukaryotes" )

    ch_versions = ch_versions.mix( DREP.out.versions.first() )

    // -- aggregate by samples
    quality_all_csv = quality.map { meta, quality_file -> quality_file }.collectFile(name: "all.csv", newLine: false)

    MODIFY_QUALITY_FILE( quality_all_csv, "aggregated_euk_quality.csv")

    aggregated_quality = MODIFY_QUALITY_FILE.out.modified_result.map { modified_csv ->
        return tuple([id: "aggregated"], modified_csv)
    }

    // -- drep MAGs --//

    combine_drep = DREP.out.dereplicated_genomes.map{ meta, drep_genomes -> drep_genomes } \
        .flatten() \
        .collect() \
        .map{ agg_genomes ->
            return tuple([id: "aggregated"], agg_genomes)
        }

    DREP_MAGS( combine_drep.join( aggregated_quality ), params.euk_drep_args_mags, 'eukaryotes' )

    ch_versions = ch_versions.mix( DREP_MAGS.out.versions.first() )

    // ---- coverage generation ----- //

    euks_depth = input.metabat_depths.collectFile(name: "euks_depth.txt.gz")

    COVERAGE_RECYCLER_EUK(
        DREP_MAGS.out.dereplicated_genomes,
        euks_depth
    )
    ch_versions = ch_versions.mix( COVERAGE_RECYCLER_EUK.out.versions.first() )

    // ---- QC generation----- //
    // -- BUSCO MAG --//
    drep_result = DREP_MAGS.out.dereplicated_genomes.map { meta, drep_genomes -> drep_genomes }.flatten()

    BUSCO( drep_result, busco_db )

    ch_versions = ch_versions.mix( BUSCO.out.versions.first() )

    BUSCO_EUKCC_QC( 
        aggregated_quality.map { meta, agg_quality_file -> agg_quality_file },
        BUSCO.out.busco_summary.collect(), 
        DREP_MAGS.out.dereplicated_genomes_list.map { meta, drep_genomes -> drep_genomes }
    )

    ch_versions = ch_versions.mix( BUSCO_EUKCC_QC.out.versions.first() )

    // ---- Taxonomy generation ----- //
    // -- BAT --//
    BAT( drep_result, cat_db_folder, cat_taxonomy_db )

    BAT_TAXONOMY_WRITER( BAT.out.bat_names.collect() )

    ch_versions = ch_versions.mix( BAT.out.versions.first() )
    ch_versions = ch_versions.mix( BAT_TAXONOMY_WRITER.out.versions.first() )

    // compress euk genomes
    GZIP_MAGS(drep_result)
    compressed_genomes = GZIP_MAGS.out.compressed.subscribe({ cluster_fasta ->
        cluster_fasta.copyTo("${params.outdir}/genomes_drep/eukaryotes/genomes/${cluster_fasta.name}")
    })
    // compress euk bins
    bins_to_compress = FILTER_QUALITY.out.qs50_filtered_genomes.map{ meta, genomes, quality -> genomes }
    GZIP_BINS(bins_to_compress.flatten())
    compressed_bins = GZIP_BINS.out.compressed.subscribe({ cluster_fasta ->
        cluster_fasta.copyTo("${params.outdir}/bins/eukaryotes/${cluster_fasta.name.split('_')[0]}/${cluster_fasta.name}")
    })

    ch_log = ch_log.mix (EUKCC_METABAT.out.progress_log )
    ch_log = ch_log.mix (EUKCC_CONCOCT.out.progress_log )
    ch_log = ch_log.mix( FILTER_QUALITY.out.progress_log )
    ch_log = ch_log.mix( DREP.out.progress_log )
    ch_log = ch_log.mix( DREP_MAGS.out.progress_log )

    emit:
    genomes                   = GZIP_MAGS.out.compressed.collect()
    stats                     = BUSCO_EUKCC_QC.out.eukcc_final_qc
    coverage                  = COVERAGE_RECYCLER_EUK.out.mag_coverage.map{ meta, coverage_file -> coverage_file }.collect()
    taxonomy                  = BAT_TAXONOMY_WRITER.out.all_bin2classification
    samtools_idxstats_metabat = LINKTABLE_METABAT.out.idxstats
    samtools_idxstats_concoct = LINKTABLE_CONCOCT.out.idxstats
    busco_short_summary       = BUSCO.out.busco_summary
    versions                  = ch_versions
    progress_log              = ch_log
}
