/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     EukCC subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EUKCC as EUKCC                             } from '../../modules/local/eukcc/main'
include { ALIGNMENT_LINKTABLE as LINKTABLE           } from '../../modules/local/align_linktable/main'

workflow EUKCC_MERGE {
    take:
    assemblies_reads_bins  // tuple( meta, assembly_file, [raw_reads], bins, depths, binner_name )
    eukcc_db

    main:

    ch_versions = Channel.empty()
    ch_log      = Channel.empty()

    /* split the inputs */
    assemblies_reads_bins.multiMap { meta, assembly, reads, bin_folder, depths, binner_name ->
        assembly: [ meta, assembly ]
        reads: [ meta, reads ]
        bin_folder: [ meta, bin_folder ]
        depths: [ meta, depths ]
        binner_name: binner_name
    }.set {
        input
    }

    LINKTABLE( input.reads.join(input.assembly).join(input.bin_folder).join(input.depths), input.binner_name) // output: tuple(meta, links.csv)

    EUKCC ( LINKTABLE.out.links_table.join( input.bin_folder ), input.binner_name, eukcc_db )

    ch_versions = ch_versions.mix( LINKTABLE.out.versions.first() )
    ch_versions = ch_versions.mix( EUKCC.out.versions.first() )
    ch_log = ch_log.mix (EUKCC.out.progress_log )


    emit:
    eukcc_csv             = EUKCC.out.eukcc_csv
    eukcc_merged_bins         = EUKCC.out.eukcc_merged_bins
    samtools_idxstats         = LINKTABLE.out.idxstats
    versions                  = ch_versions
    progress_log              = ch_log
}