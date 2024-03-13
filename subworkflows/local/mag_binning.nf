/*
 * Binning with MetaBAT2, MaxBin2 and Concoct
 */
include { METABAT2_METABAT2                     } from '../../modules/nf-core/metabat2/metabat2/main'
include { MAXBIN2                               } from '../../modules/nf-core/maxbin2/main'
include { CONVERT_DEPTHS                        } from '../../modules/local/mag/convert_depths'
include { FASTA_BINNING_CONCOCT                 } from '../nf-core/fasta_binning_concoct/main'

workflow BINNING {

    take:
    binning_input  // channel: [ val(meta), path(assembly), path(depth), [reads] ]

    main:

    ch_versions = Channel.empty()
    // optional //
    metabat_output = Channel.empty()
    concoct_output = Channel.empty()

    assembly = binning_input.map{ meta, assembly, depth, reads -> [meta, assembly] }
    depth = binning_input.map{ meta, assembly, depth, reads -> [meta, depth] }

    if ( !params.skip_maxbin2 ) {
        
        CONVERT_DEPTHS ( ch_metabat2_input )

        ch_maxbin2_input = CONVERT_DEPTHS.out.output.map { meta, assembly, depth ->
            // we provide empty reads as we don't want maxbin2 to calculate for us.
            [meta.subMap('id', 'erz', 'single_end') + [binner: 'MaxBin2'], assembly, [], depth ]
        }

        MAXBIN2 ( ch_maxbin2_input )  // input (meta, assembly, reads, depth); output can be empty folder

        maxbin_output = MAXBIN2.out.binned_fastas.map { meta, bins ->
            [ meta.subMap('id', 'erz', 'single_end'), bins ]}

        ch_versions = ch_versions.mix( CONVERT_DEPTHS.out.versions.first() )
        ch_versions = ch_versions.mix( MAXBIN2.out.versions.first() )
    }

    if ( !params.skip_metabat2 ) {
        ch_metabat2_input = assembly.join(depth).map { meta, assembly, depth ->
                                                    [meta + [binner: 'MetaBAT2'], assembly, depth]
                                                    }
        METABAT2_METABAT2 ( ch_metabat2_input )   // input: (meta, fasta, depth)

        metabat_output = METABAT2_METABAT2.out.fasta.map { meta, bins ->
            [ meta.subMap('id', 'erz', 'single_end'), bins ]}

        ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())
    }

    if ( !params.skip_concoct ) {

        assemblies_bams.map { meta, assembly, bams, bais ->
            [ meta + [binner: 'CONCOCT'], assembly, bams, bais ]
        }.multiMap { meta_extended, assembly, bams, bais ->
            bins: [ meta_extended, assembly ]
            bams: [ meta_extended, bams, bais ]
        }.set { ch_concoct_input }

        FASTA_BINNING_CONCOCT( ch_concoct_input.bins, ch_concoct_input.bams )

        concoct_output = FASTA_BINNING_CONCOCT.out.bins.map { meta, bins ->
            [ meta.subMap('id', 'erz', 'single_end'), bins ]}

        ch_versions = ch_versions.mix( FASTA_BINNING_CONCOCT.out.versions.first() )
    }

    emit:
    metabat2depths   = depth
    maxbin_bins      = maxbin_output
    concoct_bins     = concoct_output
    metabat_bins     = metabat_output
    versions         = ch_versions
}