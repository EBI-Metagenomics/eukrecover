/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EUKCC as EUKCC_CONCOCT } from '../modules/eukcc'
include { EUKCC as EUKCC_METABAT } from '../modules/eukcc'
include { LINKTABLE as LINKTABLE_CONCOCT } from '../modules/eukcc'
include { LINKTABLE as LINKTABLE_METABAT } from '../modules/eukcc'

workflow EUK_SUBWF {
    take:
        bins_concoct
        bins_metabat
        bam
        eukcc_db
    main:
        # concoct
        binner1 = channel.value("concoct")
        LINKTABLE_CONCOCT(bins_concoct, bam)
        EUKCC_CONCOCT(binner1, LINKTABLE.out.links_table, eukcc_db.first(), bins_concoct)

        # metabat2
        binner2 = channel.value("metabat2")
        LINKTABLE_METABAT(bins_metabat, bam)
        EUKCC_METABAT(binner2, LINKTABLE.out.links_table, eukcc_db.first(), bins_metabat)

        # prepare quality file
        combine_quality = EUKCC_CONCOCT.out.eukcc_csv.combine(EUKCC_METABAT.out.eukcc_csv, by: 0)
        combine_quality.view()
        //quality = EUKCC_CONCOCT.out.eukcc_csv.map{
        //    item -> item[1].combine(item[2]).text().into {
        //        file -> file.write("genome,completeness,contamination\n") }
        //        }
    emit:
        euk_quality = EUKCC_CONCOCT.out.eukcc_csv
}