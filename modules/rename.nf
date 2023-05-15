/*
 * dRep, this workflow dereplicates a set of genomes.
*/
process RENAME {

    publishDir(
        path: "${params.outdir}/drep",
        mode: 'copy',
    )

    container 'quay.io/biocontainers/drep:3.2.2--pyhdfd78af_0'

    input:
    path fasta

    output:
    path "*_moved", emit: test

    script:
    """
    mv ${fasta} ${fasta}_moved
    """

}
