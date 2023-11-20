process BINNING_REFINER {

    tag "${name} ${meta.id}"

    container 'quay.io/biocontainers/biopython:1.75'

    publishDir(
        path: "${params.outdir}/intermediate_steps/refinement/${meta.id}_binref_${name}",
        mode: params.publish_dir_mode,
        failOnError: true,
        pattern: "{meta.id}_output_${name}/Refined/*"
    )

    input:
    val(name)
    tuple val(meta), path(bin1, stageAs: "binner1/*"), path(bin2, stageAs: "binner2/*")

    output:
    tuple val(meta), path("${meta.id}_output_${name}/Refined/*"), emit: refined_bins
    path "versions.yml"                                         , emit: versions

    script:
    """
    binning_refiner.py -1 binner1 -2 binner2 -o "${meta.id}_output_${name}" -n ${name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}

process BINNING_REFINER3 {

    tag "${name} ${meta.id}"

    container 'quay.io/biocontainers/biopython:1.75'

    publishDir(
        path: "${params.outdir}/intermediate_steps/refinement/${meta.id}_binref_${name}",
        mode: 'copy',
        failOnError: true,
        pattern: "{meta.id}_output_${name}/Refined/*"
    )

    input:
    val(name)
    tuple val(meta), 
          path(bin1, stageAs: "binner1/*"),
          path(bin2, stageAs: "binner2/*"),
          path(bin3, stageAs: "binner3/*")

    output:
    tuple val(meta), path("${meta.id}_output_${name}/Refined/*"), emit: refined_bins
    path "versions.yml"                                         , emit: versions

    script:
    """
    binning_refiner.py \
    -1 binner1 \
    -2 binner2 \
    -3 binner3 \
    -o "${meta.id}_output_${name}" \
    -n ${name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}