process CONCATENATE_QUALITY_FILES {
    tag "${meta.id}"

    input:
    tuple val(meta), path(input_files)
    val output_name

    output:
    tuple val(meta), path("${meta.id}.${output_name}"), emit: concatenated_result

    script:
    """
    echo "bin\tcompleteness\tcontamination\tncbi_lng" > "${meta.id}.${output_name}"
    for i in ${input_files}; do
        tail -n +2 \${i} > help_file
        if [ -s help_file ]; then
            grep -v "completeness" \${i} >> "${meta.id}.${output_name}"
        fi
    done
    rm help_file
    """
}

process MODIFY_QUALITY_FILE {
    tag "${output_name}"

    input:
    path(quality_table_csv)
    val output_name

    output:
    path("${output_name}"), emit: modified_result

    script:
    """
    echo "genome,completeness,contamination" > ${output_name}
    grep -v "completeness" ${quality_table_csv} | cut -f1-3 | tr '\t' ',' >> ${output_name} || true
    """
}

process FILTER_QUALITY {
    tag "${meta.id}"

    label 'process_light'

    input:
    tuple val(meta), path(quality_file), path(concoct_bins, stageAs: "concoct_bins/"), path(metabat_bins, stageAs: "metabat_bins/"), path(merged_bins, stageAs: "merged_bins/")

    output:
    tuple val(meta), path("output_genomes/*"), path("quality_file.csv"), emit: qs50_filtered_genomes, optional: true
    path "progress.log"                                                , emit: progress_log

    script:
    """
    mkdir -p output_genomes
    touch quality_file.csv

    echo "Prepare drep quality"
    grep -v "completeness" ${quality_file} |\
    awk '{{if(\$2>=50 && \$2<=100 && \$3>=0 && \$3<=5){{print \$0}}}}' |\
    sort -k 2,3 -n | cut -f1 > filtered_genomes.txt || true

    echo "bins count"
    export BINS=\$(cat filtered_genomes.txt | wc -l)
    echo "\$BINS"
    if [ \$BINS -lt 2 ];
    then
        echo "No genomes"
    else
        for i in \$(ls concoct_bins | grep -w -f filtered_genomes.txt); do
            mv concoct_bins/\${i} output_genomes; done
        for i in \$(ls metabat_bins | grep -w -f filtered_genomes.txt); do
            mv metabat_bins/\${i} output_genomes; done
        for i in \$(ls merged_bins/merged_bins | grep -w -f filtered_genomes.txt); do
            mv merged_bins/merged_bins/\${i} output_genomes; done

        echo "genome,completeness,contamination" > quality_file.csv
        grep -w -f filtered_genomes.txt ${quality_file} | cut -f1-3 | tr '\\t' ',' >> quality_file.csv
    fi

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}
        concoct_bins: \$(ls concoct_bins | wc -l), metabat_bins: \$(ls metabat_bins | wc -l), merged_bins: \$(ls merged_bins/merged_bins | wc -l)
    END_LOGGING
    """
}