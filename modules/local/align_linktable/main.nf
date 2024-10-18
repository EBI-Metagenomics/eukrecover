process ALIGNMENT_LINKTABLE {

    /*
    This module aligns reads to the reference using specified arguments and produces output in the form of FASTQ.GZ, BAM and BAM.BAI files.
    It also runs:
     - eukcc linktable
    */

    label 'process_medium'

    tag "${meta.id} align to ${ref_fasta}"

    container 'quay.io/microbiome-informatics/bwa_eukcc:2.2.1_2.0'

    input:
    tuple val(meta), path(reads), path(ref_fasta), path(bins, stageAs: "bins/*"), path(depth)
    val(binner)

    output:
    tuple val(meta), path("*.links.csv"), emit: links_table
    tuple val(meta), path("*.idxstats") , emit: idxstats
    path "versions.yml"                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_args = task.ext.alignment_args

    """
    mkdir -p output

    echo " ---> index fasta"
    bwa-mem2 index ${ref_fasta}

    echo " ---> mapping files to assembly"
    bwa-mem2 mem -M \
      -t ${task.cpus} \
      ${ref_fasta} \
      ${reads} | \
    samtools view -@ ${task.cpus} ${samtools_args} - | \
    samtools sort -@ ${task.cpus} -O bam - -o output/${meta.id}_sorted.bam

    echo "samtools index sorted bam"
    samtools index -@ ${task.cpus} output/${meta.id}_sorted.bam

    echo " ---> samtools idxstats sorted bam"
    samtools idxstats --threads ${task.cpus} output/${meta.id}_sorted.bam > ${prefix}.${binner}.idxstats

    echo "linktable"
    mkdir -p bins
    export BINS=\$(ls bins | grep -v "unbinned" | wc -l)
    if [ \$BINS -eq 0 ]; then
        echo "Bins directory is empty"
        touch ${meta.id}.${binner}.links.csv
    else
        ggp_binlinks.py --ANI 99 \
        --within 1500 \
        --out ${meta.id}.${binner}.links.csv \
        --debug \
        --bindir bins \
        --bam output/${meta.id}_sorted.bam
    fi

    rm -rf output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
        biopython: \$(python3 -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
        pysam: \$(python3 -c "import pkg_resources; print(pkg_resources.get_distribution('pysam').version)")
    END_VERSIONS
    """
}