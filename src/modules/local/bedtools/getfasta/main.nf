process BEDTOOLS_GETFASTA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_transcripts"
    """
    bedtools getfasta \\
        -name \\
        -split \\
        -s \\
        -fi $fasta \\
        -bed $bed \\
        -fo ${prefix}.fa \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_transcripts"
    """
    touch ${prefix}.fa
    touch versions.yml
    """
}
