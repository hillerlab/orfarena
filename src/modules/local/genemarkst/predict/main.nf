process GENEMARKST_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
      ? 'docker://perl:5.40.4-slim-bookworm'
      : 'docker.io/library/perl:5.40.4-bookworm'}"


    input:
    tuple val(meta), path(fasta)
    path gmst

    output:
    tuple val(meta), path("*.cds.gff")  , emit: gff
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    perl gmst.pl \\
        $args \\
        $fasta \\
        --format GFF

    grep -w 'CDS' *.gff > ${prefix}.cds.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gmst:  0.0.1
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cds.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gmst:  0.0.1
    END_VERSIONS
    """
}
