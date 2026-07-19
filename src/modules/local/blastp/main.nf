process BLASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.17.0--h66d330f_0' :
        'quay.io/biocontainers/blast:2.17.0--h66d330f_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(db), val(db_name)

    output:
    tuple val(meta), path("*.blastp.txt"), emit: txt
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blastp \\
        -query $fasta \\
        -db ${db}/${db_name} \\
        -out ${prefix}.blastp.txt \\
        -evalue 1e-10 \\
        -num_threads $task.cpus \\
        -outfmt 0 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastp: \$(blastp -version 2>&1 | head -n1 | cut -d ' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.blastp.txt
    touch versions.yml
    """
}
