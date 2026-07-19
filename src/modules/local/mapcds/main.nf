process MAP_CDS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reference)
    tuple val(meta1), path(query)

    output:
    tuple val(meta1), path("*.bed")  , emit: bed
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    map_cds.rs \\
        $args \\
        -r $reference \\
        -q $query

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        map_cds: 0.0.1
    END_VERSIONS
    """

    stub:
    """
    touch output.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        map_cds: 0.0.1
    END_VERSIONS
    """
}
