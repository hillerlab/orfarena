process PY_FORMAT_PSAURON {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(csv)
    tuple val(meta1), path(ground)

    output:
    tuple val(meta), path("*.bed")  , emit: psauron

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    format.py \
        psauron \
        --psauron $csv \
        --ground $ground \
        --outdir . \
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        format: \$(format.py --version | sed 's/format.py //g')
    END_VERSIONS
    """

    stub:
    """
    touch *.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        format: \$(format.py --version | sed 's/format.py //g')
    END_VERSIONS
    """
}
