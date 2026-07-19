process DOWNLOAD_GMST {
    tag "SQANTI2 gmst: ${revision}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/git:latest' }"

    input:
    val revision

    output:
    path "gmst/*",       emit: gmst
    path "versions.yml", emit: versions

    script:
    """
    git init sqanti2

    git -C sqanti2 remote add origin \
        https://github.com/Magdoll/SQANTI2.git

    git -C sqanti2 config core.sparseCheckout true

    printf '/utilities/gmst/\\n' \
        > sqanti2/.git/info/sparse-checkout

    git -C sqanti2 fetch \
        --depth 1 \
        origin \
        ${revision}

    git -C sqanti2 checkout \
        --detach \
        FETCH_HEAD

    commit=\$(git -C sqanti2 rev-parse HEAD)

    mv sqanti2/utilities/gmst gmst
    rm -rf sqanti2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti2: "\${commit}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p gmst 
    touch gmst/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti2: "0.0.1"
    END_VERSIONS
    """
}
