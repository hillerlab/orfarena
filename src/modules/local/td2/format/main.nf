/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TD2_FORMAT — Reformat TransDecoder2 BED output into a three-column BED file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TD2_FORMAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.0--2' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.coding.td2.reformatted.bed"), emit: bed
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python - '${bed}' '${prefix}.coding.td2.reformatted.bed' <<'PY'
    import csv
    import sys
    from pathlib import Path


    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    with input_path.open("r", newline="") as input_handle, \\
         output_path.open("w", newline="") as output_handle:

        reader = csv.reader(input_handle, delimiter="\\t")
        writer = csv.writer(
            output_handle,
            delimiter="\\t",
            lineterminator="\\n",
        )

        for line_number, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#") or row[0].startswith("track"):
                continue

            if len(row) < 4:
                raise ValueError(
                    f"Line {line_number} in {input_path} contains fewer than "
                    f"four BED columns: {row}"
                )

            name = row[3]

            try:
                coordinates = name.rsplit(":", 1)[-1]
                coordinates = coordinates.split("(", 1)[0]
                coordinate_fields = coordinates.split("-")

                orf_start = coordinate_fields[0]
                orf_end = coordinate_fields[-1]

                # Validate that the extracted coordinates are integers.
                int(orf_start)
                int(orf_end)

            except (IndexError, ValueError) as error:
                raise ValueError(
                    f"Could not extract ORF coordinates from BED name "
                    f"on line {line_number}: {name!r}"
                ) from error

            writer.writerow([
                row[0],
                orf_start,
                orf_end,
            ])
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.coding.td2.reformatted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
