/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PY_FORMAT_PSAURON } from '../../../modules/local/format/psauron/main.nf'
include { PSAURON_PREDICT  } from '../../../modules/local/psauron/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 
workflow PSAURON {
    take:
      fasta      // channel: [ val(meta), path(fasta) ]
      reference  // channel: [ val(meta), path(bed) ]

    main:
      ch_versions = Channel.empty()

      PSAURON_PREDICT(fasta)

      PY_FORMAT_PSAURON(
        PSAURON_PREDICT.out.psauron,
        reference,
      )

    emit:
      psauron = PY_FORMAT_PSAURON.out.psauron
      versions    = ch_versions
}
