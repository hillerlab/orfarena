/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHOOSE                } from '../../../modules/local/choose/main.nf'
include { MAP_CDS               } from '../../../modules/local/mapcds/main.nf'
include { PY_FORMAT_TD2 } from '../../../modules/local/format/td2/main.nf'
include { TD2_LONGORFS  } from '../../../modules/local/td2/longorf/main.nf'
include { TD2_PREDICT  } from '../../../modules/local/td2/predict/main.nf'
include { TD2_FORMAT } from '../../../modules/local/td2/format/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 
workflow TD2 {
    take:
      fasta      // channel: [ val(meta), path(fasta) ]
      reference  // channel: [ val(meta), path(bed) ]

    main:
      ch_versions = Channel.empty()

      TD2_LONGORFS(fasta)
      TD2_PREDICT(fasta, TD2_LONGORFS.out.folder)

      TD2_FORMAT(
        TD2_PREDICT.out.bed,
      )

      MAP_CDS(
        reference,
        TD2_FORMAT.out.bed,
      )

      PY_FORMAT_TD2(
        TD2_PREDICT.out.bed,
        MAP_CDS.out.bed,
        reference,
      )

    emit:
      td2 = PY_FORMAT_TD2.out.td2
      versions    = ch_versions
}
