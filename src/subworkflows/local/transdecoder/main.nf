/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TRANSDECODER_LONGORF  } from '../../../modules/local/transdecoder/longorf/main.nf'
include { TRANSDECODER_PREDICT  } from '../../../modules/local/transdecoder/predict/main.nf'
include { CHOOSE                } from '../../../modules/local/choose/main.nf'
include { MAP_CDS               } from '../../../modules/local/mapcds/main.nf'
include { PY_FORMAT_TRANSDECODER } from '../../../modules/local/format/transdecoder/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 
workflow TRANSDECODER {
    take:
      fasta      // channel: [ val(meta), path(fasta) ]
      reference  // channel: [ val(meta), path(bed) ]

    main:
      ch_versions = Channel.empty()

      TRANSDECODER_LONGORF(fasta)
      TRANSDECODER_PREDICT(fasta, TRANSDECODER_LONGORF.out.folder)

      CHOOSE(
        TRANSDECODER_PREDICT.out.bed,
        "0 6 7"
      )

      MAP_CDS(
        reference,
        CHOOSE.out.output,
      )

      PY_FORMAT_TRANSDECODER(
        TRANSDECODER_PREDICT.out.bed, 
        MAP_CDS.out.bed, 
        reference
      )

    emit:
      transdecoder = PY_FORMAT_TRANSDECODER.out.transdecoder
      versions    = ch_versions
}
