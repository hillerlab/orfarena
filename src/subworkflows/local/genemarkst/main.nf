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
include { PY_FORMAT_GENEMARKST } from '../../../modules/local/format/genemarkst/main.nf'
include { GENEMARKST_PREDICT    } from '../../../modules/local/genemarkst/predict/main.nf'
include { DOWNLOAD_GMST         } from '../../../modules/local/genemarkst/download/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 
workflow GENEMARKST {
    take:
      fasta      // channel: [ val(meta), path(fasta) ]
      reference  // channel: [ val(meta), path(bed) ]
      revision   // val

    main:
      ch_versions = Channel.empty()

      DOWNLOAD_GMST(revision)
      GENEMARKST_PREDICT(fasta, DOWNLOAD_GMST.out.gmst)

      CHOOSE(
        GENEMARKST_PREDICT.out.gff,
        "0 3 4"
      )

      MAP_CDS(
        reference,
        CHOOSE.out.output,
      )

      PY_FORMAT_GENEMARKST(
        GENEMARKST_PREDICT.out.gff, 
        MAP_CDS.out.bed, 
        reference
      )

    emit:
      genemarkst = PY_FORMAT_GENEMARKST.out.genemarkst
      versions    = ch_versions
}
