#!/usr/bin/env nextflow

/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    orfarena

    Benchmark ORF prediction methods
    Authors: Alejandro Gonzales-Irribarren, Michael Hiller
    GitHub:  https://github.com/hillerlab/orfarena
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS AND MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ORFARENA as RUN } from './workflows/orfarena.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATION FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateFullRun() {
    def problems = []

    if (!params.bed)             { problems << 'missing required --bed' }
    if (!params.genome)        { problems << 'missing required --genome (genome)' }

    if (problems) {
        error "Parameter validation failed:\n  - " + problems.join('\n  - ')
        System.exit(1)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ORFARENA {
    validateFullRun()

    log.info """
    orfarena v${workflow.manifest.version}
  
    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

      Transcripts : ${params.bed} 
      Genome      : ${params.genome}
      Outdir      : ${params.outdir}
      Profile     : ${workflow.profile}
    """.stripIndent()

    ch_reference = Channel.value(file(params.bed, checkIfExists: true)).map { it -> [ [ id: it.baseName ], it ] }
    ch_blastp_db = Channel.value(params.blastp_db).map { it -> [ [ id: 'uniprot_sprot.tar.gz' ], it ] }

    RUN(
      ch_reference, 
      params.genome, 
      params.gmst_revision, 
      ch_blastp_db, 
      params.gstama_chunk_size
    )
}

workflow {
    ORFARENA()
}

workflow.onComplete {
    log.info(workflow.success
        ? "\norfarena finished. Results in: ${params.outdir}\n"
        : "\norfarena FAILED. See ${params.outdir}/pipeline_info for reports.\n")
}
