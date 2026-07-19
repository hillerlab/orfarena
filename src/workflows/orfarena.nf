/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TRANSDECODER } from '../subworkflows/local/transdecoder/main.nf'
include { GENEMARKST } from '../subworkflows/local/genemarkst/main.nf'
include { TD2 } from '../subworkflows/local/td2/main.nf'
include { TD2_STRICT } from '../subworkflows/local/td2_strict/main.nf'
include { PSAURON } from '../subworkflows/local/psauron/main.nf'
include { GSTAMA_ORF } from '../subworkflows/local/gstama/main.nf'
include { UNTAR } from '../modules/local/untar/main.nf'
include { WGET as WGET_BLASTP_DB } from '../modules/local/wget/main.nf'
include { GENOME } from '../subworkflows/local/genome/main.nf'

include { XLOCI_CDS } from '../modules/local/xloci/cds/main.nf'
include { XLOCI_EXON } from '../modules/local/xloci/exon/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 

workflow ORFARENA {
    take:
      reference         // channel: [ val(meta), path(bed) ]
      genome            // val (path)
      revision          // val (str)
      blastp_db         // channel: [ val(meta), path(blastp_db) ]
      gstama_chunk_size // val (int)

    main:
      ch_versions = Channel.empty()

      GENOME(genome)
      genome = GENOME.out.genome  

      XLOCI_CDS(genome, reference)
      XLOCI_EXON(genome, reference)

      WGET_BLASTP_DB(blastp_db)
      UNTAR(WGET_BLASTP_DB.out.outfile)

      UNTAR.out.untar
        .map { meta, db -> [ meta, db, 'uniprot_sprot' ] }
        .set { ch_blastp_db }

      TRANSDECODER(XLOCI_EXON.out.fasta, reference)
      GENEMARKST(XLOCI_EXON.out.fasta, reference, revision)
      TD2(XLOCI_EXON.out.fasta, reference)
      TD2_STRICT(XLOCI_EXON.out.fasta, reference)
      PSAURON(XLOCI_CDS.out.fasta, reference)
      GSTAMA_ORF(reference, genome, ch_blastp_db, gstama_chunk_size)

    emit:
      transdecoder = TRANSDECODER.out.transdecoder
      genemarkst   = GENEMARKST.out.genemarkst
      td2          = TD2.out.td2
      td2_strict   = TD2_STRICT.out.td2
      psauron      = PSAURON.out.psauron
      gstama       = GSTAMA_ORF.out.bed
      versions    = ch_versions
}
