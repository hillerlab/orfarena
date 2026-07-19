/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GSTAMA ORF AND NMD PREDICTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDTOOLS_GETFASTA      } from '../../../modules/local/bedtools/getfasta/main'
include { GSTAMA_ORF_SEEKER      } from '../../../modules/local/gstama/orfseeker/main'
include { BLASTP                 } from '../../../modules/local/blastp/main'
include { GSTAMA_BLASTP_PARSER   } from '../../../modules/local/gstama/blastpparser/main'
include { GSTAMA_ADD_CDS_REGIONS } from '../../../modules/local/gstama/addcdsregions/main'

workflow GSTAMA_ORF {

    take:
    bed          // channel: [ val(meta), path(bed) ] 
    fasta        // channel: [ val(meta), path(fasta) ]
    protein_db   // channel: [ val(meta), path(database files), val(database basename) ]
    chunk_size   // integer: number of ORF sequences in each BLASTP query chunk

    main:
    ch_versions = Channel.empty()

    BEDTOOLS_GETFASTA(bed, fasta)
    GSTAMA_ORF_SEEKER(BEDTOOLS_GETFASTA.out.fasta)

    GSTAMA_ORF_SEEKER.out.fasta
        .splitFasta(
            by: chunk_size,
            file: 'chunk'
        )
        .map { meta, chunk ->
            def chunk_id = (chunk.name =~ /(chunk\.\d+)/)[0][1]
            def updated_meta = meta + [
                id       : "${meta.id}.${chunk_id}",
                id_former: meta.id
            ]
            [ updated_meta, chunk ]
        }
        .set { ch_orf_chunks }

    BLASTP(ch_orf_chunks, protein_db)

    BLASTP.out.txt
        .map { meta, txt ->
            [ [ id: meta.id_former ], txt ]
        }
        .groupTuple(sort: true)
        .set { ch_blastp_by_sample }

    GSTAMA_BLASTP_PARSER(ch_blastp_by_sample)

    GSTAMA_BLASTP_PARSER.out.tsv
        .join(BEDTOOLS_GETFASTA.out.fasta)
        .join(bed)
        .map { meta, parsed, transcript_fasta, annotation_bed ->
            [ meta, parsed, annotation_bed, transcript_fasta ]
        }
        .set { ch_add_cds_regions }

    GSTAMA_ADD_CDS_REGIONS(ch_add_cds_regions)

    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)
    ch_versions = ch_versions.mix(GSTAMA_ORF_SEEKER.out.versions)
    ch_versions = ch_versions.mix(BLASTP.out.versions)
    ch_versions = ch_versions.mix(GSTAMA_BLASTP_PARSER.out.versions)
    ch_versions = ch_versions.mix(GSTAMA_ADD_CDS_REGIONS.out.versions)

    emit:
    bed      = GSTAMA_ADD_CDS_REGIONS.out.bed
    versions = ch_versions
}
