# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.1] - 2026-07-19

### Added

- Initial Nextflow pipeline scaffold with main workflow, subworkflows, and local modules
- `footprint` Rust crate for ORF footprint analysis, including CLI, BED parsing, record types, and utility functions
- `mapcds` Rust binary for CDS mapping against reference annotations
- `format.py` Python script for handling multiple gene predictor output formats
- `gmst.pl` Perl script for GeneMarkST integration
- Nextflow modules for: blastp, bedtools/getfasta, chromsize, choose, format (genemarkst, psauron, td2, transdecoder), genemarkst (download, predict), gstama (addcdsregions, blastpparser, orfseeker), gunzip, mapcds, psauron, td2 (format, longorf, predict), transdecoder (longorf, predict), ucsc/twobittofa, untar, wget, xloci (cds, exon)
- Subworkflows for: genemarkst, genome, gstama, psauron, td2, td2_strict, transdecoder
- Nextflow configuration with resource profiles and container definitions
- Benchmark and comparison table against xORF, TranslationAI, TD2, TRANSAID, PSAURON, TransDecoder, and GeneMark on the ROSCfam1.0 dog genome
- Test data including genome FASTA and mini BED files
- README with usage instructions, pipeline overview, and benchmark results
- Pipeline architecture diagrams including a detailed Mermaid flowchart covering genome preprocessing, xloci CDS/exon extraction, and all six ORF prediction pathways (TransDecoder, GeneMarkST, TD2, TD2 strict, PSAURON, GSTAMA ORF) with their respective post-processing and CDS mapping steps

### Fixed

- Removed Rust build artifacts from version control that were included in the initial commit
- Updated `.gitignore` to exclude `**/target/` directories and prevent re-inclusion of compiled binaries

[0.0.1]: https://github.com/alejandrogzi/orfarena/releases/tag/v0.0.1
