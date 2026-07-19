use clap::Parser;
use genepred::{Bed12, Bed3, GenePred, Reader};
use std::fs::File;
use std::io::Write;
use std::{collections::HashMap, path::PathBuf};

pub const SCALE: u64 = 100_000_000_000; // 100Gb

#[derive(Debug, Clone, Copy)]
enum TranscriptPosition {
    Start,
    End,
}

#[derive(Debug, Parser)]
#[command(author, version=env!("CARGO_PKG_VERSION"), about, long_about = None)]
pub struct Args {
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    #[arg(short = 'q', long = "query")]
    pub query: PathBuf,
}

fn main() {
    let args = Args::parse();
    let mut writer = std::io::BufWriter::new(File::create("output.bed").unwrap());

    let mut refs = HashMap::new();
    Reader::<Bed12>::from_mmap(args.reference)
        .unwrap()
        .records()
        .for_each(|record| {
            let record = record.unwrap();
            let name = record.name.clone().unwrap();

            refs.insert(name, record);
        });

    // let mut query = HashMap::new();
    Reader::<Bed3>::from_mmap(args.query)
        .unwrap()
        .records()
        .for_each(|record| {
            let record = record.unwrap();

            if let Some(ref_record) = refs.get_mut(&record.chrom) {
                let (cds_start, cds_end) = get_cds_from_pos(&ref_record, record.start, record.end);

                let mut extras = HashMap::new();
                extras.insert(
                    "13".as_bytes().to_vec(),
                    genepred::ExtraValue::new_scalar(record.start.to_string().into_bytes()),
                );
                extras.insert(
                    "14".as_bytes().to_vec(),
                    genepred::ExtraValue::new_scalar(record.end.to_string().into_bytes()),
                );

                ref_record.set_thick_start(Some(cds_start));
                ref_record.set_thick_end(Some(cds_end));
                ref_record.set_extras(extras);

                let line = ref_record.to_bed_with_additional_fields::<Bed12>(2);
                writer.write_all(&line).unwrap();
                writer.write_all(b"\n").unwrap();
            } else {
                eprintln!(
                    "WARN: could not find record for {}",
                    std::str::from_utf8(&record.chrom).unwrap()
                );
                return;
            }
        });
}

/// Given a position within the concatenated exonic sequence, this function
/// maps it back to its corresponding genomic coordinate.
///
/// The function iterates through the gene's exons, calculating the length of
/// each exon. It keeps a running total of the length of the exons processed so
/// far (`current_pos`). If the target position (`pos`) falls within the current
/// exon, it calculates the offset from the exon's start and returns the absolute
/// genomic coordinate. If the target position is beyond all exons, it returns `None`.
///
/// # Arguments
///
/// * `pos` - A `u64` representing the position within the concatenated exonic sequence (0-indexed).
///
/// # Returns
///
/// An `Option<u64>` which is:
/// * `Some(genomic_position)` if the position is found within the exons.
/// * `None` if the position is outside the total length of all exons.
///
/// # Panics
///
/// This function does not panic under normal conditions.
///
fn get_pos_in_exons(record: &GenePred, pos: u64) -> Option<u64> {
    get_pos_in_exons_by_kind(record, pos, TranscriptPosition::Start)
}

fn get_end_pos_in_exons(record: &GenePred, pos: u64) -> Option<u64> {
    get_pos_in_exons_by_kind(record, pos, TranscriptPosition::End)
}

fn get_pos_in_exons_by_kind(record: &GenePred, pos: u64, kind: TranscriptPosition) -> Option<u64> {
    let mut current_pos = 0;
    let mut exons = record.exons();

    exons = match record.strand() {
        Some(genepred::Strand::Forward) => exons,
        Some(genepred::Strand::Reverse) => {
            // INFO: we scale the coordinates and then sort them
            // INFO: this is equivalent to reversing the exons and reversing the coordinates
            let mut exons = exons
                .iter()
                .map(|(start, end)| (SCALE - end, SCALE - start))
                .collect::<Vec<(u64, u64)>>();

            exons.sort_by(|a, b| a.0.cmp(&b.0));
            exons
        }
        Some(genepred::Strand::Unknown) | None => {
            panic!("ERROR: unexpected strand value: {:?}", record.strand())
        }
    };

    for (exon_start, exon_end) in exons.iter() {
        let block_len = exon_end - exon_start; // INFO: exon length
        let exon_transcript_end = current_pos + block_len;

        let is_in_exon = match kind {
            TranscriptPosition::Start => pos < exon_transcript_end,
            TranscriptPosition::End => pos <= exon_transcript_end,
        };

        if is_in_exon {
            let offset = pos - current_pos;
            return Some(exon_start + offset);
        }

        current_pos += block_len;
    }

    // INFO: edge case -> dispute between 0-based and 1-based coords [always ends]
    if pos >= current_pos {
        return Some(exons.last().unwrap().1);
    }

    None
}

/// Maps the start and end positions of an Open Reading Frame (ORF) from the
/// concatenated exonic sequence to absolute genomic coordinates.
///
/// This function is crucial for converting ORF coordinates, which are relative to the
/// combined sequence of all exons, into the real genomic coordinates (chromosome, start,
/// end). ORF starts use strict exon membership, while ORF ends include exact exon
/// end boundaries because the ORF interval is half-open in transcript space.
///
/// The function also accounts for the gene's strand. For the forward strand, the
/// converted coordinates are returned directly. For the reverse strand, the coordinates
/// are flipped and adjusted based on a `SCALE` value to correctly reflect their position
/// on the reverse complement.
///
/// # Arguments
///
/// * `orf_start` - A `u64` representing the starting position of the ORF within the
///   concatenated exonic sequence.
/// * `orf_end` - A `u64` representing the ending position of the ORF within the
///   concatenated exonic sequence.
///
/// # Returns
///
/// A tuple `(u64, u64)` containing the absolute genomic start and end coordinates of the ORF.
///
/// # Panics
///
/// This function will panic if `get_pos_in_exons` returns `None` for either `orf_start`
/// or `orf_end`, indicating that the ORF coordinates are outside the exonic sequence.
///
pub fn get_cds_from_pos(record: &GenePred, orf_start: u64, orf_end: u64) -> (u64, u64) {
    let cds_start = get_pos_in_exons(record, orf_start).unwrap_or_else(|| {
        panic!("ERROR: could not map ORF_START to exonic coordinates: {orf_start} -> {record:?}")
    });

    let cds_end = get_end_pos_in_exons(record, orf_end).unwrap_or_else(|| {
        panic!("ERROR: could not map ORF_END to exonic coordinates: {orf_end} -> {record:?}")
    });

    match record.strand() {
        Some(genepred::Strand::Forward) => (cds_start, cds_end),
        Some(genepred::Strand::Reverse) => (SCALE - cds_end, SCALE - cds_start),
        Some(genepred::Strand::Unknown) | None => {
            panic!("ERROR: unexpected strand value: {:?}", record.strand())
        }
    }
}

/// Maps ORF transcript coordinates to absolute genomic coordinates (strand-aware).
///
/// This function takes ORF coordinates defined in **dataar transcript space**
/// (as produced by a CDS/ORF predictor) and returns the corresponding genomic
/// start and end coordinates of the coding region, accounting for exon structure
/// and strand orientation. For transcripts on the reverse strand, it assumes
/// coordinates have been **strand-normalized** using a fixed scale (`SCALE`) and
/// reverses them appropriately.
///
/// # Arguments
///
/// * `orf_start` - The start of the ORF in transcript (spliced) coordinates.
/// * `orf_end` - The end of the ORF in transcript (spliced) coordinates. This is exclusive.
///
/// # Returns
///
/// * `(u64, u64)` - A tuple of genomic start and end coordinates corresponding
///   to the CDS region. Returns `(start, end)` if mapping succeeds; otherwise,
///   may return `Error` if the coordinates do not map within the exon structure.
///
/// # Strand Behavior
///
/// * On the `Strand::Forward`, transcript coordinates map left-to-right across exons.
/// * On the `Strand::Reverse`, the exons are reversed and mapped using:
///   `SCALE - (exon_end - offset)`, where `SCALE` is a user-defined upper genomic bound
///   used for coordinate normalization.
///
/// # Example
///
/// ```rust, ignore
/// let mut gene = GenePred { /* initialized */ };
/// let orf_start = 15;
/// let orf_end = 45;
/// if let Some((cds_start, cds_end)) = gene.map_absolute_cds(orf_start, orf_end) {
///     println!("Genomic CDS coordinates: {} - {}", cds_start, cds_end);
/// }
/// ```
pub fn map_absolute_cds(record: &GenePred, orf_start: u64, orf_end: u64) -> (u64, u64) {
    // INFO: for Reverse this should be interpreted from 5' to 3'
    let (cds_start, cds_end) = get_cds_from_pos(record, orf_start, orf_end);

    assert!(
        cds_start > 0,
        "ERROR: {orf_start} is not bigger than 0! -> {record:?}"
    );
    assert!(
        cds_end > 0,
        "ERROR: {orf_end} is not bigget than 0! -> {record:?}"
    );

    assert!(
        cds_start >= record.start,
        "ERROR: {cds_start} is less than {}! Mapped ORF values {orf_start} - {orf_end} -> {record:?}",
        record.start
    );
    assert!(
        cds_end <= record.end,
        "ERROR: {cds_end} is greater than {}! Mapped ORF values {orf_start} - {orf_end} -> {record:?}",
        record.end
    );

    (cds_start, cds_end)
}

#[cfg(test)]
mod tests {
    use super::*;
    use genepred::{Bed12, Reader};

    #[test]
    fn test_absolute_cds_mapping_forward() {
        let data = "chr6\t8259278\t8593709\tR441_chr6__FC23#TC31#PA0#PR0#IY1000\t60\t+\t8259278\t8593709\t43,118,219\t10\t215,109,62,215,152,87,117,153,211,2165\t0,11194,114627,167692,278538,299221,313903,320308,323301,332266";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 237;
        let orf_end = 420;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start, orf_end);

        assert_eq!(predicted_cds_start, 8270494);
        assert_eq!(predicted_cds_end, 8427004);
    }

    #[test]
    fn test_absolute_cds_mapping_reverse() {
        let data = "chr7\t45482351\t45520392\tR206671_chr7__FC0#TC23#PA0#PR0#IY998\t60\t-\t45482351\t45520392\t43,118,219\t13\t1162,233,188,161,232,126,154,171,210,117,620,482,1256\t0,1321,1711,20275,21300,21646,23812,24550,24947,25522,29018,33186,36785";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 1284;
        let orf_end = 3321;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start, orf_end);
        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 45503698);
        assert_eq!(predicted_cds_end, 45515991);
    }

    #[test]
    fn test_get_cds_from_pos_forward() {
        let data = "chr12\t102521030\t102531285\tR909465_chr12__FC30#TC0#PA164#PR164#IY896\t60\t+\t102521030\t102531285\t43,118,219\t8\t399,47,94,69,159,417,476,417\t0,1113,3694,4517,6573,6895,7826,9838";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 1530;
        let orf_end = 1674;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        assert_eq!(predicted_cds_start, 102529201);
        assert_eq!(predicted_cds_end, 102530881);
    }

    #[test]
    fn test_get_pos_in_exons_forward() {
        let data = "chr12\t102521030\t102531285\tR909465_chr12__FC30#TC0#PA164#PR164#IY896\t60\t+\t102521030\t102531285\t43,118,219\t8\t399,47,94,69,159,417,476,417\t0,1113,3694,4517,6573,6895,7826,9838";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 1530;
        let orf_end = 1674;

        let predicted_cds_start = get_pos_in_exons(&record, orf_start);
        let predicted_cds_end = get_pos_in_exons(&record, orf_end);

        assert_eq!(predicted_cds_start, Some(102529201));
        assert_eq!(predicted_cds_end, Some(102530881));
    }

    #[test]
    fn test_get_cds_from_pos_reverse() {
        let data = "chr7\t45482351\t45520392\tR206671_chr7__FC0#TC23#PA0#PR0#IY998\t60\t-\t45482351\t45520392\t43,118,219\t13\t1162,233,188,161,232,126,154,171,210,117,620,482,1256\t0,1321,1711,20275,21300,21646,23812,24550,24947,25522,29018,33186,36785";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 358;
        let orf_end = 553;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 45519839);
        assert_eq!(predicted_cds_end, 45520034);
    }

    #[test]
    fn test_get_pos_in_exons_reverse() {
        let data = "chr7\t45482351\t45520392\tR206671_chr7__FC0#TC23#PA0#PR0#IY998\t60\t-\t45482351\t45520392\t43,118,219\t13\t1162,233,188,161,232,126,154,171,210,117,620,482,1256\t0,1321,1711,20275,21300,21646,23812,24550,24947,25522,29018,33186,36785";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 1284;
        let orf_end = 3321;

        let predicted_cds_end = get_pos_in_exons(&record, orf_start).unwrap();
        let predicted_cds_start = get_pos_in_exons(&record, orf_end);

        dbg!(
            SCALE - predicted_cds_start.unwrap(),
            SCALE - predicted_cds_end
        );

        assert_eq!(SCALE - predicted_cds_start.unwrap(), 45503698);
        assert_eq!(SCALE - predicted_cds_end, 45515991);
    }

    #[test]
    fn test_get_pos_in_exons_forward_refseq_tai_prediction() {
        let data = "chr1\t58713285\t58733227\tNM_009805.4\t0\t+\t58726436\t58732362\t0,0,0\t5\t374,427,106,136,975,\t0,13020,15770,17866,18967,";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 505;
        let orf_end = 1150;

        let predicted_cds_start = get_pos_in_exons(&record, orf_start).unwrap();
        let predicted_cds_end = get_pos_in_exons(&record, orf_end).unwrap();

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 58726436);
        assert_eq!(predicted_cds_end + 3, 58732362);
    }

    #[test]
    fn test_get_pos_in_exons_reverse_refseq_tai_prediction() {
        let data = "chr2\t73092800\t73214447\tNM_025942.2\t0\t-\t73093622\t73212960\t0,0,0\t11\t924,123,97,141,98,81,176,128,144,101,162,\t0,4344,6491,7273,49097,49507,63937,106600,110524,120059,121485,";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 162;
        let orf_end = 1350;

        let predicted_cds_end = get_pos_in_exons(&record, orf_start).unwrap();
        let predicted_cds_start = get_pos_in_exons(&record, orf_end);

        dbg!(
            SCALE - predicted_cds_start.unwrap(),
            SCALE - predicted_cds_end
        );

        assert_eq!(SCALE - predicted_cds_start.unwrap() - 3, 73093622);
        assert_eq!(SCALE - predicted_cds_end, 73212960);
    }

    #[test]
    fn test_get_pos_in_exons_reverse_refseq_tai_prediction_additional() {
        // 2025-09-10T11:33:47.207Z WARN  [orf::tai] WARN: translationAi predicted a non-stop ORF: "304,484,0.11812351644039154,0.1313595026731491" for GenePred { nam
        // e: "9249", chrom: "chr7", strand: Reverse, start: 99940345148, end: 99940364371, cds_start: 99940345148, cds_end: 99940364371, exons: [(99940345148, 999403
        // 45217), (99940350839, 99940350935), (99940354698, 99940354847), (99940361505, 99940361651), (99940362689, 99940362841), (99940363369, 99940363481), (999403
        // 64269, 99940364371)], introns: [(99940345218, 99940350838), (99940350936, 99940354697), (99940354848, 99940361504), (99940361652, 99940362688), (9994036284
        // 2, 99940363368), (99940363482, 99940364268)], exon_len: 826, exon_count: 7, data: "chr7\t59635629\t59654852\t9249\t60\t-\t59635629\t59654852\t43,118,219\t7
        // \t102,112,152,146,149,96,69\t0,890,1530,2720,9524,13436,19154", is_ref: false }

        let data = "chr7\t59635629\t59654852\t9249\t60\t-\t59635629\t59654852\t43,118,219\t7\t102,112,152,146,149,96,69\t0,890,1530,2720,9524,13436,19154";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 317;
        let orf_end = 437;

        let (predicted_cds_end, predicted_cds_start) =
            map_absolute_cds(&record, orf_start as u64, orf_end as u64);

        if predicted_cds_start - 3 < record.start {
            dbg!("WARN: translationAi predicted a non-stop ORF: {orf_start:?} for {?}");
        }

        if predicted_cds_end + 3 > record.end {
            dbg!("WARN: translationAi predicted a non-stop ORF: {orf_start:?} for {?}");
        }

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start - 3, 59638489);
        assert_eq!(predicted_cds_end, 59638372);
    }

    #[test]
    fn test_get_pos_in_exons_reverse_bug_prediction_additional_2() {
        let data = "scaffold_10\t60508574\t60652099\tR11900_scaffold_10__FC0#TC0#PA29#PR30#IY979#SG\t60\t-\t60508574\t60652099\t51,153,255\t4\t499,38,310,77\t0,4353,112181,143448";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 28;
        let orf_end = 748;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start as u64, orf_end as u64);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 60508750);
        assert_eq!(predicted_cds_end, 60652071);

        let orf_start = 52;
        let orf_end = 748;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start as u64, orf_end as u64);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 60508750);
        assert_eq!(predicted_cds_end, 60652047);
    }

    #[test]
    fn test_get_pos_in_exons_reverse_bug_prediction_additional_3() {
        let data = "scaffold_17\t80102749\t80103945\tR16258_scaffold_17__FC0#TC0#PA2#PR28#IY971\t0\t+\t80102749\t80103945\t43,118,219\t1\t1196\t0";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 84;
        let orf_end = 483;

        let (predicted_cds_start, predicted_cds_end) =
            map_absolute_cds(&record, orf_start as u64, orf_end as u64);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 80102833);
        assert_eq!(predicted_cds_end, 80103232);
    }

    #[test]
    fn test_get_pos_in_exons_reverse_refseq_orfipy_prediction() {
        let data = "chr5\t151561660\t151586924\tNM_001102582.1\t0\t-\t151561660\t151586906\t0,0,0\t5\t911,124,228,836,286,\t0,11135,14005,22893,24978,";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_end = 2385;

        let predicted_cds_start = get_pos_in_exons(&record, orf_end);

        dbg!(SCALE - predicted_cds_start.unwrap());

        assert_eq!(SCALE - predicted_cds_start.unwrap(), 151561660);
    }

    #[test]
    fn test_get_pos_in_exons_forward_refseq_orfipy_prediction() {
        let data = "chr1\t83116765\t83119167\tNM_001159738.1\t0\t+\t83116823\t83118717\t0,0,0\t4\t137,112,78,472,\t0,1033,1233,1930,";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 46;
        let orf_end = 349;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 83116811);
        assert_eq!(predicted_cds_end, 83118717);
    }

    #[test]
    fn test_get_pos_in_exons_transdecoder_prediction() {
        let data = "NC_037352	388166	396055	NM_001078122.2	0	+	388258	394235	0,0,0	9	172,69,52,136,125,94,105,215,948,	0,1303,1703,2322,4865,5058,5303,5854,6941,";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 92;
        let orf_end = 968;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 388258);
        assert_eq!(predicted_cds_end, 394235);
    }

    #[test]
    fn test_get_pos_in_exons_forward_translationai_prediction_wrong_start_cow() {
        let data = "NC_037328	102267973	102388511	XM_005201715.5	0	+	102272281	102388043	0,0,0	48	96,118,95,118,110,152,172,100,113,126,132,120,114,85,118,172,117,155,85,57,125,89,50,171,156,204,155,169,97,113,126,129,111,63,117,88,139,134,35,117,149,85,57,125,89,50,168,537	0,4308,5644,7568,7981,9021,10662,12452,13630,15870,18181,18399,20348,29163,32745,34854,38566,40588,41701,42678,43435,45403,48668,50006,58153,68412,70214,71796,73546,73750,73952,75545,77369,79915,81921,83444,84463,89302,90352,90474,93758,96861,97028,97687,100948,116510,117018,120001";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 96;
        let orf_end = 5535;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 102272281);
        assert_eq!(predicted_cds_end, 102388043);
    }

    #[test]
    fn test_get_pos_in_exons_forward_translationai_prediction_wrong_start_cow_additional() {
        let data = "NC_037328	127530720	127609658	XM_059888940	0	-	127533976	127537591	0,0,0	7	6871,105,11028,446,66,75,309	0,44010,46138,58325,59275,77446,78629";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 12029;
        let orf_end = 15644;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 127533976);
        assert_eq!(predicted_cds_end, 127537591);
    }

    #[test]
    fn test_get_pos_in_exons_forward_translationai_prediction_wrong_start_cow_additional_2() {
        let data = "NC_037328	69329582	69357506	XM_002684796.7	0	-	69332533	69357372	0,0,0	12	923,196,85,38,121,126,164,56,107,120,1074,186	0,2951,5733,7446,8130,8564,11264,14986,18831,20186,24458,27738";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 134;
        let orf_end = 2273;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 69332533);
        assert_eq!(predicted_cds_end, 69357372);
    }

    #[test]
    fn test_get_pos_in_exons_forward_translationai_prediction_wrong_start_cow_additional_3() {
        let data = "NC_037357	35186707	35245555	XM_024988460.2	0	+	35201872	35242092	0,0,0	22	158,79,66,42,78,90,87,84,42,87,105,99,96,90,195,195,111,81,72,103,99,2534	0,12066,12624,14433,14577,15164,16087,16869,17192,23711,26945,28746,31210,32094,38359,43644,46231,47344,49432,50523,55286,56314";

        let mut reader: Reader<Bed12> =
            Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
        let record = reader.next().unwrap().unwrap();

        let orf_start = 424;
        let orf_end = 2059;

        let (predicted_cds_start, predicted_cds_end) =
            get_cds_from_pos(&record, orf_start, orf_end);

        dbg!(predicted_cds_start, predicted_cds_end);

        assert_eq!(predicted_cds_start, 35201872);
        assert_eq!(predicted_cds_end, 35242092);
    }
}
