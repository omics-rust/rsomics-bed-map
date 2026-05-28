//! Aggregate column values from B intervals overlapping each A interval.
//!
//! Equivalent to `bedtools map -a A -b B -c COL -o OP`: for each A interval,
//! all B intervals that overlap it are found and the specified column value
//! from each B hit is aggregated using the given operation.
//!
//! Supported operations: sum, min, max, mean, count, count_distinct, median,
//! mode, collapse (comma-join), distinct (comma-join unique), absmin, absmax.
//! Null value (when no B hits) is reported as `.` matching bedtools default.
//!
//! Multiple columns and operations may be specified as comma-separated lists
//! (matching bedtools map -c 4,5 -o sum,min behaviour).
//!
//! Algorithm: B is loaded into per-chromosome Vecs (sorted by start) storing
//! field byte-offsets into a shared raw buffer — no per-field String allocation.
//! For each A record the upper-bound of overlapping B records is found with a
//! binary search; within that range only records with `end > a_start` contribute.
//! Hit fields are collected into a reusable `Vec<Vec<u8>>` per ColOp, cleared
//! between A records, giving O((N+M) log M) total with low allocation pressure.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// Aggregate operation applied to a set of column values from overlapping B records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Sum of numeric values.
    Sum,
    /// Minimum numeric value.
    Min,
    /// Maximum numeric value.
    Max,
    /// Arithmetic mean of numeric values.
    Mean,
    /// Count of overlapping B records.
    Count,
    /// Count of distinct values.
    CountDistinct,
    /// Median numeric value.
    Median,
    /// Mode (most frequent value, first if tied).
    Mode,
    /// Comma-separated list of all values.
    Collapse,
    /// Comma-separated list of distinct values (in order of first occurrence).
    Distinct,
    /// Minimum of absolute values.
    AbsMin,
    /// Maximum of absolute values.
    AbsMax,
}

impl Op {
    pub fn parse_op(s: &str) -> Option<Self> {
        match s {
            "sum" => Some(Self::Sum),
            "min" => Some(Self::Min),
            "max" => Some(Self::Max),
            "mean" => Some(Self::Mean),
            "count" => Some(Self::Count),
            "count_distinct" => Some(Self::CountDistinct),
            "median" => Some(Self::Median),
            "mode" => Some(Self::Mode),
            "collapse" => Some(Self::Collapse),
            "distinct" => Some(Self::Distinct),
            "absmin" => Some(Self::AbsMin),
            "absmax" => Some(Self::AbsMax),
            _ => None,
        }
    }
}

/// A column + operation pairing specifying one output column in the result.
#[derive(Debug, Clone)]
pub struct ColOp {
    /// 1-based column index into the B record (1 = chrom, 2 = start, 3 = end, 4+ = data).
    pub col: usize,
    pub op: Op,
}

/// A B record: integer coordinates + field boundary offsets into the shared raw buffer.
#[derive(Debug)]
struct BRecord {
    start: u64,
    end: u64,
    /// Byte offset of this record's line start in `BChrom::buf`.
    line_start: u32,
    /// For each field i, `field_ends[i]` is the exclusive byte offset of that field
    /// in `BChrom::buf`. Field 0 spans `[line_start, field_ends[0])`, etc.
    field_ends: Vec<u32>,
}

/// Raw byte buffer + sorted record index for one chromosome.
struct BChrom {
    /// Concatenated lines (no newlines).
    buf: Vec<u8>,
    /// Records sorted by start.
    records: Vec<BRecord>,
    /// Maximum (end − start) across all records. Enables lower-bound binary search.
    max_feat_size: u64,
}

impl BChrom {
    /// Return the raw bytes of 1-based column `col` from `rec`.
    fn field_bytes<'a>(&'a self, rec: &'a BRecord, col: usize) -> &'a [u8] {
        let n = rec.field_ends.len();
        if col == 0 || col > n {
            return b"";
        }
        let start = if col == 1 {
            rec.line_start as usize
        } else {
            rec.field_ends[col - 2] as usize + 1 // +1 to skip the tab separator
        };
        let end = rec.field_ends[col - 1] as usize;
        &self.buf[start..end]
    }
}

#[inline]
fn is_skippable(line: &[u8]) -> bool {
    line.is_empty()
        || line.starts_with(b"#")
        || line.starts_with(b"track")
        || line.starts_with(b"browser")
}

/// Load B BED into per-chromosome structures. Uses a single byte buffer per
/// chromosome; records store byte-offsets rather than cloned Strings.
fn load_b(path: &Path) -> Result<HashMap<String, BChrom>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::with_capacity(256 * 1024, file);

    let mut chrom_data: HashMap<String, (Vec<u8>, Vec<BRecord>, u64)> = HashMap::new();

    for line_result in reader.split(b'\n') {
        let line = line_result.map_err(RsomicsError::Io)?;
        if is_skippable(&line) {
            continue;
        }

        // Locate first three tab positions to extract chrom, start, end cheaply.
        let mut tab_pos = [0usize; 3];
        let mut ntabs = 0usize;
        for (i, &b) in line.iter().enumerate() {
            if b == b'\t' {
                if ntabs < 3 {
                    tab_pos[ntabs] = i;
                }
                ntabs += 1;
            }
        }
        if ntabs < 2 {
            continue; // malformed line
        }
        let t1 = tab_pos[0];
        let t2 = tab_pos[1];
        let t3 = if ntabs >= 3 { tab_pos[2] } else { line.len() };

        let chrom = std::str::from_utf8(&line[..t1])
            .map_err(|_| RsomicsError::InvalidInput("B: non-UTF8 chrom".into()))?;
        let start: u64 = std::str::from_utf8(&line[t1 + 1..t2])
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let end: u64 = std::str::from_utf8(&line[t2 + 1..t3])
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let (buf, records, max_feat_size) = chrom_data.entry(chrom.to_owned()).or_default();
        let feat_size = end.saturating_sub(start);
        if feat_size > *max_feat_size {
            *max_feat_size = feat_size;
        }
        let line_start = buf.len() as u32;
        buf.extend_from_slice(&line);

        // Build field_ends by scanning for tabs in this line.
        let mut field_ends: Vec<u32> = Vec::with_capacity(6);
        let mut pos: u32 = line_start;
        for &b in &line {
            pos += 1;
            if b == b'\t' {
                field_ends.push(pos - 1);
            }
        }
        field_ends.push(pos); // last field

        records.push(BRecord {
            start,
            end,
            line_start,
            field_ends,
        });
    }

    let mut result: HashMap<String, BChrom> = HashMap::with_capacity(chrom_data.len());
    for (chrom, (buf, mut records, max_feat_size)) in chrom_data {
        records.sort_unstable_by_key(|r| r.start);
        result.insert(
            chrom,
            BChrom {
                buf,
                records,
                max_feat_size,
            },
        );
    }
    Ok(result)
}

/// Aggregate a slice of byte-string values using the given operation.
fn aggregate_bytes(values: &[Vec<u8>], op: Op) -> Vec<u8> {
    if values.is_empty() {
        return b".".to_vec();
    }
    match op {
        Op::Count => values.len().to_string().into_bytes(),
        Op::CountDistinct => {
            let set: HashSet<&[u8]> = values.iter().map(Vec::as_slice).collect();
            set.len().to_string().into_bytes()
        }
        Op::Collapse => {
            let mut out = Vec::new();
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    out.push(b',');
                }
                out.extend_from_slice(v);
            }
            out
        }
        Op::Distinct => {
            let mut seen: HashSet<&[u8]> = HashSet::new();
            let mut out = Vec::new();
            let mut first = true;
            for v in values {
                if seen.insert(v.as_slice()) {
                    if !first {
                        out.push(b',');
                    }
                    out.extend_from_slice(v);
                    first = false;
                }
            }
            out
        }
        Op::Mode => {
            let mut counts: HashMap<&[u8], usize> = HashMap::new();
            for v in values {
                *counts.entry(v.as_slice()).or_insert(0) += 1;
            }
            let max_count = counts.values().copied().max().unwrap_or(0);
            values
                .iter()
                .find(|v| counts[v.as_slice()] == max_count)
                .cloned()
                .unwrap_or_else(|| b".".to_vec())
        }
        Op::Sum => {
            let s: f64 = values
                .iter()
                .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
                .sum();
            format_f64(s).into_bytes()
        }
        Op::Min => {
            let m = values
                .iter()
                .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
                .reduce(f64::min);
            m.map_or_else(|| b".".to_vec(), |x| format_f64(x).into_bytes())
        }
        Op::Max => {
            let m = values
                .iter()
                .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
                .reduce(f64::max);
            m.map_or_else(|| b".".to_vec(), |x| format_f64(x).into_bytes())
        }
        Op::Mean => {
            let nums: Vec<f64> = values
                .iter()
                .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
                .collect();
            if nums.is_empty() {
                b".".to_vec()
            } else {
                #[allow(clippy::cast_precision_loss)]
                let mean = nums.iter().sum::<f64>() / nums.len() as f64;
                format_f64(mean).into_bytes()
            }
        }
        Op::Median => {
            let mut nums: Vec<f64> = values
                .iter()
                .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
                .collect();
            if nums.is_empty() {
                return b".".to_vec();
            }
            nums.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let mid = nums.len() / 2;
            let median = if nums.len().is_multiple_of(2) {
                (nums[mid - 1] + nums[mid]) / 2.0
            } else {
                nums[mid]
            };
            format_f64(median).into_bytes()
        }
        Op::AbsMin => {
            let m = values
                .iter()
                .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
                .map(f64::abs)
                .reduce(f64::min);
            m.map_or_else(|| b".".to_vec(), |x| format_f64(x).into_bytes())
        }
        Op::AbsMax => {
            let m = values
                .iter()
                .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
                .map(f64::abs)
                .reduce(f64::max);
            m.map_or_else(|| b".".to_vec(), |x| format_f64(x).into_bytes())
        }
    }
}

/// Format a float: integer if no fractional part, else up to 6 significant digits.
fn format_f64(x: f64) -> String {
    if x.fract() == 0.0 && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        let s = format!("{x:.6}");
        let s = s.trim_end_matches('0');
        let s = s.trim_end_matches('.');
        s.to_owned()
    }
}

fn map_reader<R: io::Read>(
    reader: BufReader<R>,
    b_map: &HashMap<String, BChrom>,
    col_ops: &[ColOp],
    null: &str,
    output: &mut dyn Write,
) -> Result<()> {
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    // Reusable per-ColOp buffers of owned byte vecs (cleared between A records).
    let mut hit_values: Vec<Vec<Vec<u8>>> = vec![Vec::new(); col_ops.len()];

    for (lineno_0, line_result) in reader.split(b'\n').enumerate() {
        let line = line_result.map_err(RsomicsError::Io)?;
        if is_skippable(&line) {
            continue;
        }
        let lineno = lineno_0 + 1;

        // Parse A record chrom, start, end.
        let mut tab_pos = [0usize; 3];
        let mut ntabs = 0usize;
        for (i, &b) in line.iter().enumerate() {
            if b == b'\t' {
                if ntabs < 3 {
                    tab_pos[ntabs] = i;
                }
                ntabs += 1;
            }
        }
        if ntabs < 2 {
            return Err(RsomicsError::InvalidInput(format!(
                "A line {lineno}: fewer than 3 fields"
            )));
        }
        let t1 = tab_pos[0];
        let t2 = tab_pos[1];
        let t3 = if ntabs >= 3 { tab_pos[2] } else { line.len() };

        let chrom = std::str::from_utf8(&line[..t1])
            .map_err(|_| RsomicsError::InvalidInput(format!("A line {lineno}: non-UTF8 chrom")))?;
        let a_start: u64 = std::str::from_utf8(&line[t1 + 1..t2])
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let a_end: u64 = std::str::from_utf8(&line[t2 + 1..t3])
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        // Clear reusable hit buffers.
        for vbuf in &mut hit_values {
            vbuf.clear();
        }

        if let Some(bchrom) = b_map.get(chrom) {
            let records = &bchrom.records;
            // Upper bound: B records with start >= a_end cannot overlap.
            let upper = records.partition_point(|r| r.start < a_end);
            // Lower bound: any B record where start + max_feat_size <= a_start cannot
            // have end > a_start, so it cannot overlap. O(log N + K) per query.
            let lower = if bchrom.max_feat_size >= a_start {
                0
            } else {
                records.partition_point(|r| r.start + bchrom.max_feat_size <= a_start)
            };
            for b_rec in &records[lower..upper] {
                if b_rec.end <= a_start {
                    continue;
                }
                for (i, co) in col_ops.iter().enumerate() {
                    if co.op == Op::Count {
                        hit_values[i].push(Vec::new());
                    } else {
                        let bytes = bchrom.field_bytes(b_rec, co.col);
                        hit_values[i].push(bytes.to_vec());
                    }
                }
            }
        }

        // Write original A line + aggregated result columns.
        out.write_all(&line).map_err(RsomicsError::Io)?;
        for (i, co) in col_ops.iter().enumerate() {
            out.write_all(b"\t").map_err(RsomicsError::Io)?;
            if hit_values[i].is_empty() {
                if co.op == Op::Count {
                    out.write_all(b"0").map_err(RsomicsError::Io)?;
                } else {
                    out.write_all(null.as_bytes()).map_err(RsomicsError::Io)?;
                }
            } else {
                let result = aggregate_bytes(&hit_values[i], co.op);
                out.write_all(&result).map_err(RsomicsError::Io)?;
            }
        }
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

/// Run bedtools-map-equivalent on A file vs B file.
///
/// `col_ops` specifies which B column to aggregate and how.
/// `null` is the string to output when no B records overlap (default ".").
pub fn map(
    a_path: &Path,
    b_path: &Path,
    col_ops: &[ColOp],
    null: &str,
    output: &mut dyn Write,
) -> Result<()> {
    let b_map = load_b(b_path)?;
    let a_file = File::open(a_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", a_path.display())))?;
    map_reader(
        BufReader::with_capacity(256 * 1024, a_file),
        &b_map,
        col_ops,
        null,
        output,
    )
}

/// Same as [`map`] but reads A from stdin.
pub fn map_stdin(
    b_path: &Path,
    col_ops: &[ColOp],
    null: &str,
    output: &mut dyn Write,
) -> Result<()> {
    let b_map = load_b(b_path)?;
    map_reader(
        BufReader::with_capacity(256 * 1024, io::stdin()),
        &b_map,
        col_ops,
        null,
        output,
    )
}
