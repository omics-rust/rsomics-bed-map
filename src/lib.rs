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
//! Algorithm: both A and B are loaded fully into memory (like bedtools). B is
//! indexed per-chromosome as a sorted Vec of records with O(1) field access via
//! byte offsets into a raw buffer. For each A record, O(log N) binary search
//! gives the lower and upper bounds of candidate overlaps (using `max_feat_size`
//! to bound the left edge); the K actual overlaps are checked in O(K). Total:
//! O((N+M) + M log M + N log M + NK). Allocation-free inner loop.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// Aggregate operation applied to a set of column values from overlapping B records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    Sum,
    Min,
    Max,
    Mean,
    Count,
    CountDistinct,
    Median,
    Mode,
    Collapse,
    Distinct,
    AbsMin,
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

/// A B record: coordinates + byte offsets of fields into `BChrom::buf`.
#[derive(Debug)]
struct BRecord {
    start: u64,
    end: u64,
    /// Exclusive byte offset of each tab-separated field in `BChrom::buf`.
    field_ends: Vec<u32>,
    /// Byte offset of this line's start in `BChrom::buf`.
    line_start: u32,
}

/// Per-chromosome B data: raw buffer + sorted records + max interval width.
struct BChrom {
    buf: Vec<u8>,
    records: Vec<BRecord>,
    /// Maximum `end - start` across all records — enables lower-bound trimming.
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
            rec.field_ends[col - 2] as usize + 1
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

/// Parse tab positions [0..3) in `line`. Returns count of tabs found.
#[inline]
fn parse_tab3(line: &[u8], tab_pos: &mut [usize; 3]) -> usize {
    let mut ntabs = 0;
    for (i, &b) in line.iter().enumerate() {
        if b == b'\t' {
            if ntabs < 3 {
                tab_pos[ntabs] = i;
            }
            ntabs += 1;
        }
    }
    ntabs
}

/// Load B BED file fully into memory. Field bytes are accessed via offsets into
/// a single per-chromosome buffer — zero per-field String allocations.
fn load_b(path: &Path) -> Result<HashMap<String, BChrom>> {
    let mut raw = Vec::new();
    File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?
        .read_to_end(&mut raw)
        .map_err(RsomicsError::Io)?;

    let mut chrom_data: HashMap<String, (Vec<u8>, Vec<BRecord>, u64)> = HashMap::new();
    let mut tab_pos = [0usize; 3];

    for line in raw.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if is_skippable(line) {
            continue;
        }
        let ntabs = parse_tab3(line, &mut tab_pos);
        if ntabs < 2 {
            continue;
        }
        let t1 = tab_pos[0];
        let t2 = tab_pos[1];
        let t3 = if ntabs >= 3 { tab_pos[2] } else { line.len() };

        let chrom = std::str::from_utf8(&line[..t1])
            .map_err(|_| RsomicsError::InvalidInput("B: non-UTF8 chrom".into()))?;
        let start = fast_parse_u64(&line[t1 + 1..t2]);
        let end = fast_parse_u64(&line[t2 + 1..t3]);

        let (buf, records, max_feat) = chrom_data.entry(chrom.to_owned()).or_default();
        let feat = end.saturating_sub(start);
        if feat > *max_feat {
            *max_feat = feat;
        }
        let line_start = buf.len() as u32;
        buf.extend_from_slice(line);

        let mut field_ends: Vec<u32> = Vec::with_capacity(8);
        let mut pos: u32 = line_start;
        for &b in line {
            pos += 1;
            if b == b'\t' {
                field_ends.push(pos - 1);
            }
        }
        field_ends.push(pos);

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

/// Fast ASCII decimal parser — stops at first non-digit byte.
#[inline]
fn fast_parse_u64(s: &[u8]) -> u64 {
    let mut n = 0u64;
    for &b in s {
        if b < b'0' || b > b'9' {
            break;
        }
        n = n * 10 + u64::from(b - b'0');
    }
    n
}

/// Aggregate byte-string values for the given operation.
fn aggregate_bytes(values: &[&[u8]], op: Op) -> Vec<u8> {
    if values.is_empty() {
        return b".".to_vec();
    }
    match op {
        Op::Count => values.len().to_string().into_bytes(),
        Op::CountDistinct => {
            let set: HashSet<&[u8]> = values.iter().copied().collect();
            set.len().to_string().into_bytes()
        }
        Op::Collapse => {
            let cap = values
                .iter()
                .map(|v| v.len() + 1)
                .sum::<usize>()
                .saturating_sub(1);
            let mut out = Vec::with_capacity(cap);
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
            for &v in values {
                if seen.insert(v) {
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
            for &v in values {
                *counts.entry(v).or_insert(0) += 1;
            }
            let max_count = counts.values().copied().max().unwrap_or(0);
            values
                .iter()
                .find(|&&v| counts[v] == max_count)
                .copied()
                .unwrap_or(b".")
                .to_vec()
        }
        Op::Sum => {
            let s: f64 = values
                .iter()
                .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
                .sum();
            format_f64(s).into_bytes()
        }
        Op::Min => values
            .iter()
            .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
            .reduce(f64::min)
            .map_or_else(|| b".".to_vec(), |x| format_f64(x).into_bytes()),
        Op::Max => values
            .iter()
            .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
            .reduce(f64::max)
            .map_or_else(|| b".".to_vec(), |x| format_f64(x).into_bytes()),
        Op::Mean => {
            let nums: Vec<f64> = values
                .iter()
                .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
                .collect();
            if nums.is_empty() {
                return b".".to_vec();
            }
            #[allow(clippy::cast_precision_loss)]
            let mean = nums.iter().sum::<f64>() / nums.len() as f64;
            format_f64(mean).into_bytes()
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
            let m = if nums.len().is_multiple_of(2) {
                (nums[mid - 1] + nums[mid]) / 2.0
            } else {
                nums[mid]
            };
            format_f64(m).into_bytes()
        }
        Op::AbsMin => values
            .iter()
            .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
            .map(f64::abs)
            .reduce(f64::min)
            .map_or_else(|| b".".to_vec(), |x| format_f64(x).into_bytes()),
        Op::AbsMax => values
            .iter()
            .filter_map(|v| std::str::from_utf8(v).ok()?.parse::<f64>().ok())
            .map(f64::abs)
            .reduce(f64::max)
            .map_or_else(|| b".".to_vec(), |x| format_f64(x).into_bytes()),
    }
}

fn format_f64(x: f64) -> String {
    if x.fract() == 0.0 && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        let s = format!("{x:.6}");
        let s = s.trim_end_matches('0').trim_end_matches('.');
        s.to_owned()
    }
}

/// Process A records from raw bytes, write mapped output.
///
/// Uses `unsafe { transmute }` to extend the lifetime of `&[u8]` slices
/// borrowed from `b_map` (which lives for the whole function) into the
/// `hit_values` buffers. This is sound because:
/// - `b_map` outlives `hit_values` and `out`.
/// - `hit_values[i]` is cleared at the top of every A-record iteration
///   before any new slice is pushed, so no stale borrow from a previous
///   iteration is retained.
/// - The slices are only read by `aggregate_bytes` before the next clear.
#[allow(clippy::transmute_ptr_to_ptr)]
fn map_raw(
    a_raw: &[u8],
    b_map: &HashMap<String, BChrom>,
    col_ops: &[ColOp],
    null: &str,
    output: &mut dyn Write,
) -> Result<()> {
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut hit_values: Vec<Vec<&[u8]>> = vec![Vec::new(); col_ops.len()];
    let mut tab_pos = [0usize; 3];

    for line in a_raw.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if is_skippable(line) {
            continue;
        }
        let ntabs = parse_tab3(line, &mut tab_pos);
        if ntabs < 2 {
            continue;
        }
        let t1 = tab_pos[0];
        let t2 = tab_pos[1];
        let t3 = if ntabs >= 3 { tab_pos[2] } else { line.len() };

        let chrom = std::str::from_utf8(&line[..t1]).unwrap_or("");
        let a_start = fast_parse_u64(&line[t1 + 1..t2]);
        let a_end = fast_parse_u64(&line[t2 + 1..t3]);

        for vbuf in &mut hit_values {
            vbuf.clear();
        }

        if let Some(bchrom) = b_map.get(chrom) {
            let records = &bchrom.records;
            let lower = if bchrom.max_feat_size >= a_start {
                0
            } else {
                records.partition_point(|r| r.start + bchrom.max_feat_size <= a_start)
            };
            let upper = records.partition_point(|r| r.start < a_end);

            for b_rec in &records[lower..upper] {
                if b_rec.end <= a_start {
                    continue;
                }
                for (i, co) in col_ops.iter().enumerate() {
                    let val: &[u8] = if co.op == Op::Count {
                        b""
                    } else {
                        bchrom.field_bytes(b_rec, co.col)
                    };
                    // SAFETY: see function-level comment.
                    let val: &'static [u8] = unsafe { std::mem::transmute(val) };
                    hit_values[i].push(val);
                }
            }
        }

        out.write_all(line).map_err(RsomicsError::Io)?;
        for (i, co) in col_ops.iter().enumerate() {
            out.write_all(b"\t").map_err(RsomicsError::Io)?;
            if hit_values[i].is_empty() {
                let s = if co.op == Op::Count {
                    b"0".as_slice()
                } else {
                    null.as_bytes()
                };
                out.write_all(s).map_err(RsomicsError::Io)?;
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
pub fn map(
    a_path: &Path,
    b_path: &Path,
    col_ops: &[ColOp],
    null: &str,
    output: &mut dyn Write,
) -> Result<()> {
    let b_map = load_b(b_path)?;
    let mut a_raw = Vec::new();
    File::open(a_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", a_path.display())))?
        .read_to_end(&mut a_raw)
        .map_err(RsomicsError::Io)?;
    map_raw(&a_raw, &b_map, col_ops, null, output)
}

/// Same as [`map`] but reads A from stdin.
pub fn map_stdin(
    b_path: &Path,
    col_ops: &[ColOp],
    null: &str,
    output: &mut dyn Write,
) -> Result<()> {
    let b_map = load_b(b_path)?;
    let mut a_raw = Vec::new();
    io::stdin()
        .lock()
        .read_to_end(&mut a_raw)
        .map_err(RsomicsError::Io)?;
    map_raw(&a_raw, &b_map, col_ops, null, output)
}
