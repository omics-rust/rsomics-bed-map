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
//! Algorithm: B records loaded into per-chromosome sorted Vecs; for each A
//! record the overlap range is binary-searched for O(log N + K) queries.

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

#[derive(Debug, Clone)]
struct BRecord {
    start: u64,
    end: u64,
    /// All fields of the B line split by tab.
    fields: Vec<String>,
}

fn is_skippable(line: &str) -> bool {
    line.is_empty()
        || line.starts_with('#')
        || line.starts_with("track")
        || line.starts_with("browser")
}

/// Load B BED file into per-chromosome Vecs sorted by start.
fn load_b(path: &Path) -> Result<HashMap<String, Vec<BRecord>>> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut map: HashMap<String, Vec<BRecord>> = HashMap::new();
    for line in raw.lines() {
        if is_skippable(line) {
            continue;
        }
        let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if fields.len() < 3 {
            continue;
        }
        let chrom = fields[0].clone();
        let start: u64 = fields[1].parse().unwrap_or(0);
        let end: u64 = fields[2].parse().unwrap_or(0);
        map.entry(chrom)
            .or_default()
            .push(BRecord { start, end, fields });
    }
    for v in map.values_mut() {
        v.sort_unstable_by_key(|r| r.start);
    }
    Ok(map)
}

/// Get the value of a 1-based column from a BRecord's field list.
/// Returns empty string if column is out of range.
fn field_value(rec: &BRecord, col: usize) -> &str {
    if col == 0 || col > rec.fields.len() {
        ""
    } else {
        &rec.fields[col - 1]
    }
}

/// Aggregate a slice of string values using the given operation.
fn aggregate(values: &[&str], op: Op) -> String {
    if values.is_empty() {
        return ".".to_owned();
    }
    match op {
        Op::Count => values.len().to_string(),
        Op::CountDistinct => {
            let set: HashSet<&str> = values.iter().copied().collect();
            set.len().to_string()
        }
        Op::Collapse => values.join(","),
        Op::Distinct => {
            let mut seen = HashSet::new();
            let mut out = Vec::new();
            for v in values {
                if seen.insert(*v) {
                    out.push(*v);
                }
            }
            out.join(",")
        }
        Op::Mode => {
            let mut counts: HashMap<&str, usize> = HashMap::new();
            for v in values {
                *counts.entry(v).or_insert(0) += 1;
            }
            let max_count = counts.values().copied().max().unwrap_or(0);
            values
                .iter()
                .find(|&&v| counts[v] == max_count)
                .copied()
                .unwrap_or(".")
                .to_owned()
        }
        Op::Sum => {
            let s: f64 = values.iter().filter_map(|v| v.parse::<f64>().ok()).sum();
            format_f64(s)
        }
        Op::Min => {
            let m = values
                .iter()
                .filter_map(|v| v.parse::<f64>().ok())
                .reduce(f64::min);
            m.map(format_f64).unwrap_or(".".to_owned())
        }
        Op::Max => {
            let m = values
                .iter()
                .filter_map(|v| v.parse::<f64>().ok())
                .reduce(f64::max);
            m.map(format_f64).unwrap_or(".".to_owned())
        }
        Op::Mean => {
            let nums: Vec<f64> = values
                .iter()
                .filter_map(|v| v.parse::<f64>().ok())
                .collect();
            if nums.is_empty() {
                ".".to_owned()
            } else {
                format_f64(nums.iter().sum::<f64>() / nums.len() as f64)
            }
        }
        Op::Median => {
            let mut nums: Vec<f64> = values
                .iter()
                .filter_map(|v| v.parse::<f64>().ok())
                .collect();
            if nums.is_empty() {
                return ".".to_owned();
            }
            nums.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let mid = nums.len() / 2;
            let median = if nums.len().is_multiple_of(2) {
                (nums[mid - 1] + nums[mid]) / 2.0
            } else {
                nums[mid]
            };
            format_f64(median)
        }
        Op::AbsMin => {
            let m = values
                .iter()
                .filter_map(|v| v.parse::<f64>().ok())
                .map(|x| x.abs())
                .reduce(f64::min);
            m.map(format_f64).unwrap_or(".".to_owned())
        }
        Op::AbsMax => {
            let m = values
                .iter()
                .filter_map(|v| v.parse::<f64>().ok())
                .map(|x| x.abs())
                .reduce(f64::max);
            m.map(format_f64).unwrap_or(".".to_owned())
        }
    }
}

/// Format a float: integer if no fractional part, else up to 6 significant digits.
fn format_f64(x: f64) -> String {
    if x.fract() == 0.0 && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        // Trim trailing zeros after decimal point.
        let s = format!("{x:.6}");
        let s = s.trim_end_matches('0');
        let s = s.trim_end_matches('.');
        s.to_owned()
    }
}

fn map_reader<R: io::Read>(
    reader: BufReader<R>,
    b_map: &HashMap<String, Vec<BRecord>>,
    col_ops: &[ColOp],
    null: &str,
    output: &mut dyn Write,
) -> Result<()> {
    let mut out = BufWriter::new(output);
    let mut hit_values: Vec<Vec<&str>> = vec![Vec::new(); col_ops.len()];

    for (lineno_0, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        if is_skippable(&line) {
            continue;
        }
        let lineno = lineno_0 + 1;
        let mut fields = line.splitn(4, '\t');
        let chrom = fields
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput(format!("line {lineno}: missing chrom")))?
            .to_owned();
        let start: u64 = fields
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput(format!("line {lineno}: missing start")))?
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("line {lineno}: bad start")))?;
        let end: u64 = fields
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput(format!("line {lineno}: missing end")))?
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("line {lineno}: bad end")))?;

        // Clear hit value buffers.
        for vbuf in hit_values.iter_mut() {
            vbuf.clear();
        }

        if let Some(b_vec) = b_map.get(&chrom) {
            // Binary search for upper bound: all B with start < end.
            let limit = b_vec.partition_point(|r| r.start < end);
            for b in &b_vec[..limit] {
                if b.end <= start {
                    continue;
                }
                // B overlaps [start, end).
                for (i, co) in col_ops.iter().enumerate() {
                    // Count op doesn't need the value, but we still count hits.
                    if co.op == Op::Count {
                        hit_values[i].push("");
                    } else {
                        hit_values[i].push(field_value(b, co.col));
                    }
                }
            }
        }

        // Write A line (all columns) + aggregated columns.
        out.write_all(line.as_bytes()).map_err(RsomicsError::Io)?;
        for (i, co) in col_ops.iter().enumerate() {
            let result = if hit_values[i].is_empty() {
                if co.op == Op::Count {
                    "0".to_owned()
                } else {
                    null.to_owned()
                }
            } else {
                aggregate(&hit_values[i], co.op)
            };
            write!(out, "\t{result}").map_err(RsomicsError::Io)?;
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
    map_reader(BufReader::new(a_file), &b_map, col_ops, null, output)
}

/// Same as [`map`] but reads A from stdin.
pub fn map_stdin(
    b_path: &Path,
    col_ops: &[ColOp],
    null: &str,
    output: &mut dyn Write,
) -> Result<()> {
    let b_map = load_b(b_path)?;
    map_reader(BufReader::new(io::stdin()), &b_map, col_ops, null, output)
}
