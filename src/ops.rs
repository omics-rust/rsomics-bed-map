use std::collections::{HashMap, HashSet};

/// Aggregate operation applied to overlapping B column values.
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

/// One output column: which B field to read and how to aggregate it.
#[derive(Debug, Clone)]
pub struct ColOp {
    /// 1-based column index into the B record.
    pub col: usize,
    pub op: Op,
}

/// Aggregate byte-string values for the given operation. Returns `.` when empty.
pub(crate) fn aggregate_bytes(values: &[&[u8]], op: Op) -> Vec<u8> {
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

pub(crate) fn format_f64(x: f64) -> String {
    if x.fract() == 0.0 && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        let s = format!("{x:.6}");
        let s = s.trim_end_matches('0').trim_end_matches('.');
        s.to_owned()
    }
}
