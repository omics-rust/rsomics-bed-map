use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::index::{BChrom, fast_parse_u64, is_skippable, load_b, parse_tab3};
use crate::ops::{ColOp, Op, aggregate_bytes};

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
pub(crate) fn map_raw(
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
