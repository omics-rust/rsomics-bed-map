use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub(crate) struct BRecord {
    pub(crate) start: u64,
    pub(crate) end: u64,
    /// Exclusive end byte offset of each tab-separated field in `BChrom::buf`.
    pub(crate) field_ends: Vec<u32>,
    pub(crate) line_start: u32,
}

pub(crate) struct BChrom {
    pub(crate) buf: Vec<u8>,
    pub(crate) records: Vec<BRecord>,
    /// Maximum `end - start` across all records — bounds the binary-search left edge.
    pub(crate) max_feat_size: u64,
}

impl BChrom {
    /// Bytes of 1-based column `col` from `rec`.
    pub(crate) fn field_bytes<'a>(&'a self, rec: &'a BRecord, col: usize) -> &'a [u8] {
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
pub(crate) fn is_skippable(line: &[u8]) -> bool {
    line.is_empty()
        || line.starts_with(b"#")
        || line.starts_with(b"track")
        || line.starts_with(b"browser")
}

#[inline]
pub(crate) fn parse_tab3(line: &[u8], tab_pos: &mut [usize; 3]) -> usize {
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

#[inline]
pub(crate) fn fast_parse_u64(s: &[u8]) -> u64 {
    let mut n = 0u64;
    for &b in s {
        if !b.is_ascii_digit() {
            break;
        }
        n = n * 10 + u64::from(b - b'0');
    }
    n
}

/// Field bytes accessed via offsets into a per-chrom buffer — no per-field allocation.
pub(crate) fn load_b(path: &Path) -> Result<HashMap<String, BChrom>> {
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
