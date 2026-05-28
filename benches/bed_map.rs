use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N_RECORDS: usize = 50_000;
const CHROM_SIZE: u64 = 100_000_000;
const SEED: u64 = 0x00BE_DC10_5E51;

fn xorshift(x: &mut u64) -> u64 {
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    *x
}

fn synth_bed_with_score(path: &PathBuf, n: usize, seed: u64) {
    let chroms = ["chr1", "chr2", "chr3"];
    let mut rows: Vec<(String, u64, u64, f64)> = Vec::with_capacity(n);
    let mut rng = seed;
    for _ in 0..n {
        let chrom = chroms[(xorshift(&mut rng) % chroms.len() as u64) as usize];
        let start = xorshift(&mut rng) % (CHROM_SIZE - 1000);
        let end = start + 100 + (xorshift(&mut rng) % 900);
        let score = (xorshift(&mut rng) % 1000) as f64;
        rows.push((chrom.to_string(), start, end, score));
    }
    rows.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let f = File::create(path).expect("create bed");
    let mut w = BufWriter::new(f);
    for (c, s, e, sc) in rows {
        writeln!(w, "{c}\t{s}\t{e}\tfeature\t{sc}").unwrap();
    }
}

fn ensure_fixtures() -> (PathBuf, PathBuf) {
    let mut a = std::env::temp_dir();
    a.push(format!("rsomics-bed-map-bench-a-{N_RECORDS}.bed"));
    let mut b = std::env::temp_dir();
    b.push(format!("rsomics-bed-map-bench-b-{N_RECORDS}.bed"));
    if !a.exists() {
        synth_bed_with_score(&a, N_RECORDS, SEED);
    }
    if !b.exists() {
        synth_bed_with_score(&b, N_RECORDS, SEED ^ 0xDEAD_BEEF);
    }
    (a, b)
}

fn bench(c: &mut Criterion) {
    let (a, b) = ensure_fixtures();
    let ours = env!("CARGO_BIN_EXE_rsomics-bed-map");
    let mut group = c.benchmark_group(format!("bed_map/{N_RECORDS}"));
    group.sample_size(10);

    group.bench_function("rsomics-bed-map", |bm| {
        bm.iter(|| {
            let out = Command::new(ours)
                .arg("-b")
                .arg(&b)
                .arg("-c")
                .arg("5")
                .arg("-o")
                .arg("sum")
                .arg(&a)
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });

    if Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
    {
        group.bench_function("bedtools-map", |bm| {
            bm.iter(|| {
                let out = Command::new("bedtools")
                    .args(["map", "-a"])
                    .arg(&a)
                    .arg("-b")
                    .arg(&b)
                    .args(["-c", "5", "-o", "sum"])
                    .output()
                    .expect("bedtools run");
                assert!(
                    out.status.success(),
                    "{}",
                    String::from_utf8_lossy(&out.stderr)
                );
            });
        });
    } else {
        eprintln!("bedtools not on PATH — skipping upstream comparison");
    }

    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
