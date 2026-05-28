use std::path::Path;
use std::process::Command;

use rsomics_bed_map::{ColOp, Op, map};

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn sum_col5() {
    // regionA [0,100): overlaps feat1 [10,50) score=10 and feat2 [40,90) score=20 → sum=30
    // regionB [200,400): overlaps feat3 [250,350) score=30 and feat4 [300,450) score=40 → sum=70
    // regionC chr2 [0,500): overlaps feat5 [100,200) score=50 and feat6 [300,400) score=60 → sum=110
    let a = golden("a.bed");
    let b = golden("b.bed");
    let col_ops = vec![ColOp {
        col: 5,
        op: Op::Sum,
    }];
    let mut out = Vec::new();
    map(&a, &b, &col_ops, ".", &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = result.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(lines.len(), 3, "expected 3 output lines: {result}");
    let cols0: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(
        cols0.last().copied(),
        Some("30"),
        "regionA sum should be 30: {}",
        lines[0]
    );
    let cols1: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(
        cols1.last().copied(),
        Some("70"),
        "regionB sum should be 70: {}",
        lines[1]
    );
    let cols2: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(
        cols2.last().copied(),
        Some("110"),
        "regionC sum should be 110: {}",
        lines[2]
    );
}

#[test]
fn count_op() {
    let a = golden("a.bed");
    let b = golden("b.bed");
    let col_ops = vec![ColOp {
        col: 4,
        op: Op::Count,
    }];
    let mut out = Vec::new();
    map(&a, &b, &col_ops, ".", &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = result.lines().filter(|l| !l.is_empty()).collect();
    let counts: Vec<&str> = lines
        .iter()
        .map(|l| l.split('\t').next_back().unwrap_or("."))
        .collect();
    assert_eq!(counts, vec!["2", "2", "2"], "counts: {result}");
}

#[test]
fn no_overlap_null() {
    use std::io::Write;
    use tempfile::NamedTempFile;
    let mut fa = NamedTempFile::new().unwrap();
    let mut fb = NamedTempFile::new().unwrap();
    writeln!(fa, "chr1\t0\t100\tregion").unwrap();
    writeln!(fb, "chr1\t200\t300\tfeat\t99").unwrap();
    let col_ops = vec![ColOp {
        col: 5,
        op: Op::Sum,
    }];
    let mut out = Vec::new();
    map(fa.path(), fb.path(), &col_ops, ".", &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    assert!(
        result.trim_end().ends_with('.'),
        "null should be '.': {result}"
    );
}

#[test]
fn mean_op() {
    let a = golden("a.bed");
    let b = golden("b.bed");
    let col_ops = vec![ColOp {
        col: 5,
        op: Op::Mean,
    }];
    let mut out = Vec::new();
    map(&a, &b, &col_ops, ".", &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = result.lines().filter(|l| !l.is_empty()).collect();
    // regionA: mean(10, 20) = 15
    let cols0: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(
        cols0.last().copied(),
        Some("15"),
        "regionA mean should be 15: {}",
        lines[0]
    );
}

#[test]
fn bedtools_compat() {
    let bedtools = Command::new("bedtools").arg("--version").output();
    if bedtools.is_err() || !bedtools.unwrap().status.success() {
        eprintln!("bedtools not available — skipping compat test");
        return;
    }

    let a = golden("a.bed");
    let b = golden("b.bed");

    // Test sum of col 5 (bedtools uses -c 5 -o sum).
    let col_ops = vec![ColOp {
        col: 5,
        op: Op::Sum,
    }];
    let mut ours = Vec::new();
    map(&a, &b, &col_ops, ".", &mut ours).unwrap();
    let ours_str = String::from_utf8(ours).unwrap();

    let bt = Command::new("bedtools")
        .args(["map", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .args(["-c", "5", "-o", "sum"])
        .output()
        .expect("bedtools map failed");
    let bt_str = String::from_utf8(bt.stdout).unwrap();

    let mut ours_lines: Vec<&str> = ours_str.lines().filter(|l| !l.is_empty()).collect();
    let mut bt_lines: Vec<&str> = bt_str.lines().filter(|l| !l.is_empty()).collect();
    ours_lines.sort_unstable();
    bt_lines.sort_unstable();

    assert_eq!(
        ours_lines, bt_lines,
        "output differs from bedtools map -c 5 -o sum"
    );
}
