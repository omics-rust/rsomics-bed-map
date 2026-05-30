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
//! Algorithm: O((N+M) + M log M + N log M + NK). Both files fully loaded;
//! B indexed per-chromosome as a sorted Vec with O(1) field access via byte
//! offsets into a raw buffer. `max_feat_size` bounds the binary-search left
//! edge. Allocation-free inner loop.

mod driver;
mod index;
mod ops;

pub use driver::{map, map_stdin};
pub use ops::{ColOp, Op};
