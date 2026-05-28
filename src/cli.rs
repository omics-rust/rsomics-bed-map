use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_map::{ColOp, Op, map, map_stdin};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-map", disable_help_flag = true)]
pub struct Cli {
    /// Input BED file A (default: stdin)
    input: Option<PathBuf>,
    /// B BED file to map over A intervals
    #[arg(short = 'b', long)]
    b: PathBuf,
    /// Comma-separated 1-based column index(es) in B to aggregate (e.g. 4 or 4,5)
    #[arg(short = 'c', long, default_value = "5")]
    columns: String,
    /// Comma-separated operation(s) to apply: sum,min,max,mean,count,count_distinct,
    /// median,mode,collapse,distinct,absmin,absmax
    #[arg(short = 'o', long, default_value = "sum")]
    operations: String,
    /// Null/no-hit output string (default: ".")
    #[arg(short = 'n', long, default_value = ".")]
    null: String,
    /// Output path (default: stdout)
    #[arg(long = "out")]
    output: Option<PathBuf>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let cols: Vec<usize> = self
            .columns
            .split(',')
            .map(|s| {
                s.trim()
                    .parse::<usize>()
                    .map_err(|_| RsomicsError::InvalidInput(format!("invalid column index: {s:?}")))
            })
            .collect::<Result<Vec<_>>>()?;

        let ops: Vec<Op> = self
            .operations
            .split(',')
            .map(|s| {
                Op::parse_op(s.trim())
                    .ok_or_else(|| RsomicsError::InvalidInput(format!("unknown operation: {s:?}")))
            })
            .collect::<Result<Vec<_>>>()?;

        if cols.len() != ops.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "number of columns ({}) must match number of operations ({})",
                cols.len(),
                ops.len()
            )));
        }

        let col_ops: Vec<ColOp> = cols
            .into_iter()
            .zip(ops)
            .map(|(col, op)| ColOp { col, op })
            .collect();

        let mut stdout_lock;
        let mut file_out;
        let out: &mut dyn io::Write = if let Some(ref p) = self.output {
            file_out = std::fs::File::create(p).map_err(RsomicsError::Io)?;
            &mut file_out
        } else {
            stdout_lock = io::stdout().lock();
            &mut stdout_lock
        };

        match self.input {
            Some(ref p) => map(p.as_path(), &self.b, &col_ops, &self.null, out),
            None => map_stdin(&self.b, &col_ops, &self.null, out),
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Aggregate column values from B intervals overlapping each A interval (bedtools map equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["-b <B> [OPTIONS] [INPUT]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('b'),
                long: "b",
                aliases: &[],
                value: Some("<FILE>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "B BED file to map over A intervals",
                why_default: None,
            },
            FlagSpec {
                short: Some('c'),
                long: "columns",
                aliases: &[],
                value: Some("<INT[,INT...]>"),
                type_hint: Some("String"),
                required: false,
                default: Some("5"),
                description: "Comma-separated 1-based column index(es) in B to aggregate",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "operations",
                aliases: &[],
                value: Some("<OP[,OP...]>"),
                type_hint: Some("String"),
                required: false,
                default: Some("sum"),
                description: "Comma-separated operations: sum,min,max,mean,count,count_distinct,median,mode,collapse,distinct,absmin,absmax",
                why_default: None,
            },
            FlagSpec {
                short: Some('n'),
                long: "null",
                aliases: &[],
                value: Some("<STR>"),
                type_hint: Some("String"),
                required: false,
                default: Some("."),
                description: "Output string when no B records overlap A",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "out",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdout"),
                description: "Output path",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Sum scores of B peaks overlapping each A region",
            command: "rsomics-bed-map -b peaks.bed -c 5 -o sum regions.bed",
        },
        Example {
            description: "Report count and mean of B column 4 per A interval",
            command: "rsomics-bed-map -b features.bed -c 4,4 -o count,mean regions.bed",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
