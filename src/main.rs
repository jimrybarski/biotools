use anyhow::{bail, Context, Result, anyhow};
use bio::alignment::pairwise::Aligner;
use bio::alignment::{Alignment, AlignmentOperation};
use bio::alphabets::dna::revcomp;
use bio::seq_analysis::gc::gc_content as compute_gc_content;
use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(version, about="Simple bioinformatics tools for sequence analysis and manipulation", long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    #[command(about = "Converts a nucleic acid sequence to its reverse complement.")]
    ReverseComplement {
        #[arg(help = "RNA/DNA sequence")]
        seqs: Vec<String>,
    },
    #[command(about = "Computes the length of a sequence.")]
    Length {
        #[arg(help = "DNA/RNA/protein sequence")]
        seq: Vec<String>,
    },
    #[command(about = "Computes the GC content of a nucleic acid sequence.")]
    GCContent {
        #[arg(help = "RNA/DNA sequence")]
        seqs: Vec<String>,
    },
    #[command(about = "Performs a local pairwise alignment of two sequences.")]
    PairwiseLocal {
        #[arg(help = "DNA/RNA sequence")]
        seqs: Vec<String>,
        #[arg(long, help = "Gap open penalty", default_value_t = 2)]
        gap_open: i32,
        #[arg(long, help = "Gap extend penalty", default_value_t = 1)]
        gap_extend: i32,
        #[arg(long, help = "Hide start/end coordinates of aligned segments")]
        hide_coords: bool,
        #[arg(
            long,
            help = "Aligns the first sequence and its reverse complement with the second sequence and returns the superior alignment."
        )]
        try_rc: bool,
    },
    #[command(about = "Performs a semiglobal pairwise alignment of two sequences.")]
    PairwiseSemiglobal {
        #[arg(help = "DNA/RNA sequence")]
        seqs: Vec<String>,
        #[arg(long, help = "Gap open penalty", default_value_t = 2)]
        gap_open: i32,
        #[arg(long, help = "Gap extend penalty", default_value_t = 1)]
        gap_extend: i32,
        #[arg(long, help = "Hide start/end coordinates of aligned segments")]
        hide_coords: bool,
        #[arg(
            long,
            help = "Aligns the first sequence and its reverse complement with the second sequence and returns the superior alignment."
        )]
        try_rc: bool,
    },
    #[command(about = "Performs a global pairwise alignment of two sequences.")]
    PairwiseGlobal {
        #[arg(help = "DNA/RNA sequence")]
        seqs: Vec<String>,
        #[arg(long, help = "Gap open penalty", default_value_t = 2)]
        gap_open: i32,
        #[arg(long, help = "Gap extend penalty", default_value_t = 1)]
        gap_extend: i32,
        #[arg(long, help = "Hide start/end coordinates of aligned segments")]
        hide_coords: bool,
        #[arg(
            long,
            help = "Aligns the first sequence and its reverse complement with the second sequence and returns the superior alignment."
        )]
        try_rc: bool,
    },
}

extern crate bio;

enum AlignmentCommand {
    Local,
    Global,
    Semiglobal,
}

fn build_reverse_complement(user_input: Vec<String>) -> Result<String> {
    let reversed_complements: Vec<String> = user_input
        .into_iter()
        .rev()
        .map(|sequence| {
            let seq = sequence.into_bytes();
            String::from_utf8(revcomp(seq)).context("Failed to build reverse complement")
        })
        .collect::<Result<Vec<_>, anyhow::Error>>()?;
    Ok(reversed_complements.join(" "))
}

fn get_string_length(seqs: Vec<String>) -> Result<String> {
    Ok(seqs.join("")
        .chars()
        .filter(|ch| *ch != '-')
        .filter(|ch| *ch != ' ')
        .collect::<Vec<_>>()
        .len()
        .to_string())
}

fn confirm_valid_nucleic_acid(seq: &str) -> Result<()> {
    for (i, c) in seq.chars().enumerate() {
        if !matches!(c, 'A' | 'C' | 'G' | 'T' | 'U' | 'a' | 'c' | 'g' | 't' | 'u') {
            return Err(anyhow!("Invalid/ambiguous base: '{c}' at position {i}"));
        }
    }
    Ok(())
}

fn gc_content(seqs: Vec<String>) -> Result<String> {
    let seq = seqs.join("").replace(" ", "").replace("-", "");
    confirm_valid_nucleic_acid(&seq)?;
    let gc = compute_gc_content(seq.as_bytes());
    Ok(format!("{:.16}", gc))
}

fn format_alignment(alignment: Alignment, a: String, b: String, hide_coords: bool, comment: String) -> String {
    let (a_raw_start, b_raw_start) = if hide_coords {
        ("".to_string(), "".to_string())
    } else {
        (alignment.xstart.to_string(), alignment.ystart.to_string())
    };

    let (a_raw_end, b_raw_end) = if hide_coords {
        ("".to_string(), "".to_string())
    } else {
        (alignment.xend.to_string(), alignment.yend.to_string())
    };

    let start_length = std::cmp::max(a_raw_start.len(), b_raw_start.len());
    let end_length = std::cmp::max(a_raw_end.len(), b_raw_end.len());

    let (mut a_pretty, mut alignment_string, mut b_pretty, a_end, b_end) = if hide_coords {
        (
            "".to_string(),
            "".to_string(),
            "".to_string(),
            "".to_string(),
            "".to_string(),
        )
    } else {
        (
            format!("{:>width$} ", a_raw_start, width = start_length),
            format!("{:>width$} ", "".to_string(), width = start_length),
            format!("{:>width$} ", b_raw_start, width = start_length),
            format!("{:>width$}", a_raw_end, width = end_length),
            format!("{:>width$}", b_raw_end, width = end_length),
        )
    };

    let mut a_index = alignment.xstart;
    let mut b_index = alignment.ystart;

    for op in alignment.operations {
        match op {
            AlignmentOperation::Match => {
                let c = &a[a_index..a_index + 1];
                a_pretty.push_str(c);
                a_index += 1;
                alignment_string.push('|');
                let d = &b[b_index..b_index + 1];
                b_pretty.push_str(d);
                b_index += 1;
            }
            AlignmentOperation::Del => {
                a_pretty.push('-');
                alignment_string.push(' ');
                let d = &b[b_index..b_index + 1];
                b_pretty.push_str(d);
                b_index += 1;
            }
            AlignmentOperation::Ins => {
                let c = &a[a_index..a_index + 1];
                a_pretty.push_str(c);
                a_index += 1;
                alignment_string.push(' ');
                b_pretty.push('-');
            }
            AlignmentOperation::Subst => {
                let c = &a[a_index..a_index + 1];
                a_pretty.push_str(c);
                a_index += 1;
                alignment_string.push('.');
                let d = &b[b_index..b_index + 1];
                b_pretty.push_str(d);
                b_index += 1;
            }
            AlignmentOperation::Xclip(n) => {
                for _ in 0..n {
                    a_pretty.push('-');
                    b_pretty.push(' ');
                    b_pretty.push(' ');
                }
            }
            AlignmentOperation::Yclip(n) => {
                for _ in 0..n {
                    a_pretty.push(' ');
                    b_pretty.push(' ');
                    b_pretty.push('-');
                }
            }
        }
    }
    format!("{a_pretty} {a_end}{comment}\n{alignment_string}\n{b_pretty} {b_end}")
}

fn run_alignment(
    alignment_command: &AlignmentCommand,
    aligner: &mut Aligner<impl Fn(u8, u8) -> i32>,
    a_bytes: &[u8],
    b_bytes: &[u8],
) -> Alignment {
    match alignment_command {
        AlignmentCommand::Local => aligner.local(a_bytes, b_bytes),
        AlignmentCommand::Semiglobal => aligner.semiglobal(a_bytes, b_bytes),
        AlignmentCommand::Global => aligner.global(a_bytes, b_bytes),
    }
}

fn pairwise(
    alignment_command: AlignmentCommand,
    seqs: Vec<String>,
    gap_open_score: i32,
    gap_extend_score: i32,
    hide_coords: bool,
    try_rc: bool,
) -> Result<String> {
    if seqs.len() != 2 {
        bail!("Pairwise comparison needs exactly two sequences");
    }
    let gap_open_score = -gap_open_score;
    let gap_extend_score = -gap_extend_score;

    let a = seqs[0].clone();
    let b = seqs[1].clone();
    let a_bytes = a.as_bytes();
    let b_bytes = b.as_bytes();
    let score = |a: u8, b: u8| if a.to_ascii_uppercase() == b.to_ascii_uppercase() { 1i32 } else { -1i32 };

    let mut aligner =
        Aligner::with_capacity(a.len(), b.len(), gap_open_score, gap_extend_score, &score);
    let alignment = run_alignment(&alignment_command, &mut aligner, a_bytes, b_bytes);
    let (alignment, a, comment) = if try_rc {
        let a_rc_bytes = revcomp(a.as_bytes());
        let a_rc = String::from_utf8(a_rc_bytes.clone())?;
        let alignment_rc = run_alignment(&alignment_command, &mut aligner, &a_rc_bytes, b_bytes);
        if alignment.score >= alignment_rc.score {
            (alignment, a, "".to_string())
        } else {
            (alignment_rc, a_rc, " RC".to_string())
        }
    } else {
        (alignment, a, "".to_string())
    };

    let pretty_alignment = format_alignment(alignment, a, b, hide_coords, comment);
    Ok(pretty_alignment)
}

fn abort(error_message: &str) -> ! {
    eprintln!("biotools error: {}", error_message);
    std::process::exit(47)
}

fn main() -> Result<()> {
    let args = Args::parse();

    let output = match args.command {
        Commands::ReverseComplement { seqs } => build_reverse_complement(seqs),
        Commands::Length { seq } => get_string_length(seq),
        Commands::GCContent { seqs } => gc_content(seqs),
        Commands::PairwiseLocal {
            seqs,
            gap_open,
            gap_extend,
            hide_coords,
            try_rc,
        } => pairwise(
            AlignmentCommand::Local,
            seqs,
            gap_open,
            gap_extend,
            hide_coords,
            try_rc,
        ),
        Commands::PairwiseSemiglobal {
            seqs,
            gap_open,
            gap_extend,
            hide_coords,
            try_rc,
        } => pairwise(
            AlignmentCommand::Semiglobal,
            seqs,
            gap_open,
            gap_extend,
            hide_coords,
            try_rc,
        ),
        Commands::PairwiseGlobal {
            seqs,
            gap_open,
            gap_extend,
            hide_coords,
            try_rc,
        } => pairwise(
            AlignmentCommand::Global,
            seqs,
            gap_open,
            gap_extend,
            hide_coords,
            try_rc,
        ),
    };

    match output {
        Ok(text) => {
            println!("{}", text);
            Ok(())
        }
        Err(e) => abort(&format!("Error: {:?}", e)),
    }
}
