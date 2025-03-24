use anyhow::{anyhow, bail, Context, Result};
use bio::alignment::pairwise::Aligner;
use bio::alignment::{Alignment, AlignmentOperation};
use bio::alphabets::dna::revcomp;
use bio::seq_analysis::gc::gc_content as rustbio_gc_content;
use clap::{Parser, Subcommand};
use std::cmp;

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
        #[arg(long, help = "Maximum width of aligned characters", default_value_t = 60)]
        line_width: usize,
        #[arg(long, help = "Use zero-based coordinates")]
        use_0_based_coords: bool,
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
        #[arg(long, help = "Maximum width of aligned characters", default_value_t = 60)]
        line_width: usize,
        #[arg(long, help = "Use zero-based coordinates")]
        use_0_based_coords: bool,
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
        #[arg(long, help = "Maximum width of aligned characters", default_value_t = 60)]
        line_width: usize,
        #[arg(long, help = "Use zero-based coordinates")]
        use_0_based_coords: bool,
    },
}

extern crate bio;

enum AlignmentCommand {
    Local,
    Global,
    Semiglobal,
}

struct DisplayOptions {
    hide_coords: bool,
    try_rc: bool,
    line_width: usize,
    use_0_based_coords: bool,
}

fn build_reverse_complement(seqs: Vec<String>) -> Result<String> {
    let reversed_complements: Vec<String> = seqs
        .into_iter()
        .rev()
        .map(|sequence| {
            let seq = sequence.into_bytes();
            String::from_utf8(revcomp(seq)).context("Failed to build reverse complement")
        })
        .collect::<Result<Vec<_>, anyhow::Error>>()?;
    Ok(reversed_complements.join(" "))
}

fn get_seq_length(seqs: Vec<String>) -> Result<String> {
    Ok(seqs
        .join("")
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

fn compute_gc_content(seqs: Vec<String>) -> Result<f32> {
    let seq = seqs.join("").replace(" ", "").replace("-", "");
    confirm_valid_nucleic_acid(&seq)?;
    Ok(rustbio_gc_content(seq.as_bytes()))
}

fn gc_content(seqs: Vec<String>) -> Result<String> {
    let gc = compute_gc_content(seqs)?;
    Ok(format!("{:.16}", gc))
}

struct AlignmentDisplayLine {
    a_alignment: String,
    b_alignment: String,
    a_start: usize,
    a_end: usize,
    b_start: usize,
    b_end: usize,
    alignment_string: String,
}

fn make_display_lines(
    alignment: Alignment,
    a_seq: String,
    b_seq: String,
    line_width: usize,
) -> (Vec<AlignmentDisplayLine>, usize) {
    let mut a_start = alignment.xstart;
    let mut b_start = alignment.ystart;
    let mut a_index = alignment.xstart;
    let mut b_index = alignment.ystart;
    let mut display_lines: Vec<AlignmentDisplayLine> = vec![];

    for op_chunk in alignment.operations.chunks(line_width) {
        let mut a_alignment = "".to_string();
        let mut b_alignment = "".to_string();
        let mut alignment_string = "".to_string();

        for op in op_chunk {
            match op {
                AlignmentOperation::Match => {
                    let a_char = &a_seq[a_index..a_index + 1];
                    a_alignment.push_str(a_char);
                    a_index += 1;

                    let b_char = &b_seq[b_index..b_index + 1];
                    b_alignment.push_str(b_char);
                    b_index += 1;

                    alignment_string.push('|');
                }
                AlignmentOperation::Del => {
                    a_alignment.push('-');
                    let b_char = &b_seq[b_index..b_index + 1];
                    b_alignment.push_str(b_char);
                    b_index += 1;

                    alignment_string.push(' ');
                }
                AlignmentOperation::Ins => {
                    let a_char = &a_seq[a_index..a_index + 1];
                    a_alignment.push_str(a_char);
                    a_index += 1;

                    b_alignment.push('-');

                    alignment_string.push(' ');
                }
                AlignmentOperation::Subst => {
                    let a_char = &a_seq[a_index..a_index + 1];
                    a_alignment.push_str(a_char);
                    a_index += 1;

                    let b_char = &b_seq[b_index..b_index + 1];
                    b_alignment.push_str(b_char);
                    b_index += 1;

                    alignment_string.push('.');
                }
                AlignmentOperation::Xclip(n) => {
                    for _ in 0..*n {
                        a_alignment.push('-');
                        b_alignment.push(' ');
                        alignment_string.push(' ');
                    }
                }
                AlignmentOperation::Yclip(n) => {
                    for _ in 0..*n {
                        a_alignment.push(' ');
                        b_alignment.push('-');
                        alignment_string.push(' ');
                    }
                }
            }
        }
        let display_line = AlignmentDisplayLine {
            a_alignment,
            b_alignment,
            a_start,
            a_end: a_index,
            b_start,
            b_end: b_index,
            alignment_string,
        };
        a_start = a_index;
        b_start = b_index;
        display_lines.push(display_line);
    }
    (display_lines, a_index)
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
    opts: DisplayOptions,
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
    let score = |a: u8, b: u8| {
        if a.to_ascii_uppercase() == b.to_ascii_uppercase() {
            1i32
        } else {
            -1i32
        }
    };

    let mut aligner =
        Aligner::with_capacity(a.len(), b.len(), gap_open_score, gap_extend_score, &score);
    let alignment = run_alignment(&alignment_command, &mut aligner, a_bytes, b_bytes);
    let (alignment, a, a_is_rc) = if opts.try_rc {
        let a_rc_bytes = revcomp(a.as_bytes());
        let a_rc = String::from_utf8(a_rc_bytes.clone())?;
        let alignment_rc = run_alignment(&alignment_command, &mut aligner, &a_rc_bytes, b_bytes);
        if alignment.score >= alignment_rc.score {
            (alignment, a, false)
        } else {
            (alignment_rc, a_rc, true)
        }
    } else {
        (alignment, a, false)
    };

    let (display_lines, a_end) = make_display_lines(alignment, a, b, opts.line_width);
    let pretty_alignment = format_display_lines(
        &display_lines,
        opts.hide_coords,
        a_end,
        a_is_rc,
        opts.use_0_based_coords,
    );
    Ok(pretty_alignment)
}

fn format_display_lines(
    display_lines: &[AlignmentDisplayLine],
    hide_coords: bool,
    final_a_end: usize,
    a_is_rc: bool,
    use_0_based_coordinates: bool,
) -> String {
    let mut output: Vec<String> = vec![];

    // If we have numbers that span an order of magnitude, they won't have the same width (e.g.
    // "10" occupies two colums while "100" occupies three). We need to pad the coordinates on
    // the left of the alignments with spaces so that each alignment string is lined up with
    // all the others. Thus, we need to find the widest number first.
    let mut max_coordinate = 0;
    for line in display_lines {
        let (a_start, _, b_start) =
            calculate_formatted_coordinates(line, final_a_end, a_is_rc, use_0_based_coordinates);

        max_coordinate = cmp::max(
            max_coordinate,
            [a_start, b_start].into_iter().max().unwrap(),
        );
    }
    let width = max_coordinate.to_string().len();

    for line in display_lines {
        let mut line_output = "".to_string();

        if hide_coords {
            line_output.push_str(&line.a_alignment);
            line_output.push('\n');
            line_output.push_str(&line.alignment_string);
            line_output.push('\n');
            line_output.push_str(&line.b_alignment);
            output.push(line_output);
        } else {
            let (a_start, a_end, b_start) = calculate_formatted_coordinates(
                line,
                final_a_end,
                a_is_rc,
                use_0_based_coordinates,
            );

            let a_start_text = format!("{:>width$}", a_start, width = width);
            let b_start_text = format!("{:>width$}", b_start, width = width);
            let alignment_string_start_text = format!("{:>width$}", "".to_string(), width = width);

            line_output
                .push_str(format!("{} {} {}\n", a_start_text, line.a_alignment, a_end).as_str());
            line_output.push_str(
                format!(
                    "{} {}\n",
                    alignment_string_start_text, line.alignment_string
                )
                .as_str(),
            );
            line_output
                .push_str(format!("{} {} {}", b_start_text, line.b_alignment, line.b_end).as_str());
            output.push(line_output);
        }
    }
    output.join("\n\n")
}

fn calculate_formatted_coordinates(
    line: &AlignmentDisplayLine,
    final_a_end: usize,
    a_is_rc: bool,
    use_0_based_coordinates: bool,
) -> (usize, usize, usize) {
    let mut a_start = if a_is_rc {
        final_a_end - line.a_start
    } else {
        line.a_start
    };

    let mut a_end = if a_is_rc {
        final_a_end - line.a_end + 1
    } else {
        line.a_end
    };

    let mut b_start = line.b_start;
    if !use_0_based_coordinates {
        b_start += 1;
        if !a_is_rc {
            a_start += 1;
        }
    } else if a_is_rc {
        a_end -= 1;
    }
    (a_start, a_end, b_start)
}
fn abort(error_message: &str) -> ! {
    eprintln!("biotools error: {}", error_message);
    std::process::exit(47)
}

fn main() -> Result<()> {
    let args = Args::parse();

    let output = match args.command {
        Commands::ReverseComplement { seqs } => build_reverse_complement(seqs),
        Commands::Length { seq } => get_seq_length(seq),
        Commands::GCContent { seqs } => gc_content(seqs),
        Commands::PairwiseLocal {
            seqs,
            gap_open,
            gap_extend,
            hide_coords,
            try_rc,
            line_width,
            use_0_based_coords,
        } => {
            let display_opts = DisplayOptions {
                hide_coords,
                try_rc,
                line_width,
                use_0_based_coords,
            };
            pairwise(
                AlignmentCommand::Local,
                seqs,
                gap_open,
                gap_extend,
                display_opts,
            )
        }
        Commands::PairwiseSemiglobal {
            seqs,
            gap_open,
            gap_extend,
            hide_coords,
            try_rc,
            line_width,
            use_0_based_coords,
        } => {
            let display_opts = DisplayOptions {
                hide_coords,
                try_rc,
                line_width,
                use_0_based_coords,
            };

            pairwise(
                AlignmentCommand::Semiglobal,
                seqs,
                gap_open,
                gap_extend,
                display_opts,
            )
        }
        Commands::PairwiseGlobal {
            seqs,
            gap_open,
            gap_extend,
            hide_coords,
            try_rc,
            line_width,
            use_0_based_coords,
        } => {
            let display_opts = DisplayOptions {
                hide_coords,
                try_rc,
                line_width,
                use_0_based_coords,
            };
            pairwise(
                AlignmentCommand::Global,
                seqs,
                gap_open,
                gap_extend,
                display_opts,
            )
        }
    };

    match output {
        Ok(text) => {
            println!("{}", text);
            Ok(())
        }
        Err(e) => abort(&format!("Error: {:?}", e)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let seqs = vec!["GATTACA".to_string()];
        let rc = build_reverse_complement(seqs).unwrap();
        assert_eq!("TGTAATC".to_string(), rc);
    }

    #[test]
    fn test_reverse_complement_with_pairwise_cruft() {
        let seqs = vec!["GATT".to_string(), "ACA-TTGA".to_string()];
        let rc = build_reverse_complement(seqs).unwrap();
        assert_eq!("TCAA-TGT AATC".to_string(), rc);
    }

    #[test]
    fn test_length() {
        let seqs = vec!["GATTACA".to_string()];
        let length = get_seq_length(seqs).unwrap();
        assert_eq!(length, 7.to_string());
    }

    #[test]
    fn test_length_spaces() {
        let seqs = vec!["GAT".to_string(), "C".to_string(), "TACA".to_string()];
        let length = get_seq_length(seqs).unwrap();
        assert_eq!(length, 8.to_string());
    }

    #[test]
    fn test_length_gaps() {
        let seqs = vec!["GAT-CT ACA".to_string()];
        let length = get_seq_length(seqs).unwrap();
        assert_eq!(length, 8.to_string());
    }

    #[test]
    fn test_gc_content() {
        let seqs = vec!["GGG".to_string(), "GAA-TA".to_string()];
        let gc = gc_content(seqs).unwrap();
        assert_eq!(gc, "0.5000000000000000");
    }

    #[test]
    fn test_compute_gc_content_0() {
        let seqs = vec!["AT".to_string(), "TTAA".to_string()];
        let gc = compute_gc_content(seqs).unwrap();
        assert_eq!(gc, 0.0);
    }

    #[test]
    fn test_compute_gc_content_25() {
        let seqs = vec!["GG".to_string(), "TTTAAA".to_string()];
        let gc = compute_gc_content(seqs).unwrap();
        assert_eq!(gc, 0.25);
    }

    #[test]
    fn test_compute_gc_content_50() {
        let seqs = vec!["GGG".to_string(), "GAA-TA".to_string()];
        let gc = compute_gc_content(seqs).unwrap();
        assert_eq!(gc, 0.5);
    }

    #[test]
    fn test_compute_gc_content_75() {
        let seqs = vec!["GGTT".to_string(), "CCGG".to_string()];
        let gc = compute_gc_content(seqs).unwrap();
        assert_eq!(gc, 0.75);
    }

    #[test]
    fn test_compute_gc_content_100() {
        let seqs = vec!["GG".to_string(), "GGG".to_string()];
        let gc = compute_gc_content(seqs).unwrap();
        assert_eq!(gc, 1.0);
    }

    #[test]
    fn test_pairwise_local() {
        let opts = DisplayOptions {
            hide_coords: false,
            line_width: 60,
            try_rc: false,
            use_0_based_coords: true,
        };
        let actual = pairwise(
            AlignmentCommand::Local,
            vec!["ACAGT".to_string(), "ACGT".to_string()],
            2,
            1,
            opts
        )
        .unwrap();
        let expected = "3 GT 5\n  ||\n2 GT 4";
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_pairwise_semiglobal() {
        let opts = DisplayOptions {
            hide_coords: false,
            line_width: 60,
            try_rc: false,
            use_0_based_coords: true,
        };
        let actual = pairwise(
            AlignmentCommand::Semiglobal,
            vec!["ACAGT".to_string(), "ACGT".to_string()],
            2,
            1,
            opts
        )
        .unwrap();
        let expected = "0 ACAGT 5\n  || ||\n0 AC-GT 4";
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_pairwise_global() {
        let opts = DisplayOptions {
            hide_coords: false,
            line_width: 60,
            try_rc: false,
            use_0_based_coords: true,
        };
        let actual = pairwise(
            AlignmentCommand::Global,
            vec!["GGGGCCCCGGGGACAGT".to_string(), "ACGT".to_string()],
            2,
            1,
            opts
        )
        .unwrap();
        let expected = "0 GGGGCCCCGGGGACAGT 17\n              || ||\n0 ------------AC-GT 4";
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_pairwise_semiglobal_hide_coords() {
        let opts = DisplayOptions {
            hide_coords: true,
            line_width: 60,
            try_rc: false,
            use_0_based_coords: true,
        };
        let actual = pairwise(
            AlignmentCommand::Semiglobal,
            vec!["ACAGT".to_string(), "ACGT".to_string()],
            2,
            1,
            opts
        )
        .unwrap();
        let expected = "ACAGT\n|| ||\nAC-GT";
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_pairwise_semiglobal_tryrc() {
        let opts = DisplayOptions {
            hide_coords: false,
            line_width: 60,
            try_rc: true,
            use_0_based_coords: true,
        };
        let actual = pairwise(
            AlignmentCommand::Semiglobal,
            vec!["TGTAATC".to_string(), "GGCGATTACAATGACA".to_string()],
            2,
            1,
            opts
        )
        .unwrap();
        let expected = "7 GATTACA 0\n  |||||||\n3 GATTACA 10";
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_pairwise_semiglobal_high_gap_penalties() {
        let opts = DisplayOptions {
            hide_coords: false,
            line_width: 60,
            try_rc: false,
            use_0_based_coords: true,
        };
        let actual = pairwise(
            AlignmentCommand::Semiglobal,
            vec!["ACGT".to_string(), "ACAAAAGT".to_string()],
            5,
            5,
            opts
        )
        .unwrap();
        let expected = "0 ACGT 4\n  |.||\n4 AAGT 8";
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_pairwise_semiglobal_zero_gap_penalties() {
        let opts = DisplayOptions {
            hide_coords: false,
            line_width: 60,
            try_rc: true,
            use_0_based_coords: true,
        };
        let actual = pairwise(
            AlignmentCommand::Semiglobal,
            vec!["ACGT".to_string(), "ACAAAAGT".to_string()],
            0,
            0,
            opts
        )
        .unwrap();
        let expected = "0 AC----GT 4\n  ||    ||\n0 ACAAAAGT 8";
        assert_eq!(actual, expected);
    }
}
