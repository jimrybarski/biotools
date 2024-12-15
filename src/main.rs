extern crate bio;
use anyhow::{bail, Context, Result};
use bio::alignment::pairwise::Aligner;
use bio::alignment::{Alignment, AlignmentOperation};
use bio::alphabets::dna::revcomp;
use std::env;

const GAP_OPEN_SCORE: i32 = -2;
const GAP_EXTEND_SCORE: i32 = -1;

enum AlignmentCommand {
    Local,
    Global,
    Semiglobal,
}

enum Command {
    ReverseComplement,
    StringLength,
    GCContent,
    PairwiseLocal,
    PairwiseSemiglobal,
    PairwiseGlobal,
}

fn parse_command(potential_command: &str) -> Result<Command> {
    match potential_command {
        "reverse-complement" => Ok(Command::ReverseComplement),
        "length" => Ok(Command::StringLength),
        "gc-content" => Ok(Command::GCContent),
        "pairwise-local" => Ok(Command::PairwiseLocal),
        "pairwise-semiglobal" => Ok(Command::PairwiseSemiglobal),
        "pairwise-global" => Ok(Command::PairwiseGlobal),
        _ => bail!("Invalid command"),
    }
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

fn get_string_length(user_input: Vec<String>) -> Result<String> {
    Ok(user_input[0]
        .chars()
        .filter(|ch| *ch != '-')
        .filter(|ch| *ch != ' ')
        .collect::<Vec<_>>()
        .len()
        .to_string())
}

fn gc_content(user_input: Vec<String>) -> Result<String> {
    let mut gc: f64 = 0.0;
    let mut at: f64 = 0.0;
    for chunk in user_input {
        for ch in chunk.to_uppercase().chars() {
            match ch {
                'A' | 'T' => {
                    at += 1.0;
                }
                'C' | 'G' => {
                    gc += 1.0;
                }
                _ => {}
            }
        }
    }
    let ratio = gc / (at + gc);
    Ok(format!("{:.16}", ratio))
}

fn format_alignment(alignment: Alignment, a: String, b: String) -> String {
    let mut alignment_string = String::new();
    let mut a_pretty = String::new();
    let mut a_index = alignment.xstart;
    let mut b_pretty = String::new();
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
    format!("{a_pretty}\n{alignment_string}\n{b_pretty}")
}

fn pairwise(alignment_command: AlignmentCommand, user_input: Vec<String>) -> Result<String> {
    if user_input.len() != 2 {
        bail!("Pairwise comparison needs exactly two sequences");
    }
    let a = user_input[0].clone();
    let b = user_input[1].clone();
    let a_bytes = a.as_bytes();
    let b_bytes = b.as_bytes();
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    
    let mut aligner = Aligner::with_capacity(a.len(), b.len(), GAP_OPEN_SCORE, GAP_EXTEND_SCORE, &score);

    let alignment = match alignment_command {
        AlignmentCommand::Local => aligner.local(a_bytes, b_bytes),
        AlignmentCommand::Semiglobal => aligner.semiglobal(a_bytes, b_bytes),
        AlignmentCommand::Global => aligner.global(a_bytes, b_bytes),
    };

    let pretty_alignment = format_alignment(alignment, a, b);
    Ok(pretty_alignment)
}

fn abort(error_message: &str) -> ! {
    eprintln!("biotools error: {}", error_message);
    std::process::exit(47)
}

fn main() -> Result<()> {
    let user_input: Vec<String> = env::args().skip(2).collect();
    if let Some(command) = env::args().nth(1) {
        let output = match parse_command(&command) {
            Ok(Command::ReverseComplement) => build_reverse_complement(user_input),
            Ok(Command::StringLength) => get_string_length(user_input),
            Ok(Command::GCContent) => gc_content(user_input),
            Ok(Command::PairwiseLocal) => pairwise(AlignmentCommand::Local, user_input),
            Ok(Command::PairwiseSemiglobal) => pairwise(AlignmentCommand::Semiglobal, user_input),
            Ok(Command::PairwiseGlobal) => pairwise(AlignmentCommand::Global, user_input),
            Err(e) => Err(e),
        };
        match output {
            Ok(text) => {
                println!("{}", text);
                Ok(())
            }
            Err(e) => abort(&format!("Error: {:?}", e)),
        }
    } else {
        abort("Missing command!");
    }
}
