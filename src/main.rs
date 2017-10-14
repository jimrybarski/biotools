#![feature(test)]

extern crate bio;
extern crate test;
use std::env;
use bio::alphabets::dna::{revcomp,complement};
use std::string::FromUtf8Error;

enum Command {
    Complement,
    ReverseComplement,
    StringLength,
    GCContent
}

enum Error {
    InvalidCommand
}

fn parse_command(potential_command: &str) -> Result<Command, Error>{
    match potential_command {
        "rc" => Ok(Command::ReverseComplement),
        "c" => Ok(Command::Complement),
        "len" => Ok(Command::StringLength),
        "gc" => Ok(Command::GCContent),
        _ => Err(Error::InvalidCommand)
    }
}

fn build_reverse_complement(user_input: &str) -> Result<String, FromUtf8Error> {
    String::from_utf8(revcomp(user_input
                              .trim()
                              .replace(" ", "")
                              .as_bytes()))
}

fn build_complement(user_input: &str) -> Result<String, FromUtf8Error> {
    String::from_utf8(user_input.trim()
                                .replace(" ", "")
                                .bytes()
                                .map(complement)
                                .collect())
}

fn get_string_length(user_input: &str) -> Result<String, FromUtf8Error> {
    Ok(user_input.len().to_string())
}

fn gc_content(user_input: &str) -> Result<String, FromUtf8Error> {
    let mut gc: f64 = 0.0;
    let mut at: f64 = 0.0;
    for ch in user_input
              .to_uppercase()
              .chars() {
        match ch {
            'A' | 'T' => { at += 1.0; },
            'C' | 'G' => { gc += 1.0; },
            _ => {}
        }
    }
    let ratio = gc / (at + gc);
    Ok(format!("{:.16}", ratio))
}

fn abort(error_message: &str) {
    println!("seqtools-cli error: {}", error_message);
    std::process::exit(1)
}

fn main() {
    let user_input: String = env::args()
                             .skip(2)
                             .collect();
    if let Some(command) = env::args().nth(1) {
        let output = match parse_command(&command) {
            Ok(Command::Complement) => { build_complement(&user_input) },
            Ok(Command::ReverseComplement) => { build_reverse_complement(&user_input) },
            Ok(Command::StringLength) => { get_string_length(&user_input) },
            Ok(Command::GCContent) => { gc_content(&user_input) },
            Err(_) => { return abort("Invalid command!") }
        };
        match output {
            Ok(text) => { println!("{}", text); },
            Err(_) => { return abort("Invalid sequence!") }
        }
    } else {
        abort("Missing command!");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[test]
    fn test_build_complement() {
        let complement = build_complement("AAAACGT").unwrap();
        assert_eq!(complement, String::from("TTTTGCA"));
    }

    #[test]
    fn test_build_reverse_complement() {
        let complement = build_reverse_complement("AAAACGT").unwrap();
        assert_eq!(complement, "ACGTTTT");
    }

    #[test]
    fn test_gc_content() {
        let ratio = gc_content("ATATGCGC").unwrap();
        assert_eq!(ratio, "0.5000000000000000");
    }

    #[test]
    fn test_gc_content2() {
        let ratio = gc_content("ATATTTTTA").unwrap();
        assert_eq!(ratio, "0.0000000000000000");
    }

    #[test]
    fn test_gc_content3() {
        let ratio = gc_content("GGGCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGCGCGCGCGGCGC").unwrap();
        assert_eq!(ratio, "1.0000000000000000");
    }

    #[test]
    fn test_gc_content4() {
        let ratio = gc_content("GATTACA").unwrap();
        assert_eq!(ratio, "0.2857142857142857");
    }
    // The benchmarks are designed to take around 1 microsecond to run in V1.20
    //
    #[bench]
    #[allow(unused)]
    fn bench_complement(b: &mut Bencher) {
        b.iter(|| {
            let complement = build_complement("AAAACGTGGGGGGATCGACGACACA").unwrap();
        });
    }

    #[bench]
    #[allow(unused)]
    fn bench_reverse_complement(b: &mut Bencher) {
        b.iter(|| {
            let complement = build_reverse_complement("AAAACGTGGGGGGATCGACGACACA").unwrap();
        });
    }
}