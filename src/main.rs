#![feature(test)]

extern crate bio;
extern crate test;
use std::env;
use bio::alphabets::dna::{revcomp,complement};
use std::string::FromUtf8Error;

enum Command {
    Complement,
    ReverseComplement
}

enum Error {
    InvalidCommand
}

fn parse_command(potential_command: &str) -> Result<Command, Error>{
    match potential_command {
        "rc" => Ok(Command::ReverseComplement),
        "c" => Ok(Command::Complement),
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

    #[bench]
    #[allow(unused)]
    fn bench_complement(b: &mut Bencher) {
        b.iter(|| {
            let complement = build_complement("AAAACGTGGGGGGATCGACGACACACTGGATAATATATACGACTAGCCTACGATCGATCGT").unwrap();
        });
    }

    #[bench]
    #[allow(unused)]
    fn bench_reverse_complement(b: &mut Bencher) {
        b.iter(|| {
            let complement = build_reverse_complement("AAAACGTGGGGGGATCGACGACACACTGGATAATATATACGACTAGCCTACGATCGATCGT").unwrap();
        });
    }
}