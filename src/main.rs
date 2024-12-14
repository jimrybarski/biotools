extern crate bio;
use bio::alphabets::dna::revcomp;
use std::env;
use std::string::FromUtf8Error;

enum Command {
    ReverseComplement,
    StringLength,
    GCContent,
}

enum Error {
    InvalidCommand,
}

fn parse_command(potential_command: &str) -> Result<Command, Error> {
    match potential_command {
        "reverse-complement" => Ok(Command::ReverseComplement),
        "length" => Ok(Command::StringLength),
        "gc-content" => Ok(Command::GCContent),
        _ => Err(Error::InvalidCommand),
    }
}

fn build_reverse_complement(user_input: &str) -> Result<String, FromUtf8Error> {
    String::from_utf8(revcomp(user_input.trim().as_bytes()))
}

fn get_string_length(user_input: &str) -> Result<String, FromUtf8Error> {
    Ok(user_input
        .chars()
        .filter(|ch| *ch != '-')
        .filter(|ch| *ch != ' ')
        .collect::<Vec<_>>().len().to_string())
}

fn gc_content(user_input: &str) -> Result<String, FromUtf8Error> {
    let mut gc: f64 = 0.0;
    let mut at: f64 = 0.0;
    for ch in user_input.to_uppercase().chars() {
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
    let ratio = gc / (at + gc);
    Ok(format!("{:.16}", ratio))
}

fn abort(error_message: &str) {
    println!("biotools error: {}", error_message);
    std::process::exit(1)
}

fn main() {
    let user_input: String = env::args().skip(2).collect();
    if let Some(command) = env::args().nth(1) {
        let output = match parse_command(&command) {
            Ok(Command::ReverseComplement) => build_reverse_complement(&user_input),
            Ok(Command::StringLength) => get_string_length(&user_input),
            Ok(Command::GCContent) => gc_content(&user_input),
            Err(_) => return abort("Invalid command!"),
        };
        match output {
            Ok(text) => {
                println!("{}", text);
            }
            Err(_) => abort("Invalid sequence!"),
        }
    } else {
        abort("Missing command!");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn test_gc_content_with_gaps() {
        let ratio = gc_content("GATT-ACA").unwrap();
        assert_eq!(ratio, "0.2857142857142857");
    }
}

