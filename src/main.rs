extern crate bio;
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

fn build_reverse_complement(user_input: String) -> Result<String, FromUtf8Error> {
    String::from_utf8(revcomp(user_input
                              .trim()
                              .replace(" ", "")
                              .as_bytes()))
}

fn build_complement(user_input: String) -> Result<String, FromUtf8Error> {
    String::from_utf8(user_input.trim()
                                .replace(" ", "")
                                .bytes()
                                .map(complement)
                                .collect())
}

fn abort(error_message: &str) {
    println!("{}", error_message);
    std::process::exit(1)
}

fn main() {
    let user_input: String = env::args()
                             .skip(2)
                             .collect();
    if let Some(command) = env::args().nth(1) {
        let output = match parse_command(&command) {
            Ok(Command::Complement) => {
                build_complement(user_input)
            },
            Ok(Command::ReverseComplement) => {
                build_reverse_complement(user_input)
            }
            Err(_) => { return abort("seqtools-cli: Invalid command!") }
        };
        match output {
            Ok(text) => { println!("{}", text); }
            Err(_) => { return abort("seqtools-cli: Invalid sequence!") }
        }
    }
}
