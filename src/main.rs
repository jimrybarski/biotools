extern crate bio;
use std::env;
use bio::alphabets;


fn main() {
    let user_input = match env::args().nth(1) {
        Some(x) => x,
        None => { println!("You need to enter a sequence."); std::process::exit(1); }
    };
    let rc = alphabets::dna::revcomp(user_input.trim().replace(" ", "").as_bytes());
    let output = match String::from_utf8(rc) {
        Ok(x) => x,
        Err(_) => { println!("Invalid sequence."); std::process::exit(1); }
    };

    println!("{}", output);
    std::process::exit(0);
}
