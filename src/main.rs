extern crate bio;
use std::env;
use bio::alphabets::dna::revcomp;
use std::string::FromUtf8Error;


fn reverse_complement(string: String) -> Result<String, FromUtf8Error> {
    String::from_utf8(revcomp(string
                              .trim()
                              .replace(" ", "")
                              .as_bytes()))
}


fn main() {
  let rc: String = env::args()
                   .skip(1)
                   .collect();
        
  match reverse_complement(rc) {
    Ok(o) => { println!("{}", o); },
    Err(_) => { std::process::exit(1); }
  }
}
