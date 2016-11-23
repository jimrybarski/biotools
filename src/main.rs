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
  let rc = env::args()
           .nth(1)
           .ok_or("Please specify input".to_owned())
           .and_then(|user_input| reverse_complement(user_input)
                                  .map_err(|err| err.to_string()));

  match rc {
    Ok(result) => { println!("{}", result); },
    Err(err) => { println!("{}", err); }
  }
}
