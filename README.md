# Biotools

![Build](https://github.com/jimrybarski/biotools/actions/workflows/build.yml/badge.svg) ![Tests](https://github.com/jimrybarski/biotools/actions/workflows/tests.yml/badge.svg)

Installation: `cargo install biotools`  

Build: `git clone git@github.com:jimrybarski/biotools.git && cd biotools && cargo build --release`  

You can add `--help` after any subcommand for subcommand-specific options.  

```
Simple bioinformatics tools for sequence analysis and manipulation

Usage: biotools <COMMAND>

Subcommands:
  reverse-complement   Converts a nucleic acid sequence to its reverse complement.
  length               Computes the length of a sequence.
  gc-content           Computes the GC content of a nucleic acid sequence.
  pairwise-local       Performs a local pairwise alignment of two sequences.
  pairwise-semiglobal  Performs a semiglobal pairwise alignment of two sequences.
  pairwise-global      Performs a global pairwise alignment of two sequences.
  help                 Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```
