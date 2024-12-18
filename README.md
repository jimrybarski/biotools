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

### Reverse complement

There are no options. Spaces and dashes are allowed to permit usage with gap-containing pairwise alignments.

```bash
$ biotools reverse-complement GATTACA
TGTAATC
$ biotools reverse-complement GATT ACA-TTGA
TCAA-TGT AATC
```

### Length

There are no options. Spaces and dashes are allowed to permit usage with gap-containing pairwise alignments.
```
$ biotools length GATTACA
7
$ biotools length GAT  C   TACA
8
$ biotools length GAT-CTACA
8
```

### GC content

There are no options. Spaces and dashes are allowed to permit usage with gap-containing pairwise alignments.
```
$ biotools gc-content GGG GAA-TA
0.5000000000000000
```

### Pairwise alignment

There are three pairwise alignment commands, for local, semiglobal and global alignments.

```
$ biotools pairwise-local ACAGT ACGT
3 GT 5
  ||
2 GT 4

$ biotools pairwise-semiglobal ACAGT ACGT
0 ACAGT 5
  || ||
0 AC-GT 4

$ biotools pairwise-global GGGGCCCCGGGGACAGT ACGT
0 GGGGCCCCGGGGACAGT 17
              || ||
0 ------------AC-GT  4
```

Coordinates can be disabled:

```
biotools pairwise-semiglobal ACAGT ACGT --hide-coords
ACAGT
|| ||
AC-GT
```

You can have biotools try both the forward sequence and reverse complement and then pick the one with the best alignment score. The text `RC` will be displayed the right of the first sequence if the reverse complement was better.

```
$ biotools pairwise-semiglobal TGTAATC GGCGATTACAATGACA
0 TGTAATC  7
  |..|||.
6 TACAATG 13

$ biotools pairwise-semiglobal TGTAATC GGCGATTACAATGACA --try-rc
0 GATTACA  7 RC
  |||||||
3 GATTACA 10
```

You can adjust the gap penalties (substitutions are always penalized a value of 1). These are given as positive numbers.

```
$ biotools pairwise-semiglobal ACGT ACAAAAGT --gap-open 5 --gap-extend 5
0 ACGT 4
  |.||
4 AAGT 8

$ biotools pairwise-semiglobal ACGT ACAAAAGT --gap-open 0 --gap-extend 0
0 AC----GT 4
  ||    ||
0 ACAAAAGT 8
```
