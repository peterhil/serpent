# Serpent


## Explore DNA data with Serpent

Serpent is an exploration into DNA and RNA sequences, nucleotide
bases, codons, amino acids and genome data.


## Tools provided

### Work with FASTA files and sequences

* `serpent cat`: concatenate and print FASTA files
* `serpent find`: find FASTA files in directories
* `serpent find -s`: find and print FASTA sequences in files and directories

### Convert data

* `serpent encode`: Convert data into different encoded representations
* `serpent decode`: Map codons into numbers 0...64

### Analyse and plot FASTA data visually

* `serpent ac`: print and plot autocorrelation on DNA and RNA sequences
* `serpent fft`: plot FFTs on DNA and RNA sequences
* `serpent hist`: plot histogram statistics
* `serpent image`: visualise DNA and RNA data as images
* `serpent seq`: plot sequence count statistics

### Statistics

* `serpent codons`: Print codon statistics
* `serpent pep`: Print peptide statistics

See `serpent -h` for all subcommands and `serpent <subcommand> -h` for options!


## Motivation

I have wanted to explore DNA data in order to to learn and maybe
invent some compression algorithms for DNA data for about two decades.
