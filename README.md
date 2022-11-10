# split2sample

This program is forked and based on the splitp Rust SPLiT-seq read pre-processing script.
Initial implementation was performed by Rob P. Idea, logics and later implementations came from Stefan L.

### Usage

The `split2sample` program takes several arguments.  The usage can be printed 
from the command line using `split2sample -h`.

```
split2samples 0.1.0
Stefan L. <stefan.lang@med.lu.se>, Rob P. <rob@cs.umd.edu>
Split a pair of BD rhapsody fastq files (R1 and R2) into sample specific fastq pairs

USAGE:
    split2samples --reads <READS> --file <FILE> --specie <SPECIE> --outpath <OUTPATH>

OPTIONS:
    -f, --file <FILE>          the input R2 samples file
    -h, --help                 Print help information
    -o, --outpath <OUTPATH>    the outpath
    -r, --reads <READS>        the input R1 reads file
    -s, --specie <SPECIE>      the specie of the library [mouse, human]
    -V, --version              Print version information

```

This debug branch is meant to implement simple a multi processor approach:
Instead of implementing this in the script the script will check for the exisance of files before finiching a analysis.

The analysis is split up in two parts:
1. the fastq files get scanned for sample tags. For each sample tag the cell ids connected to this sample tag are stored.
2. the fastq file gets processed again identifying all cells and writing them out into sample specific fastq files.

The way I want to make this muliti processor is to split the program into two run modes:
1. identifyCells will identify the cell<->sample conectios and write them to disk
2. splitFastq will use the before collected cell<->sample data and split the fastq.

Applying this to smaller fastq files should speed up the process.

I am currently splitting the fastq files using fastqsplitter. But this is not part of the rust script.

```
split2samples -r Data_R1.fastq.gz -f Data_R2.fastq.gz -o output -s mouse
```

### Installation

```
git clone https://github.com/stela2502/split2samples
cd split2samples
cargo build
cp target/debug/split2samples /usr/bin
``` 

To run the test data (a tiny bit of a real dataset):

```
./target/debug/split2samples -r testData/testData_R1.fastq.gz -f testData/testData_R2.fastq.gz -o testData/output -s mouse
```

Or the bigger test data with 1e+5 reads:
```
./target/debug/split2samples -r testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse
```
### Limitations / differences

This program is totally untested and under heavy development.
This is only the first draft - let's see where this heads to.
