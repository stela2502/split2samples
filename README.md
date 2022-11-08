# split2sample

This program is forked and based on the splitp Rust SPLiT-seq read pre-processing script.
Main implementation was performed by Rob P. Idea and logics came from Stefan L.

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

**Please take note** that `split2samples` will create 13 files in the output folder - 12 sample specific 
and one ambig.fastq.gz containing all reads where the sample id could not be identified. 
Please note that not all of these fastq files contain data.

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
### Limitations / differences

This program is totally untested and under heavy development.
This is only the first draft - let's see where this heads to.
