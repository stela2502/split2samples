# This is seriousely broken!

try to really understand what is going on there:

```
 target/release/quantifyRhapsody -r  .\testData\OneSingleCell.369083.R2.fastq.gz -f .\testData\OneSingleCell.369083.R1.fastq.gz -o testData/output_one -s mouse  -e testData/genes.fasta
 -a testData/MyAbSeqPanel.fasta -m 0
```

This should have quite some Igha reads. And I do not get why they do not show. I have manually checked these sequences:

* GCGTGTATCA
* ACTCTGGAA
* ACTCTGGAAACAGGGT

Why does that not show in the results?! This is total bullshit what I have done here?!?

# split2sample

This program is forked and based on the splitp Rust SPLiT-seq read pre-processing script.
Initial implementation was performed by Rob P. Idea, logics and later implementations came from Stefan L.

# Usage

The `split2sample` program takes several arguments.  The usage can be printed 
from the command line using `split2sample -h`.

```
./target/debug/split2samples -h
split2samples 0.1.0
Stefan L. <stefan.lang@med.lu.se>, Rob P. <rob@cs.umd.edu>
Split a pair of BD rhapsody fastq files (R1 and R2) into sample specific fastq pairs

USAGE:
    split2samples --reads <READS> --file <FILE> --specie <SPECIE> --outpath <OUTPATH> --mode <MODE>

OPTIONS:
    -f, --file <FILE>          the input R2 samples file
    -h, --help                 Print help information
    -m, --mode <MODE>          the mode of the program [fastqSplit, cellIdent, sampleSplit]
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

## Splitting a set of two fastq files

into 10 (20) fastq files. Each of the 10 fastq files gets 100 fastq entries. Afterwads the script switches to the next ansd so on.
I have tried to copy the function of fastqsplitter. I have not tested this part.

```
split2samples -m fastqSplit  -r Data_R1.fastq.gz -f Data_R2.fastq.gz -o output -s mouse
```

This will (hopefully) create a total of 20 fastq files in the outpath

## Identifying sample -> cell links

```
split2samples -m cellIdent -r Data_R1.fastq.gz -f Data_R2.fastq.gz -o output -s mouse
```

This will create a set of 12 Data_R1.fastq.gz.sample[1-12].ints.txt files that are necessary for the next step.

## Split the Fastq files into sample specififc fastq files

```
split2samples -m sampleSplit -r Data_R1.fastq.gz -f Data_R2.fastq.gz -o output -s mouse
```

This will create a set of 24 (12 +12) fastq files - two for each sample.


# Installation

```
git clone https://github.com/stela2502/split2samples
cd split2samples
cargo build --release
cp target/release/split2samples /usr/bin
cp target/release/demux10x /usr/bin
cp target/release/quantifyRhapsody /usr/bin
``` 

Do not forget the --release while building the tool. 
The test case for quantifyRhapsody would finish in 7 sec instead of ~0.5 sec (x14!)
using a AMD Ryzen 7 5700X processor and a SSD as mass storage.

To run the test data (a tiny bit of a real dataset):

```
./target/debug/split2samples -r testData/testData_R1.fastq.gz -f testData/testData_R2.fastq.gz -o testData/output -s mouse
```

Or the bigger test data with 1e+5 reads:
```
#./target/release/split2samples -m fastqSplit -r testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse
./target/release/split2samples -m cellIdent -r testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse
./target/release/split2samples -m sampleSplit -r testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse
# or the two together
./target/release/split2samples -m analysis -r testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse
```

Test case for the quantifyRhapsody:

```
target/release/quantifyRhapsody -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse  -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 30
```


# Limitations / differences

This program is totally untested and under heavy development.
This is only the first draft - let's see where this heads to.


./target/release/split2samples -m cellIdent -r testData/output_1e5/fastqSplit.1.1e5_mRNA_S1_R1_001.fastq.gz -f testData/output_1e5/fastqSplit.1.1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse


