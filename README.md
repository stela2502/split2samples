# Rustody - a tool to quickly analyze targeted BD-Rhapsody sequencings

This tool replaces the official 7-Bridges BD analysis programs first steps.

The final output from this tool is a sparse matrix of gene expression values, antibody tags and in addition a dense matrix with read counts for the sample reads.

I found the sample table most helpful in the detection of populations containing duplicate cells.

The output from here can easily be read into any single cell analysis package for downstream analysis like Seurat or Scanpy.

# Installation

You need the Rust compiler: https://www.rust-lang.org/tools/install


Then you can clone this repo and complie the code (example for a Linux system).
But it also compiles on Windows. I just never use that for actual work.


```
git clone https://github.com/stela2502/split2samples
cd split2samples
cargo build --release
cp target/release/split2samples /usr/bin
cp target/release/demux10x /usr/bin
cp target/release/quantifyRhapsody /usr/bin
cp target/release/quantifyRhapsody_noUMI /usr/bin
cp target/release/bd_cell_id_2_seq /usr/bin
cp target/release/bd_get_single_cell /usr/bin
``` 

Do not forget the --release while building the tool. 
The test case for quantifyRhapsody would finish in 7 sec instead of ~0.5 sec (x14!)
using a AMD Ryzen 7 5700X processor and a SSD as mass storage.


## Testing

To run the test data (a tiny bit of a real dataset):

```
target/release/quantifyRhapsody -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse  -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 30 -v v1
```

This will produce an output consisting of two cells. and it should run super fast.


# Usage

The `quantifyRhapsody` program takes several arguments.  The usage can be printed 
from the command line using `quantifyRhapsody -h`.

```
./target/release/quantifyRhapsody -h

USAGE:
    quantifyRhapsody.exe --reads <READS> --file <FILE> --specie <SPECIE> --outpath <OUTPATH> --expression <EXPRESSION> --antybody <ANTYBODY> --min-umi <MIN_UMI>

OPTIONS:
    -a, --antibody <ANTIBODY>        the fasta database containing the antibody tags
    -e, --expression <EXPRESSION>    the fasta database containing the genes
    -f, --file <FILE>                the input R2 fastq file
    -h, --help                       Print help information
    -m, --min-umi <MIN_UMI>          the minimum reads (sample + genes + antibody combined)
    -o, --outpath <OUTPATH>          the outpath
    -r, --reads <READS>              the input R1 fastq file
    -s, --specie <SPECIE>            the specie of the library [mouse, human]
    -V, --version                    Print version information
```

You see - this is the one I compiled on Windows 11.


# Additional Programs

There are several other programs in this package:

 1. split2samples will split the BD Rhapsody fastq files into sample spceific fastq files. This script is older and ~4 times slower in creating just the fastq files when compared to quantifyRhapsody quantifying the data.
 2. demux10x is a small spin off that actually processes 10x single cell data and searches for a set fasta entries.
 3. bd_cell_id_2_seq BD Rhapsody cells do get an ID in the results. If you want to get the sequences coding for this cells you can use this program
 4. bd_get_single_cell will select only one single cell from the fastq files.


# Limitations / differences

This program is totally untested and under heavy development.
This is only the first draft - let's see where this heads to.


