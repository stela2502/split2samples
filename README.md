[![Rust](https://github.com/stela2502/Rustody/actions/workflows/rust.yml/badge.svg?branch=main2)](https://github.com/stela2502/Rustody/actions/workflows/rust.yml)

# Rustody - a tool to quickly analyze targeted BD-Rhapsody sequencings

This tool replaces the official 7-Bridges BD analysis programs first steps.

The final output from this tool is a sparse matrix of gene expression values, antibody tags and in addition a dense matrix with read counts for the sample reads.

I found the sample table most helpful in the detection of populations containing duplicate cells.

The output from here can easily be read into any single cell analysis package for downstream analysis like Seurat or Scanpy.

You can inspect the state of the program using this [deatiled comparison between the seven bridges BD Rhapsody pipeline and the Rhapsody output here]( ./testData/BD_results/CombinedAnalysis_scanpy_v1.2.1.ipynb).

## News

Mapper has significantly improved: both the false positive as well as false negative rate had improved.
The improvement was possible by using a needleman-wunsch inspired algorithm.
The 32 bp matches are now tolerant to not only bp mismatches, but also insertions and deletions.

Compared to the mere bp replacement matching we e.g. find almost 10x more reads from the Ighm locus in the test data.
And none of the reads I have seen so far looked like a not Ighm transcript. All I checked were also mapped to the Ighm transcripts using NCBI BLAST (I only checked the strange looking ones).

I have added the Ighm reads that were detected with both settings to this repository.


quantify_rhapsody has finally gotten a muti processor upgrade: quantify_rhapsody_multi.
I have not tested it out completetly now, but am confident it works correctly. (Final last words - I know). Just the PCR duplicates are not collected correctly as they can now only be measured in each chunk of the data. But the UMIs is all this tool does measure.
So even if the PCR duplicates are not counted correctly they will nevertheless be excluded from the final data.

The news [can be read here](./News.md).

# Installation

You need the Rust compiler: https://www.rust-lang.org/tools/install


Then you can clone this repo and complie the code (example for a Linux system).
But it also compiles on Windows. I just never use that for actual work.


```
git clone https://github.com/stela2502/Rustody
cd Rustody
cargo build --release
cp target/release/split2samples /usr/bin
cp target/release/quantify_rhapsody /usr/bin
cp target/release/quantify_rhapsody_multi /usr/bin
cp target/release/bd_cell_id_2_seq /usr/bin
cp target/release/bd_get_single_cell /usr/bin
cp target/release/get_n_cell_reads /usr/bin
cp target/release/int_2_seq /usr/bin
cp target/release/create_index /usr/bin
cp target/release/create_index_te /usr/bin
``` 

Do not forget the --release while building the tool. 
The test case for quantify_rhapsody would finish in 55 sec instead of 3 sec (> x15!)
using a AMD Ryzen 7 5700X processor and a SSD as mass storage.


## Testing

The test data is a tiny bit of an (unpublished) real data set.
It is included here solely to present this mapping procedure.

To run the test data (a tiny bit of a real dataset):

```
target/release/quantify_rhapsody -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v1
```

Latest results:

```
Analysis will stop after having processed 18446744073709551615 fastq entries containing a cell info

init models
the log file: Mapping_log.txt
Changing the expression start gene id to 462
After indexing all fastq files we have the following indices:
the mRNA index:
I have 16929 kmers for 462 genes with 0% duplicate entries
gene names like 'Rpl36a'
gene_ids range from Some(0) to Some(461)

the sample id index:
I have 370 kmers for 12 genes with 0.010398613% duplicate entries
gene names like 'Sample12'
gene_ids range from Some(467) to Some(478)

and the antibodies index:
I have 97 kmers for 5 genes with 0.3445946% duplicate entries
gene names like 'IgM'
gene_ids range from Some(462) to Some(466)

Writing index version 5
with kmer_len 32
And a total of 16929 data entries
   0.10 mio reads (69.86% with cell_id, 66.12% with gene_id)                                                                                                                                                       

Writing outfiles ...
filtering cells
Dropping cell with too little counts (n=36109)
75 cells have passed the cutoff of 10 umi counts per cell
writing gene expression
sparse Matrix: 75 cell(s), 104 gene(s) and 520 entries written to path Ok("testData/output_1e5/BD_Rhapsody_expression"); 
Writing Antibody counts
sparse Matrix: 75 cell(s), 5 gene(s) and 85 entries written to path Ok("testData/output_1e5/BD_Rhapsody_antibodies"); 
Writing samples table
dense matrix: 75 cell written

Summary:
total      reads  : 100000 reads
no cell ID reads  : 16060 reads (16.06% of total)
no gene ID reads  : 3736 reads (3.74% of total)
N's or too short  : 14080 reads (14.08% of total)
cellular reads    : 69860 reads (69.86% of total)

collected read counts:
expression reads  : 47845 reads (68.49% of cellular)
antibody reads    : 17304 reads (24.77% of cellular)
sample reads      : 975 reads (1.40% of cellular)

reported umi counts:
expression reads  : 654 reads (0.94% of cellular)
antibody reads    : 235 reads (0.34% of cellular)
sample reads      : 12 reads (0.02% of cellular)

PCR duplicates or bad cells: 68959 reads (98.71% of cellular)

timings:
   overall run time 0 h 0 min 3 sec 158 millisec
   file-io run time 0 h 0 min 3 sec 125 millisec
single-cpu run time 0 h 0 min 0 sec 0 millisec
 multi-cpu run time 0 h 0 min 0 sec 25 millisec


Cell->Sample table written to "testData/output_1e5/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 3 sec 158milli sec

```


Or test the multi processor version:

```
target/release/quantify_rhapsody_multi -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v1
```

```
...

timings:
   overall run time 0 h 0 min 0 sec 975 millisec
   file-io run time 0 h 0 min 0 sec 195 millisec
single-cpu run time 0 h 0 min 0 sec 71 millisec
 multi-cpu run time 0 h 0 min 0 sec 479 millisec

```

[An example analysis of this data is available here:]( ./testData/BD_results/CombinedAnalysis_scanpy_v1.2.1.ipynb).


# Usage

## quantfyRhapsody

The `quantifyRhapsody` program takes several arguments.  The usage can be printed 
from the command line using `quantifyRhapsody -h`.

```
target/release/quantify_rhapsody  -h
rustody 1.2.1
Stefan L. <stefan.lang@med.lu.se>
Quantifies a DB Rhapsody experiment and creates sparse matrix outfiles. You need quite long R1 and
R2 reads for this! (>70R1 and >70R2 \[v1\] and 52 bp reads for v2.96 and v2.384)

USAGE:
    quantify_rhapsody [OPTIONS] --reads <READS> --file <FILE> --specie <SPECIE> --outpath <OUTPATH> --min-umi <MIN_UMI> --version <VERSION>

OPTIONS:
    -a, --antibody <ANTIBODY>          the fasta database containing the antibody tags
    -e, --expression <EXPRESSION>      the fasta database containing the genes
    -f, --file <FILE>                  the input R2 samples file
        --gene-kmers <GENE_KMERS>      minimal sequencing quality [default: 32]
    -h, --help                         Print help information
    -i, --index <INDEX>                a pre-defined index folder produced by the cerateIndex scipt
    -m, --min-umi <MIN_UMI>            the minimum (UMI) reads per cell (sample + genes + antibody
                                       combined)
        --max-reads <MAX_READS>        Optional: end the analysis after processing <max_reads> cell
                                       fastq entries [default: 18446744073709551615]
        --min-quality <MIN_QUALITY>    minimal sequencing quality [default: 25]
    -o, --outpath <OUTPATH>            the outpath
    -r, --reads <READS>                the input R1 reads file
    -s, --specie <SPECIE>              the specie of the library [mouse, human]
    -v, --version <VERSION>            the version of beads you used v1, v2.96 or v2.384

```

The multiprocessor version has only one more option:
```
    -n, --num-threads <NUM_THREADS>    how many threads to use to analyze this (default max
                                       available)
```

This project compiles and runs both on Linux and Windows11 (MaxOS - untested).

Both R1 and R2 files can also be comma separated lists of fastq files.
As always the same order must be preserved in both lists.

The benefit between running all files separately is that the umis are controlled for all fastq files and therefore a 'better' result is obtained than processing the fastq files separately and merging the resulting data.

```
target/release/quantify_rhapsody -r  testData/1e5_mRNA_S1_R1_001.fastq.gz,testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz,testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v1
```

output of that call (version 1.2.1).

```
dense matrix: 75 cell written

Summary:
total      reads  : 200000 reads
no cell ID reads  : 32120 reads (16.06% of total)
no gene ID reads  : 7472 reads (3.74% of total)
N's or too short  : 28160 reads (14.08% of total)
cellular reads    : 139720 reads (69.86% of total)

collected read counts:
expression reads  : 95690 reads (68.49% of cellular)
antibody reads    : 34608 reads (24.77% of cellular)
sample reads      : 1950 reads (1.40% of cellular)

reported umi counts:
expression reads  : 654 reads (0.47% of cellular)
antibody reads    : 235 reads (0.17% of cellular)
sample reads      : 12 reads (0.01% of cellular)

PCR duplicates or bad cells: 138819 reads (99.36% of cellular)

timings:
   overall run time 0 h 0 min 6 sec 181 millisec
   file-io run time 0 h 0 min 6 sec 145 millisec
single-cpu run time 0 h 0 min 0 sec 0 millisec
 multi-cpu run time 0 h 0 min 0 sec 28 millisec


Cell->Sample table written to "testData/output_1e5/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 6 sec 182milli sec

```


## createIndex

BD Rhapsody now also supports whole genome transcriptomics. For this setting it becomes vital to have a mapping index that is capable of storing the genome information.
In order to do this a 2bit based 8mer numeric initial mapper has been combined with 32bp full length match to build up a 40bp full length match.

The creatIndex produces these indices for quantifyRhapsody. It only needs a geneome fasta file and a genome gtf or gff annotation file:

```
Rustody 1.0.0
Stefan L. <stefan.lang@med.lu.se>
Create a binary indes for the quantify_rhapsody program. Make sure to not cat any gzipped files for
Rust. At least here we will only read from the first stream in a catted gz file!

USAGE:
    create_index [OPTIONS]

OPTIONS:
    -a, --antibody <ANTIBODY>          [default: testData/MyAbSeqPanel.fasta]
    -f, --file <FILE>                  the fasta genome data [default: testData/testGenes.fa.gz]
    -g, --gtf <GTF>                    the gff/gtf gene information table (ONE gzipped stream - do
                                       NOT cat these files!) [default: testData/testGenes.gtf.gz]
        --gene-kmers <GENE_KMERS>      the mapping kmer length [default: 32]
        --geneid <GENEID>              the string to check for gene id (default gene_id) [default:
                                       gene_id]
        --genename <GENENAME>          the string to check for gene level names (default gene_name)
                                       [default: gene_name]
    -h, --help                         Print help information
    -n, --num-threads <NUM_THREADS>    how many threads to use to analyze this (default 1)
    -o, --outpath <OUTPATH>            the outpath [default: testData/mapperTest]
        --text <text>                  create text outfile instead of binary [default: false]
        --transcript <TRANSCRIPT>      the string to check for transcript levels names
                                       (transcript_id) [default: transcript_id]
    -V, --version                      Print version information
```

```
./target/release/create_index -n 12 -f testData/mapperTest/Juan_genes.fa.gz -g testData/mapperTest/Juan_genes.fixed.gtf.gz -o testData/mapperTest/index > testData/mapperTest/mRNA.fa
```

**Performance**

```
We created this fast_mapper object:
I have 45846 kmers for 2481 genes with 0% duplicate entries
gene names like 'ENSMUST00000156250.1'
gene_ids range from Some(0) to Some(2480)

Writing index version 5
with kmer_len 32
And a total of 45846 data entries
A text version of the index was written to testData/mapperTest/index - you can simply remove that after an optional inspection.
   overall run time 0 h 0 min 2 sec 249 millisec
   file-io run time 0 h 0 min 0 sec 196 millisec
single-cpu run time 0 h 0 min 0 sec 359 millisec
 multi-cpu run time 0 h 0 min 1 sec 570 millisec

finished in 0 h 0 min 2 sec 250 milli sec
```

**File size**
```
 39K index.1.gene.txt
6,1M index.1.Index

(22M index.1.Index.txt
5,8M mRNA.fa)
```

**And does that work?**


```
target/release/quantify_rhapsody_multi -r testData/cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f testData/cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o testData/BD_results/Rustody_S1 -s mouse  -e testData/2276_20220531_chang_to_rpl36a_amplicons.fasta -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96


...

Writing outfiles ...
filtering cells
Dropping cell with too little counts (n=54)
67 cells have passed the cutoff of 200 umi counts per cell
writing gene expression
sparse Matrix: 67 cell(s), 305 gene(s) and 4199 entries written to path Ok("testData/BD_results/Rustody_S1/BD_Rhapsody_expression"); 
Writing Antibody counts
sparse Matrix: 67 cell(s), 3 gene(s) and 139 entries written to path Ok("testData/BD_results/Rustody_S1/BD_Rhapsody_antibodies"); 
Writing samples table
dense matrix: 67 cell written

Summary:
total      reads  : 500000 reads
no cell ID reads  : 17436 reads (3.49% of total)
no gene ID reads  : 5893 reads (1.18% of total)
N's or too short  : 152 reads (0.03% of total)
cellular reads    : 482412 reads (96.48% of total)

collected read counts:
expression reads  : 472057 reads (97.85% of cellular)
antibody reads    : 3498 reads (0.73% of cellular)
sample reads      : 964 reads (0.20% of cellular)

reported umi counts:
expression reads  : 195010 reads (40.42% of cellular)
antibody reads    : 2866 reads (0.59% of cellular)
sample reads      : 834 reads (0.17% of cellular)

PCR duplicates or bad cells: 283702 reads (58.81% of cellular)

timings:
   overall run time 0 h 0 min 3 sec 452 millisec
   file-io run time 0 h 0 min 0 sec 722 millisec
single-cpu run time 0 h 0 min 0 sec 31 millisec
 multi-cpu run time 0 h 0 min 2 sec 539 millisec


Cell->Sample table written to "testData/BD_results/Rustody_S1/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 3 sec 453milli sec

```

The same as the result of that?

```
target/release/quantify_rhapsody_multi -r testData/cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f testData/cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o testData/BD_results/Rustody_S1_index -s mouse  -i testData/mapperTest/index/ -e testData/addOn.fa -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96

...

Writing outfiles ...
filtering cells
Dropping cell with too little counts (n=55)
62 cells have passed the cutoff of 200 umi counts per cell
writing gene expression
sparse Matrix: 62 cell(s), 257 gene(s) and 2647 entries written to path Ok("testData/BD_results/Rustody_S1_index/BD_Rhapsody_expression"); 
Writing Antibody counts
sparse Matrix: 62 cell(s), 3 gene(s) and 136 entries written to path Ok("testData/BD_results/Rustody_S1_index/BD_Rhapsody_antibodies"); 
Writing samples table
dense matrix: 62 cell written

Summary:
total      reads  : 500000 reads
no cell ID reads  : 17436 reads (3.49% of total)
no gene ID reads  : 180719 reads (36.14% of total)
N's or too short  : 152 reads (0.03% of total)
cellular reads    : 482412 reads (96.48% of total)

collected read counts:
expression reads  : 297231 reads (61.61% of cellular)
antibody reads    : 3498 reads (0.73% of cellular)
sample reads      : 964 reads (0.20% of cellular)

reported umi counts:
expression reads  : 129870 reads (26.92% of cellular)
antibody reads    : 2863 reads (0.59% of cellular)
sample reads      : 828 reads (0.17% of cellular)

PCR duplicates or bad cells: 348851 reads (72.31% of cellular)

timings:
   overall run time 0 h 0 min 4 sec 941 millisec
   file-io run time 0 h 0 min 0 sec 741 millisec
single-cpu run time 0 h 0 min 0 sec 22 millisec
 multi-cpu run time 0 h 0 min 4 sec 10 millisec


Cell->Sample table written to "testData/BD_results/Rustody_S1_index/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 4 sec 948milli sec

```

In version 1.2.2 quantify_rhapsody_multi found 5 more cells and 48 more genes using the fasta data. So there is more work needed until this index does do it's job! If we consider that the genomic index also identifies 74 unspliced transcripts is identifies 122 genes less than the fasta based index. A more detailed analysis has to follow. Does the index idex the correct data?

The program by default creates transcipt specific indices. Which is kind of counter intuitive and could be changed later on if necessary.

And what about a genomic index combining Mus_musculus.GRCm39.dna.toplevel.fa.gz and Mus_musculus.GRCm39.104.subset.gtf.gz - I need to check what this subset does actually mean...


```

target/release/quantify_rhapsody_multi -r testData/cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f testData/cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o testData/BD_results/Rustody_S1_index_genomic -s mouse mouse/GRCm39 -e testData/addOn.fa -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96

Writing outfiles ...
filtering cells
Dropping cell with too little counts (n=55)
65 cells have passed the cutoff of 200 umi counts per cell
writing gene expression
sparse Matrix: 65 cell(s), 295 gene(s) and 3519 entries written to path Ok("testData/BD_results/Rustody_S1_index_genomic/BD_Rhapsody_expression"); 
Writing Antibody counts
sparse Matrix: 65 cell(s), 3 gene(s) and 137 entries written to path Ok("testData/BD_results/Rustody_S1_index_genomic/BD_Rhapsody_antibodies"); 
Writing samples table
dense matrix: 65 cell written

Summary:
total      reads  : 500000 reads
no cell ID reads  : 17436 reads (3.49% of total)
no gene ID reads  : 149388 reads (29.88% of total)
N's or too short  : 152 reads (0.03% of total)
cellular reads    : 482412 reads (96.48% of total)

collected read counts:
expression reads  : 328562 reads (68.11% of cellular)
antibody reads    : 3498 reads (0.73% of cellular)
sample reads      : 964 reads (0.20% of cellular)

reported umi counts:
expression reads  : 143662 reads (29.78% of cellular)
antibody reads    : 2864 reads (0.59% of cellular)
sample reads      : 828 reads (0.17% of cellular)

PCR duplicates or bad cells: 335058 reads (69.45% of cellular)

timings:
   overall run time 0 h 0 min 10 sec 816 millisec
   file-io run time 0 h 0 min 0 sec 745 millisec
single-cpu run time 0 h 0 min 0 sec 22 millisec
 multi-cpu run time 0 h 0 min 4 sec 978 millisec


Cell->Sample table written to "testData/BD_results/Rustody_S1_index_genomic/SampleCounts.tsv"

```

That one only identified 16 unspliced transcripts and therefore only 26 genes less than the fasta based index. There is clearly something incorrect with the test index.


I hope you found my little easter egg: you can combine --expression and --index options in the program call. This way you can have a global index and nevertheless add project specific transcriptions like dtTomato in this specific case.


# Speed comparisons to local BD software installation

That is something one does not need to talk about. THis software is incredibly much faster than the BD pipeline.
But Rustody does also not create a bam file.


# Additional Programs

There are several other programs in this package:

## create_index

Instead of creating the index each time from the fasta files this program can generate the index from a fasta.gz file paired up with a gtf.gz or gff.gz file.
This way intron/exon boundaries can be addressed.

```

./target/release/create_index -f testData/mapperTest/Juan_genes.fa.gz -g testData/mapperTest/Juan_genes.fixed.gtf.gz -o testData/mapperTest/index > testData/mapperTest/index/mRNA.fa

```

```
gtf mode
 total first keys 30546
 total second keys 57140
 total single gene per second key 56241
 total multimapper per second key 899
Writing index version 2
with kmer_len 16
And a total of 100863 data entries
finished in 0 h 0 min 0 sec 322 milli sec
```

```
./target/release/quantify_rhapsody -a testData/MyAbSeqPanel.fasta -i testData/mapperTest/index -r testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz  -o testData/mapperTest/ -s mouse -v v1 -m 30
```

```
Summary:
total      reads  : 500000 reads
no cell ID reads  : 0 reads
no gene ID reads  : 1203 reads
N's or too short  : 3028 reads
cellular reads    : 495769 reads (99.15% of total)
expression reads  : 495723 reads (99.14% of total)
antibody reads    : 37 reads (0.01% of total)
sample tag reads  : 9 reads (0.00% of total)
pcr duplicates    : 320518 reads (64.65% of usable)

quantify_rhapsody finished in 0h 0min 2 sec 28milli sec
```

 ## split2samples 

 This will split the BD Rhapsody fastq files into sample spceific fastq files. This script is older and ~4 times slower in creating just the fastq files when compared to quantifyRhapsody quantifying the data.

```
target/release/split2samples -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_split -s mouse -v v1 -m fastqSplit
```
 ## demux10x 

 This is a small spin off that actually processes 10x single cell data and searches for a set fasta entries. (no test data included)

 This package has gotten it's own reop: https://github.com/stela2502/demux10xHTO


 ## bd_cell_id_2_seq 

 BD Rhapsody cells do get an ID in the results. If you want to get the sequences coding for one cell you can use this program:
 
 ```
target/release/bd_cell_id_2_seq -i 3857748 -v v1
The sequence is:
ACCAAGGAC---TTGGAGGTA---CAACCTCCA
```


## bd_get_single_cell 

This will select only one single cell from the fastq files.

```
target/release/bd_get_single_cell -i 3857747 -v v1 -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_getCell
writing all reads from the cell 3857747
[100/?] ⠁                                                                                                                                                 I found 14 reads for the cell 3857747
```

# Limitations / differences

TCR / BRC analyses are not supported at the moment.
The multiprocessor tool breaks if very view cells are detected. Try using less prcessors if this problem occures. A minimum of num_theads +1 cells is required for the export process to work. 


# specific test cases

## Cd3e should not be in this cell - why is it detected

Version 0.1.0 and likely previouse version do detect Cd3e in this cell, but it should not be expressed.

```
target/release/quantify_rhapsody_multi -r  testData/OneSingleCell.66.R2.fastq.gz -f testData/OneSingleCell.66.R1.fastq.gz  -o testData/output_one_cell -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v2.96

zcat testData/output_one_cell/BD_Rhapsody_expression/features.tsv.gz | grep Cd3e
```

This problem has been fixed. Cd3e is no longer "expressed" in this cell.

And this problem re-surfaced later on after I spent so much time making the mapping more sensitive :-D - good I have a negative control here!


## Zbtb16 expression

Comparing my results to the results from the DB pipeline showed Zbtb16 as the gene with the lowest correlation between the tools.
To debug this I selected two cells from the whole data which one had detected expression in both tools and one was only detected in Rustody.
I had hoped my mapping was more sensitive, but it turned out it was these stragne PolA containing R2 reads (again).

I have updated the PolyA detection in the newest version and also updated the report to be more phony of why a read was filtered out.