[![Rust](https://github.com/stela2502/Rustody/actions/workflows/rust.yml/badge.svg)](https://github.com/stela2502/Rustody/actions/workflows/rust.yml)

# Rustody - a tool to quickly analyze targeted BD-Rhapsody sequencings

This tool replaces the official 7-Bridges BD analysis programs first steps.

The final output from this tool is a sparse matrix of gene expression values, antibody tags and in addition a dense matrix with read counts for the sample reads.

I found the sample table most helpful in the detection of populations containing duplicate cells.

The output from here can easily be read into any single cell analysis package for downstream analysis like Seurat or Scanpy.

You can inspect the state of the program using this [deatiled comparison between the seven bridges BD Rhapsody pipeline and the Rhapsody output here]( ./testData/BD_results/CombinedAnalysis_scanpy.ipynb).

## News

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

To run the test data (a tiny bit of a real dataset):

```
target/release/quantify_rhapsody -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v1
```

Latest results:

```
Summary:
total      reads  : 100000 reads
no cell ID reads  : 16355 reads (16.35% of total)
no gene ID reads  : 1217 reads (1.22% of total)
N's or too short  : 14080 reads (14.08% of total)
cellular reads    : 69565 reads (69.56% of total)
expression reads  : 338 reads (0.49% of cellular)
antibody reads    : 15 reads (0.02% of cellular)
sample   reads    : 0 reads (0.00% of cellular)
mapped reads      : 353 reads (0.51% of cellular)
```


Or test the multi processor version:

```
target/release/quantify_rhapsody_multi -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v1
```

Or to even validate this data you can run a R test script (which requires the R::Seurat package) like that:

```
Rscript Rtest/TestExample.R
```


# Usage

## quantfyRhapsody

The `quantifyRhapsody` program takes several arguments.  The usage can be printed 
from the command line using `quantifyRhapsody -h`.

```
target/release/quantify_rhapsody  -h
Rustody 1.0.0
Stefan L. <stefan.lang@med.lu.se>
Quantifies a DB Rhapsody experiment and creates sparse matrix outfiles. You need quite long R1 and
R2 reads for this! (>70R1 and >70R2 [v1] and 52 bp reads for v2.96 and v2.384)

USAGE:
    quantify_rhapsody.exe [OPTIONS] --reads <READS> --file <FILE> --specie <SPECIE> --outpath <OUTPATH> --expression <EXPRESSION> --antibody <ANTIBODY> --min-umi <MIN_UMI> --version <VERSION>

OPTIONS:
    -a, --antibody <ANTIBODY>          the fasta database containing the antibody tags
    -e, --expression <EXPRESSION>      the fasta database containing the genes
    -f, --file <FILE>                  the input R2 samples file
        --gene-kmers <GENE_KMERS>      minimal sequencing quality [default: 32]
    -h, --help                         Print help information
    -m, --min-umi <MIN_UMI>            the minimum reads per cell (sample + genes + antibody
                                       combined)
        --max-reads <MAX_READS>        Optional: end the analysis after processing <max_reads> cell
                                       fastq entries [default: 18446744073709551615]
        --min-quality <MIN_QUALITY>    minimal sequencing quality [default: 25]
    -o, --outpath <OUTPATH>            the outpath
    -r, --reads <READS>                the input R1 reads file
    -s, --specie <SPECIE>              the specie of the library [mouse, human]
    -v, --version <VERSION>            the version of beads you used v1, v2.96 or v2.384

```

This was both compiled and run under Windows11.


Both R1 and R2 files can also be comma separated lists of fastq files.
As always the same order must be preserved in both lists.

The benefit between running both files separately is that the umis are controlled for both fastq files and therefore a 'better' result is obtained than running the fastq files separately and merging the resulting data.

```
target/release/quantify_rhapsody -r  testData/1e5_mRNA_S1_R1_001.fastq.gz,testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz,testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v1
```


## createIndex

BD Rhapsody now also supports whole genome transcriptomics. For this setting it becomes vital to have a mapping index that is capable of storing the genome information.
In order to do this a 2bit based 8mer numeric initial mapper has been combined with 32bp full length match to build up a 40bp full length match.

The creatIndex produces these indices for quantifyRhapsody. It only needs a geneome fasta file and a genome gtf or gff annotation file:

```
target/release/create_index -h
Rustody 1.0.0
Stefan L. <stefan.lang@med.lu.se>
Just run a test case for speed reproducability and simplicity

USAGE:
    create_index [OPTIONS]

OPTIONS:
    -a, --antibody <ANTIBODY>        [default: testData/MyAbSeqPanel.fasta]
    -f, --file <FILE>                the fasta genome data [default: testData/testGenes.fa.gz]
    -g, --gtf <GTF>                  the gff gene information table [default:
                                     testData/testGenes.gtf.gz]
        --gene-kmers <GENE_KMERS>    thje mapping kmer length [default: 32]
    -h, --help                       Print help information
    -o, --outpath <OUTPATH>          the outpath [default: testData/mapperTest]
    -V, --version                    Print version information
```

```
./target/release/create_index -f testData/mapperTest/Juan_genes.fa.gz -g testData/mapperTest/Juan_genes.fixed.gtf.gz -o testData/mapperTest/index

```
**Runtime**

```
real    50m50,368s
user    8m31,537s
sys     42m18,564s
```

**File size**
```
Will come once the tool creates one :-D
```


# Speed comparisons to local BD software installation

The BD rhapsody software is available for Mac and Linux, not for Windows. Whereas this rust program here also compiles and runs on Windows.

System used for the speed comparisons:

```
Operating System: Ubuntu 22.04.1 LTS
          Kernel: Linux 5.15.0-58-generic
model name      : AMD Ryzen 5 3600X 6-Core Processor
total mem       : 64Gb
storage         : Samsung SSD 970 EVO Plus 2TB
```

The test dataset consists of 500.000 reads with cell information.

## time for a BD analysis 500k reads S2

```
# in the testData folder:
cwl-runner --singularity 1.11.1.cwl S2_1.11.1.yml

user    16m47,165s
user    16m47,870s
user    16m48,953s
```

OK I assume we get the point ~ 17min.

## Rustody commands any version (up to now)

```
target/release/quantify_rhapsody -r testData/cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f testData/cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o testData/BD_results/Rustody_S1 -s mouse  -e testData/2276_20220531_chang_to_rpl36a_amplicons.fasta -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96 --gene-kmers 16

target/release/quantify_rhapsody -r testData/cells.1.Rhapsody_SV_index2_S2_R1_001.fastq.gz -f testData/cells.1.Rhapsody_SV_index2_S2_R2_001.fastq.gz -o testData/BD_results/Rustody_S2 -s mouse  -e testData/2276_20220531_chang_to_rpl36a_amplicons.fasta -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96 --gene-kmers 16
```

## time for a Rustody analysis 500k reads S2

```
 time target/release/quantify_rhapsody -r testData/cells.1.Rhapsody_SV_index2_S2_R1_001.fastq.gz -f testData/cells.1.Rhapsody_SV_index2_S2_R2_001.fastq.gz -o Rustody_S2 -s mouse  -e testData/2276_20220531_chang_to_rpl36a_amplicons.fasta -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96 --gene-kmers 16

init models
the log file: Mapping_log.txt
   0.50 mio reads (99.35% usable; 76.08% PCR dupl. [usable] [37.52% for last batch])                                               

Writing outfiles ...
sparse Matrix: 34 cell(s), 450 gene(s) and 4348 entries written (88 cells too view umis) to path Ok("Rustody_S2/BD_Rhapsody_expression"); 
sparse Matrix: 34 cell(s), 4 gene(s) and 20 entries written (88 cells too view umis) to path Ok("Rustody_S2/BD_Rhapsody_antibodies"); 
dense matrix: 34 cell written - 88 cells too view umis

Summary:
total      reads  : 500000 reads
no cell ID reads  : 0 reads
no gene ID reads  : 3568 reads
N's or too short  : 3228 reads
cellular reads    : 493204 reads (98.64% of total)
expression reads  : 493237 reads (98.65% of total)
antibody reads    : 74 reads (0.01% of total)
sample tag reads  : 34 reads (0.01% of total)
pcr duplicates    : 375215 reads (76.08% of usable)

Cell->Sample table written to "Rustody_S2/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 1 sec 208milli sec


real    0m1,213s
user    0m1,197s
sys     0m0,016s
```

### version 1.2.0

```
time ../target/release/quantify_rhapsody -r cells.1.Rhapsody_SV_index2_S2_R1_001.fastq.gz -f cells.1.Rhapsody_SV_index2_S2_R2_001.fastq.gz -o Rustody_S2 -s mouse  -e 2276_20220531_chang_to_rpl36a_amplicons.fasta -a MyAbSeqPanel.fasta -m 200 -v v2.96 --gene-kmers 16
Analysis will stop after having processed 18446744073709551615 fastq entries containing a cell info

init models
the log file: Mapping_log.txt
   0.50 mio reads (99.34% usable; 77.13% PCR dupl. [usable] [37.54% for last batch])                                                                                                                               

Writing outfiles ...
sparse Matrix: 34 cell(s), 283 gene(s) and 2706 entries written (87 cells too view umis) to path Ok("Rustody_S2/BD_Rhapsody_expression"); 
sparse Matrix: 34 cell(s), 2 gene(s) and 2 entries written (87 cells too view umis) to path Ok("Rustody_S2/BD_Rhapsody_antibodies"); 
dense matrix: 34 cell written - 87 cells too view umis

Summary:
total      reads  : 500000 reads
no cell ID reads  : 0 reads
no gene ID reads  : 10066 reads
N's or too short  : 3228 reads
cellular reads    : 486706 reads (97.34% of total)
expression reads  : 486680 reads (97.34% of total)
antibody reads    : 2 reads (0.00% of total)
sample tag reads  : 24 reads (0.00% of total)
pcr duplicates    : 375396 reads (77.13% of usable)

Cell->Sample table written to "Rustody_S2/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 2 sec 481milli sec


real  0m2,486s
user  0m2,482s
sys 0m0,004s
```

### version 1.2.1

```
../target/release/quantify_rhapsody -r cells.1.Rhapsody_SV_index2_S2_R1_001.fastq.gz -f cells.1.Rhapsody_SV_index2_S2_R2_001.fastq.gz -o Rustody_S2 -s mouse  -e 2276_20220531_chang_to_rpl36a_amplicons.fasta -a MyAbSeqPanel.fasta -m 200 -v v2.96 --gene-kmers 16
Analysis will stop after having processed 18446744073709551615 fastq entries containing a cell info

init models
the log file: Mapping_log.txt
   0.50 mio reads (99.35% usable; 77.16% PCR dupl. [usable] [38.04% for last batch])

Writing outfiles ...
sparse Matrix: 34 cell(s), 286 gene(s) and 2739 entries written (87 cells too view umis) to path Ok("Rustody_S2\\BD_Rhapsody_expression");

sparse Matrix: 34 cell(s), 3 gene(s) and 29 entries written (87 cells too view umis) to path Ok("Rustody_S2\\BD_Rhapsody_antibodies");
dense matrix: 34 cell written - 87 cells too view umis

Summary:
total      reads  : 500000 reads
no cell ID reads  : 0 reads
no gene ID reads  : 3723 reads
N's or too short  : 3228 reads
cellular reads    : 493049 reads (98.61% of total)
expression reads  : 492391 reads (98.48% of total)
antibody reads    : 80 reads (0.02% of total)
sample tag reads  : 578 reads (0.12% of total)
pcr duplicates    : 380424 reads (77.16% of usable)

Cell->Sample table written to "Rustody_S2/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 2 sec 71milli sec
```

But if I use 32 bp for the kmers instead of 16 - I loose almost halve of the antibody reads: 54 instead of 80. That is a bug and I need to track that down. At the same time I gain 9 expression reads and in total 6 reported genes. Need to really look into that!

So \> 4 sec. I do not think it makes sense to calculate the difference here.

The two results from this small example are compared [here]( ./testData/BD_results/CombinedAnalysis_scanpy.ipynb).

And the performance on a real dataset?

## Rustody on S2

The full fastq file for this data has in total 781.213.661 reads.

```
Writing outfiles ...
sparse Matrix: 18925 cell(s) and 432 gene(s) and 1903979 entries written (460277 cells too view umis) to path Ok("Sample2_Rustody/BD_Rhapsody_expression");
sparse Matrix: 18925 cell(s) and 4 gene(s) and 56573 entries written (460277 cells too view umis) to path Ok("Sample2_Rustody/BD_Rhapsody_antibodies");
dense matrix: 18925 cell written - 460277 cells too view umis

Summary:
total      reads  : 781213661 reads
no cell ID reads  : 227281056 reads
bad quality       : 6436860 reads
N's or too short  : 96963249 reads
cellular reads    : 456969356 reads (58.49% of total)
expression reads  : 334600069 reads (42.83% of total)
antibody reads    : 100959465 reads (12.92% of total)
sample tag reads  : 21599984 reads (2.76% of total)
pcr duplicates    : 384288196 reads (84.09% of usable)

Cell->Sample table written to "Sample2_Rustody/SampleCounts.tsv"

quantify_rhapsody finished in 1h 36min 8 sec 767milli sec
```

### Outfiles

```
5,5M feb 10 12:12 ./Sample2_Rustody/BD_Rhapsody_expression/matrix.mtx.gz
2,3K feb 10 12:12 ./Sample2_Rustody/BD_Rhapsody_expression/features.tsv.gz
 55K feb 10 12:12 ./Sample2_Rustody/BD_Rhapsody_expression/barcodes.tsv.gz
219K feb 10 12:12 ./Sample2_Rustody/BD_Rhapsody_antibodies/matrix.mtx.gz
  85 feb 10 12:12 ./Sample2_Rustody/BD_Rhapsody_antibodies/features.tsv.gz
 55K feb 10 12:12 ./Sample2_Rustody/BD_Rhapsody_antibodies/barcodes.tsv.gz
1,1M feb 10 12:12 ./Sample2_Rustody/SampleCounts.tsv
72K feb 10 12:13 ./Sample2_Rustody/Mapping_log.txt
```

Little less than 2h. Let's check how much time BD's version does need...



## And BD software for the S2 sample


Repeated runs on my desktop did faile after up to 38h runtime.
Therefore I report the seven bridges run time here: 7 hours, 39 minutes.
That is significantly faster than the 38h on my system, but it nevertheless is ~4x slower than my Rhapsody implementation on a single core.

Why single core? The system IO does affect the Rust implementation. 
The splitting of large gzipped fastq files can easiliy take more than halve an hour and  given the fast processing time
the additional work to implement mutiprocessor capabilities for Rustody seams unnecessary.


## The last run on my desktop:

```
user    2283m4,866s
```
Maximum memory requirement over this time: ~ 22.03 Gb.

More than 38 hours the BD system stopped with an error.

# Additional Programs

There are several other programs in this package:

## create_index

Instead of creating the index each time from the fasta files this program can generate the index from a fasta.gz file paired up with a gtf.gz or gff.gz file.
This way intron/exon boundaries can be addressed.

```

./target/release/create_index -f testData/mapperTest/Juan_genes.fa.gz -g testData/mapperTest/Juan_genes.fixed.gtf.gz -o testData/mapperTest/index --gene-kmers 16 > testData/mapperTest/index/mRNA.fa

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

 ```
 target/release/demux10x  -h
Rustody 0.1.0
Stefan L. <stefan.lang@med.lu.se>
Split a pair of BD rhapsody fastq files (R1 and R2) into sample specific fastq pairs

USAGE:
    demux10x --reads <READS> --file <FILE> --bc <BC> --outpath <OUTPATH> --min-umi <MIN_UMI>

OPTIONS:
    -b, --bc <BC>              the barcodes table name<tab>bc
    -f, --file <FILE>          the input R2 genes file
    -h, --help                 Print help information
    -m, --min-umi <MIN_UMI>    the minimum reads (sample + genes + antybody combined)
    -o, --outpath <OUTPATH>    the outpath
    -r, --reads <READS>        the input R1 reads file
    -V, --version              Print version information
```

 ## bd_cell_id_2_seq 

 BD Rhapsody cells do get an ID in the results. If you want to get the sequences coding for one cell you can use this program:
 
 ```
target/release/bd_cell_id_2_seq -i 3857748 -v v1
The sequence is:
Ok("ACCAAGGAC")
Ok("TTGGAGGTA")
Ok("AATTCGGCG")
```


## bd_get_single_cell 

This will select only one single cell from the fastq files.

```
target/release/bd_get_single_cell -i 3857748 -v v1 -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_getCell
writing all reads from the cell 3857748
[100/?] ⠁                                                                                                                                                                                                  I found 14 reads for the cell 3857748
```

# Limitations / differences

TCR / BRC analyses are not supported at the moment.
The multiprocessor tool breaks if very view cells are detected. Try using less prcessors is this problem occures. A minimum of num_theads +1 cell is required for the export process to work. 


# specific test cases

## Cd3e should not be in this cell - why is it detected

Version 0.1.0 and likely previouse version do detect Cd3e in this cell, but it should not be expressed.

```
target/release/quantify_rhapsody_mulit -r  testData/OneSingleCell.66.R2.fastq.gz -f testData/OneSingleCell.66.R1.fastq.gz  -o testData/output_one_cell -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v2.96

zcat testData/output_one_cell/BD_Rhapsody_expression/features.tsv.gz | grep Cd3e
```

```
init models
the log file: Mapping_log.txt
After indexing all fastq files we have the following index:
I have 7855 kmers and 467 genes
Writing index version 3
with kmer_len 32
And a total of 8063 data entries

Parsing file pair 1

I am using 12 cpus
Starting with data collection
   0.03 mio reads (99.98% with cell info, 89.64% with gene match

Writing outfiles ...
filtering cells
Dropping cell with too little counts (n=0)
sparse Matrix: 1 cell(s), 136 gene(s) and 136 entries written (0 cells too view umis) to path Ok("testData/output_one_cell/BD_Rhapsody_expression");
Writing Antibody counts
sparse Matrix: 1 cell(s), 3 gene(s) and 3 entries written (0 cells too view umis) to path Ok("testData/output_one_cell/BD_Rhapsody_antibodies");
Writing GeneEXpression counts
dense matrix: 1 cell written

Summary:
total      reads  : 32972 reads
no cell ID reads  : 0 reads (0.00% of total)
no gene ID reads  : 3411 reads (10.35% of total)
N's or too short  : 6 reads (0.02% of total)
cellular reads    : 32966 reads (99.98% of total)
expression reads  : 24125 reads (73.18% of cellular)
antibody reads    : 4275 reads (12.97% of cellular)
sample   reads    : 1155 reads (3.50% of cellular)
unique reads      : 29555 reads (89.65% of cellular)

pca duplicates or bad cells: 3411 reads (10.35% of cellular)

timings:
   overall run time 0 h 0 min 0 sec 90 millisec
   file-io run time 0 h 0 min 0 sec 47 millisec
single-cpu run time 0 h 0 min 0 sec 4 millisec
 multi-cpu run time 0 h 0 min 0 sec 24 millisec


Cell->Sample table written to "testData/output_one_cell/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 0 sec 91milli sec
```

This should not return a match. As actuall Cd3e is not detected in 1.0.0.

But even in 1.0.01 there are a view cells where the expression of that gene is detected.

```
target/release/quantify_rhapsody -f  testData/OneSingleCell.11521.R1.fastq.gz -r testData/OneSingleCell.11521.R2.fastq.gz  -o testData/output_one_cell -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v2.96

zcat testData/output_one_cell/BD_Rhapsody_expression/features.tsv.gz | grep Cd3e
```

But this match can obviousely not be blamed on this software as the R1 read that causes this match is a 100% match using BLAST:
```
Query: None Query ID: lcl|Query_19287 Length: 71


>Mus musculus strain C57BL/6J chromosome 9, GRCm39
Sequence ID: NC_000075.7 Length: 124359700
Range 1: 44910123 to 44910193

Score:132 bits(71), Expect:2e-29, 
Identities:71/71(100%),  Gaps:0/71(0%), Strand: Plus/Minus

Query  1         ACAGGTCCTGCCCCATTTATAGATCCTGGCCCAGCCCCTGCCACAGGTGCCTCTCCAGAT  60
                 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  44910193  ACAGGTCCTGCCCCATTTATAGATCCTGGCCCAGCCCCTGCCACAGGTGCCTCTCCAGAT  44910134

Query  61        TTCCCCTTAGA  71
                 |||||||||||
Sbjct  44910133  TTCCCCTTAGA  44910123
```
