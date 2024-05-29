[![Rust](https://github.com/stela2502/Rustody/actions/workflows/rust.yml/badge.svg?branch=main2)](https://github.com/stela2502/Rustody/actions/workflows/rust.yml)

# Rustody - a tool to quickly analyze targeted BD-Rhapsody sequencings

This tool replaces the official 7-Bridges BD analysis programs first steps.

The final output from this tool is a sparse matrix of gene expression values, antibody tags and in addition a dense matrix with read counts for the sample reads.

I found the sample table most helpful in the detection of populations containing duplicate cells.

The output from here can easily be read into any single cell analysis package for downstream analysis like Seurat or Scanpy.

You can inspect the state of the program using this [deatiled comparison between the seven bridges BD Rhapsody pipeline and the Rhapsody output here]( ./testData/BD_results/CombinedAnalysis_scanpy-better_mapping.ipynb).

## News

### 2.2.0

Strange - the whole system was kind of out of order - I had dropped back to version 2.1.3 before I canged my Cargo.toml info.
Anyhow - the exported sam files are now usable - samtools does not complain about them any more.

Features: Targeted remapping of areas around Deletions and Insertions. I have high hopes for this to be working with e.g. ChrM mapping of cell phates. But I am not there at the moment.

### 2.2.1

The new mapping strategy seams to work. Slow but working.
I now can create SAM files, too! New binaries genomic_mapper and index_fq_gene_mapper; the first mapping using the index the second creates. This index does not need a gtf file to create.


### 2.2.0

Implemeting a new mapping strategy. Instead of storing all the 8+32bp combinations for the transcripts, just store a transcript_id and a position in the transcript and store all transcripts in the index.
Work in progress!

### 2.1.0

I finally found fishy reads and had to further improve on the mapping.
Quantify_rhapsody_multi now gives the user access to these new options:

``--min-matches 1`` I one 8bp intitial match (100% identity) and one 32 bp relaxed match of at least 5 tries identify exactly one gene this read will be tagged as coming from that gene. This value is also used to filter if multiple genes are detected, but there the nw-val is more important.
``--highest-humming-val 0.9`` The humming value is used for a quick filtering of really useless reads. It is calculated as absolute difference between all trimers of both the target and the search 32bp fragment. This value is divided by the total length of the comparison. 0.9 is the default value and is very inclusive.
``--highest-nw-val 0.3`` The nw value is a Needleman-Wunsch inspired value. All 32bp fragments that pass the humming test, the initial table of the Needleman Wunsch algorithm is calculated and the final value from that comparison is again divided by the total length of the initial (max) 32bp fragments. In the end both the amount of passing matches as well as the mean nw value of a read will be used to identify the matching gene.
Values above 0.3 will lead to a lot of false postives.

In addition simple DNA fragments are now excluded from the index - like simple repeats or long stretches of one base. They did lead to mapping of really shady reads e.g. containing polA sequences from some random genes (NCBI blastn confirmed).

### 2.0.0

Transposable elements like Sine and Line repeats can now be indexed and mapped - both in BD Rhapsody as well as in 10X expression datasets. Analysis of 10x data is very slow at the moment. Problem is currently (08.02.2024) worked on. 

Mapper has significantly improved: both the false positive as well as false negative rate had improved.
The improvement was possible by using a needleman-wunsch inspired algorithm.
The 32 bp matches are now tolerant to not only bp mismatches, but also insertions and deletions.

Compared to the mere bp replacement matching we e.g. find almost 10x more reads from the Ighm locus in the test data.
And none of the reads I have seen so far looked like a not Ighm transcript. All I checked were also mapped to the Ighm transcripts using NCBI BLAST (I only checked the strange looking ones).

I have added the Ighm reads that were detected with both settings to this repository.


quantify_rhapsody has finally gotten a muti processor upgrade: quantify_rhapsody_multi.
As the old single processer version is no longer actively used the new multi version does now even find more reads.

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
cp target/release/quantify_gene_mapper /usr/bin
cp target/release/bd_cell_id_2_seq /usr/bin
cp target/release/bd_get_single_cell /usr/bin
cp target/release/genomic_mapper /usr/bin
cp target/release/get_n_cell_reads /usr/bin
cp target/release/int_2_seq /usr/bin
cp target/release/check_read_gene_mapper /usr/bin
cp target/release/check_read_fast_mapper /usr/bin
cp target/release/create_gene_mapper_index /usr/bin
cp target/release/create_index /usr/bin
cp target/release/create_index_te /usr/bin
cp target/release/te_analysis /usr/bin
cp target/release/create_gene_mapper_index /usr/bin
``` 

With 10x data we need also the 10x whitelists and therefore we need to set a environment variable.
This is not complicated, but needs admin rights if not working in an docker or singularity environment:

Make sure to compile the source before you run these install scripts. They do not compile the source.

**Linux or other Unix likes:**

This would need sudo right if not working in a docker or singularity image:
```
bash install_resources_unix.sh
```

**Windows (10 or 11)**
```
.\install_resources_windows.ps1 -ExecutionPolicy Bypass
```


Do not forget the --release while building the tool. 
The test case for quantify_rhapsody would finish in 55 sec instead of 3 sec (> x15!)
using a AMD Ryzen 7 5700X processor and a SSD as mass storage.

This info is a little outdated - quantify_rhapsody would need ~ 5 to 6 sec in version 2.0.0, but the \_multi version needs only 1 sec!


## Testing

The test data is a tiny bit of an (unpublished) real data set.
It is included here solely to present this mapping procedure.

To run the test data (a tiny bit of a real dataset):

### the main mappers

```
target/release/quantify_rhapsody_multi -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v1
```
Or the newer version of a mapper:
```
target/release/quantify_gene_mapper -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5_gm -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v1
```

Bummer - the gene_mapper version of my tool is way slower (50x!!) than the fast_mapper one - but the results should be lot more reliable!

```
tis 14 maj 2024 14:21:05 CEST

writing gene expression
sparse Matrix: 73 cell(s), 99 gene(s) and 428 entries written to path Ok("testData/output_1e5/BD_Rhapsody_expression"); 
Writing Antibody counts
sparse Matrix: 73 cell(s), 5 gene(s) and 106 entries written to path Ok("testData/output_1e5/BD_Rhapsody_antibodies"); 
Writing samples table
dense matrix: 73 cell written

Summary:
cellular   reads  : 68556 reads (68.56% of total)
no cell ID reads  : 17366 reads (17.37% of total)
no gene ID reads  : 0 reads (0.00% of total)
filtered   reads  : 14078 reads (14.08% of total)
 ->  multimapper  : 0 reads (0.00% of total)
 -> bad qualiity  : 10696 reads (10.70% of total)
 ->    too short  : 3381 reads (3.38% of total)
 ->          N's  : 1 reads (0.00% of total)

total      reads  : 100000 reads

collected read counts:
expression reads  : 44166 reads (64.42% of cellular)
antibody reads    : 19429 reads (28.34% of cellular)
sample reads      : 862 reads (1.26% of cellular)

reported UMI counts:
expression reads  : 540 UMIs (0.79% of cellular)
antibody reads    : 329 UMIs (0.48% of cellular)
sample reads      : 9 UMIs (0.01% of cellular)

PCR duplicates or bad cells: 67678 reads (98.72% of cellular)

timings:
   overall run time 0 h 0 min 25 sec 518 millisec
   file-io run time 0 h 0 min 0 sec 248 millisec
single-cpu run time 0 h 0 min 0 sec 98 millisec
 multi-cpu run time 0 h 0 min 25 sec 130 millisec


Cell->Sample table written to "testData/output_1e5/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 25 sec 519milli sec
```

With this new gene_mapper model you now get a sam file together with your fasta data:

```
samtools view -hb testData/output_1e5_gm/AlignedReadsOfInterest.sam > testData/output_1e5_gm/AlignedReadsOfInterest.bam
samtools faidx testData/output_1e5_gm/GenesOfInterest.fasta
samtools sort testData/output_1e5_gm/AlignedReadsOfInterest.bam > testData/output_1e5_gm/AlignedReadsOfInterest.sorted.bam
samtools index  testData/output_1e5_gm/AlignedReadsOfInterest.sorted.bam
```


Or want to test the 10x version of the tool (which is a LOT slower!):
Here the most interesting is the sam file that will also be produced in the newest version: ``testData/10x/results/AlignedReadsOfInterest.sam``. I hope this will help to look into chrM cell tagging.
```
./target/release/quantify_gene_mapper -e testData/chrM_GRCm39.primary_assembly.genome.fa.gz -r testData/10x/1k_mouse_kidney_CNIK_3pv3_S1_L004_R1_4m.fastq.gz  -f testData/10x/1k_mouse_kidney_CNIK_3pv3_S1_L004_R2_4m.fastq.gz -o testData/10x/results/ --exp 10x --specie human --min-umi 20 --version "Single Cell 3' v3" --report4genes chrM --highest-nw-val 0.25  --chunk-size 10000
```

[An example analysis of the fast_mapper based analys data is available here:]( ./testData/BD_results/CombinedAnalysis_scanpy_v1.2.1.ipynb).

### the new genomic mapper

The main output from this mapper is a sam file with all mapping sequences.
The idea behind this is of cause to look into mutations. Hope to get this mapper to really find interesting stuf.

```
# ceate the index
target/release/index_fq_gene_mapper -f testData/chrM_GRCm39.primary_assembly.genome.fa.gz --outpath testData/ChrM_MouseIdx
# use the index

./target/release/genomic_mapper -i testData/ChrM_MouseIdx -r testData/10x/1k_mouse_kidney_CNIK_3pv3_S1_L004_R1_4m.fastq.gz  -f testData/10x/1k_mouse_kidney_CNIK_3pv3_S1_L004_R2_4m.fastq.gz -o testData/10x/results/ --exp 10x --specie human --min-umi 20 --version "Single Cell 3' v3" --highest-nw-val 0.25  --chunk-size 10000
```


# Usage

## quantfyRhapsody

The `quantifyRhapsody` program takes several arguments.  The usage can be printed 
from the command line using `quantifyRhapsody -h`.

```
target/release/quantify_rhapsody_multi  -h
rustody 1.2.5
Stefan L. <stefan.lang@med.lu.se>
Quantifies a DB Rhapsody experiment and creates sparse matrix outfiles. You need quite long R1 and
R2 reads for this! (>70R1 and >70R2 \[v1\] and 52 bp reads for v2.96 and v2.384)

USAGE:
    quantify_rhapsody_multi [OPTIONS] --reads <READS> --file <FILE> --specie <SPECIE> --outpath <OUTPATH> --min-umi <MIN_UMI> --version <VERSION>

OPTIONS:
    -a, --antibody <ANTIBODY>          the fasta database containing the antibody tags
        --chunk-size <CHUNK_SIZE>      how many sequences should be analyzed in one chunk [default:
                                       1000000]
    -e, --expression <EXPRESSION>      the fasta database containing the genes
        --exp <EXP>                    this is a BD rhapsody or a 10x expression experiment?
                                       [default: bd]
    -f, --file <FILE>                  the input R2 samples file
        --gene-kmers <GENE_KMERS>      minimal sequencing quality [default: 32]
    -h, --help                         Print help information
    -i, --index <INDEX>                a pre-defined index folder produced by the cerateIndex scipt
    -m, --min-umi <MIN_UMI>            the minimum (UMI) reads per cell (sample + genes + antibody
                                       combined)
        --max-reads <MAX_READS>        Optional: end the analysis after processing <max_reads> cell
                                       fastq entries [default: 18446744073709551615]
        --min-quality <MIN_QUALITY>    minimal sequencing quality [default: 25]
    -n, --num-threads <NUM_THREADS>    how many threads to use to analyze this (default max
                                       available)
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
target/release/quantify_rhapsody_multi -r  testData/1e5_mRNA_S1_R1_001.fastq.gz,testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz,testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 10 -v v1
```


## createIndex

BD Rhapsody now also supports whole genome transcriptomics. For this setting it becomes vital to have a mapping index that is capable of storing the genome information.
In order to do this a 2bit based 8mer numeric initial mapper has been combined with 32bp full length match to build up a 40bp full length match.

The creatIndex produces these indices for quantifyRhapsody. It only needs a geneome fasta file and a genome gtf or gff annotation file:

```
target/release/create_index -h
rustody 1.0.0
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
                                       (transcript_id) [default: gene_name]
    -V, --version                      Print version information

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

At the moment this tool would create a transcript specific index, but report only gene names likely creating a broken index. Use the ``--transcript gene_name`` to prevent the tool from miss-behaving.

```

./target/release/create_gene_mapper_index -g testData/KI270728.1.gtf.gz -f testData/KI270728.1.fa.gz -o testData/mapperTest/index 

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

# Analyze the BD test data Rhap-VDJ-Demo

The data can be found on the the BD public data share:
http://bd-rhapsody-public.s3-website-us-east-1.amazonaws.com/Rhapsody-Demo-Data-Inputs/VDJDemo/


The VDJ part of this data is not handled at the moment using Rustody and currently there is no time plan or idea of how to implement that.

Analyzing this data with the latest version (devel branch commit 63ea23f) using this command:
```
/home/med-sal/git_Projects/Rustody/target/release/quantify_rhapsody_multi -r /mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L004_R1_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L003_R1_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L001_R1_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L002_R1_001.fastq.gz  -f /mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L004_R2_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L002_R2_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L001_R2_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L003_R2_001.fastq.gz -o /mnt/data2/RhapsodyTest/VDJ_v1_example/rustify_testData_result -s human -e /mnt/data2/RhapsodyTest/VDJ_v1_example/BD_Rhapsody_Immune_Response_Panel_Hs.fasta -m 200 -v 'v1' 1>&2 > /mnt/data2/RhapsodyTest/VDJ_v1_example/rustify_testData_result.log
```

I get these measurements:
```
I am using 12 cpus


Writing outfiles ...
filtering cells
Dropping cell with too little counts (n=143511)
2766 cells have passed the cutoff of 200 umi counts per cell.


writing gene expression
sparse Matrix: 2766 cell(s), 367 gene(s) and 414279 entries written to path Ok("/mnt/data2/RhapsodyTest/VDJ_v1_example/rustify_testData_result/BD_Rhapsody_expression"); 
Writing Antibody counts
Writing samples table
dense matrix: 2766 cell written

Summary:
cellular   reads  : 4832305 reads (62.61% of total)
no cell ID reads  : 2371673 reads (30.73% of total)
no gene ID reads  : 0 reads (0.00% of total)
filtered   reads  : 514614 reads (6.67% of total)
 ->  multimapper  : 0 reads (0.00% of total)
 -> bad qualiity  : 511637 reads (6.63% of total)
 ->    too short  : 1009 reads (0.01% of total)
 ->          N's  : 1968 reads (0.03% of total)

total      reads  : 7718592 reads

collected read counts:
expression reads  : 4727781 reads (97.84% of cellular)
antibody reads    : 0 reads (0.00% of cellular)
sample reads      : 0 reads (0.00% of cellular)

reported UMI counts:
expression reads  : 2748802 UMIs (56.88% of cellular)
antibody reads    : 0 UMIs (0.00% of cellular)
sample reads      : 0 UMIs (0.00% of cellular)

PCR duplicates or bad cells: 2083503 reads (43.12% of cellular)

timings:
   overall run time 0 h 1 min 30 sec 621 millisec
   file-io run time 0 h 0 min 16 sec 246 millisec
single-cpu run time 0 h 0 min 1 sec 532 millisec
 multi-cpu run time 0 h 1 min 11 sec 484 millisec


Cell->Sample table written to "/mnt/data2/RhapsodyTest/VDJ_v1_example/rustify_testData_result/SampleCounts.tsv"

quantify_rhapsody finished in 0h 1min 30 sec 662milli sec
```

## New version with updated mapper is way faster - still

The new version has the first time some control over the mapping process:

``--min-matches 1`` I one 8bp intitial match (100% identity) and one 32 bp relaxed match of at least 5 tries identify exactly one gene this read will be tagged as coming from that gene.
``--highest-humming-val 0.9`` The humming value is used for a quick filtering of really useless reads. It is calculated as absolute difference between all trimers of both the target and the search 32bp fragment. This value is divided by the total length of the comparison. 0.9 is the default value and is very inclusive.
``--highest-nw-val 0.5`` The nw value is a Needleman-Wunsch inspired value. All 32bp fragments that pass the humming test, the initial table of the Needleman Wunsch algorithm is calculated and the final value from that comparison is again divided by the total length of the initial (max) 32bp fragments. In the end both the amount of pssing matches as well as the mean nw value of a read will be used to identify the matching gene.

The new version also is more restrictive with index creation. It does not allow for too simple 8bp + 32bp combinations like these:
```
Useless oligo(s) detected! TGTGTGTG, AGTGTGAGTGTGAGCGAGAGGGTGAGTGTGGT for gene CCL19
Useless oligo(s) detected! TTTTTTGT, TTGTTTGTTTTGTTTTGTTTGTTGTTTGTTGT for gene CD9
Useless oligo(s) detected! TATATATT, TATATTTTTAAAATATTTATTTATTTATTTAT for gene CSF2
Useless oligo(s) detected! GTGTGTGT, GTGTGTGTGTGTGTGTGTGTGTGTATGACTAA for gene FASLG
Useless oligo(s) detected! TTTTTTAA, GTCTATGTTTTAAAATAATATGTAAATTTTTC for gene FLT3
Useless oligo(s) detected! AAAAAATA, AAATAAATAAATAAACAAATAAAAAATT for gene IL18
Useless oligo(s) detected! CTTTTTTA, AATATAAAAATGGGTGTTATTTT for gene LAMP1
Useless oligo(s) detected! TTTTTTTA, AGAAAAAAAAGAGAAATGAATAAAGAATCTAC for gene LIF
Useless oligo(s) detected! TTTTTTTT, TAAGAAAAAAAAGAGAAATGAATAAAGAATCT for gene LIF
Useless oligo(s) detected! CAAAAAAA, TTAAATTATTTATTTATGGAGGATGGAGAGAG for gene LTA
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAAAAAAAAAAAAAAAATTTAT for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAAAAAAAATTTATGTATATAA for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAATTTATGTATATAAAAGTTGCA for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAAAAAAAAAAAAATTTATGTA for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAAAAATTTATGTATATAAAAG for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAAAAAAAAAAAAAAAAAATTT for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAAAAAAAAAATTTATGTATAT for gene SLC25A37
Useless oligo(s) detected! AAAAAAAA, AAAAAAAAAAAAATTTATGTATATAAAAGTTG for gene SLC25A37
Useless oligo(s) detected! TGTGTGTT, TTTTCTTTTTCTTTCTTTTTATTTTTTTTGAA for gene TBX21
Useless oligo(s) detected! GTTTTTTC, TTTTTTTTGTTTTGTTTTTTTTTTTTTTTTTT for gene THBS1
Useless oligo(s) detected! TTTTTTTT, GTTTTGTTTTTTTTTTTTTTTTTTTTTGCTTT for gene THBS1
Useless oligo(s) detected! ACACACAC, AACATCACAATGACACACACATCACACACACA for gene TNFSF14
Useless oligo(s) detected! ACACACAT, TACACACACATCACAATGACAAACACAACATT for gene TNFSF14
Useless oligo(s) detected! ACACACAC, ACAACATCACAATGACACACACATCACACACA for gene TNFSF14
Useless oligo(s) detected! ACACACAT, CACACACACATCACAATGACAAACACACAACA for gene TNFSF14
Useless oligo(s) detected! TGTGTGTA, TATATATATATATATGTTTATGTATATATGTG for gene VEGFA
Useless oligo(s) detected! TATATATA, TATATATGTTTATGTATATATGTGATTCTGAT for gene VEGFA
```

You can check that using the map_one_sequence tool as it will also index your fastq files and report the dropped fragments.


For this comparison we set the default value to a little more restrictive values (nw value of 0.5 leads to really useless matches!).

Of casue all paths need to be changed to accomodate your situation.
```
 /home/med-sal/git_Projects/Rustody/target/release/quantify_rhapsody_multi -r /mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L004_R1_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L003_R1_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L001_R1_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L002_R1_001.fastq.gz  -f /mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L004_R2_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L002_R2_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L001_R2_001.fastq.gz,/mnt/data2/RhapsodyTest/VDJ_v1_example/RhapVDJDemo-mRNA_S5_L003_R2_001.fastq.gz -o /mnt/data2/RhapsodyTest/VDJ_v1_example/rustify_testData_result -s human -e /mnt/data2/RhapsodyTest/VDJ_v1_example/BD_Rhapsody_Immune_Response_Panel_Hs.fasta -m 200 -v 'v1' --min-matches 2 --highest-humming-val 0.6 --highest-nw-val 0.2
```

And here comes the result using my (old) 12 core AMD machine.

```
Writing outfiles ...
filtering cells
Dropping cell with too little counts (n=143818)
2767 cells have passed the cutoff of 200 umi counts per cell.


writing gene expression
sparse Matrix: 2767 cell(s), 368 gene(s) and 418656 entries written to path Ok("/mnt/data2/RhapsodyTest/VDJ_v1_example/rustify_testData_result/BD_Rhapsody_expression"); 
Writing Antibody counts
No genes to report on - no data written to path Some("/mnt/data2/RhapsodyTest/VDJ_v1_example/rustify_testData_result/BD_Rhapsody_antibodies")
Writing samples table
dense matrix: 2767 cell written

Summary:
cellular   reads  : 4832305 reads (62.61% of total)
no cell ID reads  : 2371673 reads (30.73% of total)
no gene ID reads  : 0 reads (0.00% of total)
filtered   reads  : 514614 reads (6.67% of total)
 ->  multimapper  : 0 reads (0.00% of total)
 -> bad qualiity  : 511637 reads (6.63% of total)
 ->    too short  : 1009 reads (0.01% of total)
 ->          N's  : 1968 reads (0.03% of total)

total      reads  : 7718592 reads

collected read counts:
expression reads  : 4745007 reads (98.19% of cellular)
antibody reads    : 0 reads (0.00% of cellular)
sample reads      : 0 reads (0.00% of cellular)

reported UMI counts:
expression reads  : 2761260 UMIs (57.14% of cellular)
antibody reads    : 0 UMIs (0.00% of cellular)
sample reads      : 0 UMIs (0.00% of cellular)

PCR duplicates or bad cells: 2071045 reads (42.86% of cellular)

timings:
   overall run time 0 h 0 min 33 sec 738 millisec
   file-io run time 0 h 0 min 15 sec 103 millisec
single-cpu run time 0 h 0 min 1 sec 584 millisec
 multi-cpu run time 0 h 0 min 16 sec 124 millisec


Cell->Sample table written to "/mnt/data2/RhapsodyTest/VDJ_v1_example/rustify_testData_result/SampleCounts.tsv"

quantify_rhapsody finished in 0h 0min 33 sec 739milli sec
```

```
cp target/release/split2samples ~/sens05_home/bin
cp target/release/quantify_rhapsody ~/sens05_home/bin
cp target/release/quantify_rhapsody_multi ~/sens05_home/bin
cp target/release/bd_cell_id_2_seq ~/sens05_home/bin
cp target/release/bd_get_single_cell ~/sens05_home/bin
cp target/release/genomic_mapper ~/sens05_home/bin
cp target/release/get_n_cell_reads ~/sens05_home/bin
cp target/release/int_2_seq ~/sens05_home/bin
cp target/release/create_index ~/sens05_home/bin
cp target/release/create_index_te ~/sens05_home/bin
cp target/release/create_gene_mapper_index ~/sens05_home/bin
cp target/release/te_analysis ~/sens05_home/bin
```

