# Rustody - a tool to quickly analyze targeted BD-Rhapsody sequencings

This tool replaces the official 7-Bridges BD analysis programs first steps.

The final output from this tool is a sparse matrix of gene expression values, antibody tags and in addition a dense matrix with read counts for the sample reads.

I found the sample table most helpful in the detection of populations containing duplicate cells.

The output from here can easily be read into any single cell analysis package for downstream analysis like Seurat or Scanpy.

You can inspect the state of the program using this [deatiled comparison between the seven bridges BD Rhapsody pipeline and the Rhapsody output here]( ./testData/BD_results/CombinedAnalysis_scanpy.ipynb).

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
cp target/release/quantify_rhapsody /usr/bin
cp target/release/bd_cell_id_2_seq /usr/bin
cp target/release/bd_get_single_cell /usr/bin
cp target/release/get_n_cell_reads /usr/bin
``` 

Do not forget the --release while building the tool. 
The test case for quantify_rhapsody would finish in 55 sec instead of 3 sec (> x15!)
using a AMD Ryzen 7 5700X processor and a SSD as mass storage.


## Testing

To run the test data (a tiny bit of a real dataset):

```
target/release/quantify_rhapsody -r  testData/1e5_mRNA_S1_R1_001.fastq.gz -f testData/1e5_mRNA_S1_R2_001.fastq.gz -o testData/output_1e5 -s mouse  -e testData/genes.fasta -a testData/MyAbSeqPanel.fasta -m 30 -v v1
```

Or to even validate this data you can run a R test script (which requires the R::Seurat package) like that:

```
Rscript Rtest/TestExample.R
```


# Usage

The `quantifyRhapsody` program takes several arguments.  The usage can be printed 
from the command line using `quantifyRhapsody -h`.

```
target/release/quantify_rhapsody  -h
Rustody 0.3.4
Stefan L. <stefan.lang@med.lu.se>
Quantifies a DB Rhapsody experiment and creates sparse matrix outfiles. You need quite long R1 and
R2 reads for this! (>70R1 and >70R2 [v1] and 52 bp reads for v2.96 and v2.384)

USAGE:
    quantify_rhapsody [OPTIONS] --reads <READS> --file <FILE> --specie <SPECIE> --outpath <OUTPATH> --expression <EXPRESSION> --antibody <ANTIBODY> --min-umi <MIN_UMI> --version <VERSION>

OPTIONS:
    -a, --antibody <ANTIBODY>        the fasta database containing the antibody tags
    -e, --expression <EXPRESSION>    the fasta database containing the genes
    -f, --file <FILE>                the input R2 samples file
    -h, --help                       Print help information
    -m, --min-umi <MIN_UMI>          the minimum reads per cell (sample + genes + antibody combined)
        --max-reads <MAX_READS>      Optional: end the analysis after processing <max_reads> cell
                                     fastq entries [default: 18446744073709551615]
    -o, --outpath <OUTPATH>          the outpath
    -r, --reads <READS>              the input R1 reads file
    -s, --specie <SPECIE>            the specie of the library [mouse, human]
    -v, --version <VERSION>          the version of beads you used v1, v2.96 or v2.384

```

You see - this is the one I compiled on Windows 11.


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

## time for a BD analysis S2

```
# in the testData folder:
cwl-runner --singularity 1.11.1.cwl S2_1.11.1.yml

user    16m47,165s
user    16m47,870s
user    16m48,953s
```

OK I assume we get the point ~ 17min.

## time for a Rustody analysis S2

```
 time ../target/release/quantify_rhapsody -r cells.1.Rhapsody_SV_index2_S2_R1_001.fastq.gz -f cells.1.Rhapsody_SV_index2_S2_R2_001.fastq.gz -o Rustody_S2 -s mouse  -e 2276_20220531_chang_to_rpl36a_amplicons.fasta -a MyAbSeqPanel.fasta -m 200 -v v2.96

user    0m3,983s
user    0m4,040s
user    0m4,020s
```

So ~ 4 sec. I do not think it makes sense to calculate the difference here.

How does it look for a real dataset?

## Rustody on S2

In total 781.213.661 reads.
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

Little less than 2h. Let's check how much time BD's version doe need...

## And BD software for the S2 sample

```
user    215m25,699s
```

Almost 4 hours. So Rhapsody takes about halve the time to get comparable results.

Ha - but that run was a faliure: ``OSError: [Errno 28] No space left on device``. I had ~ 50Gb on my tmp drive and 357Gb on the data drive. Is that too little?!

Restarted the system 2023 02 11 ~ 9:30

# Additional Programs

There are several other programs in this package:

 1. split2samples will split the BD Rhapsody fastq files into sample spceific fastq files. This script is older and ~4 times slower in creating just the fastq files when compared to quantifyRhapsody quantifying the data.
 2. demux10x is a small spin off that actually processes 10x single cell data and searches for a set fasta entries.
 3. bd_cell_id_2_seq BD Rhapsody cells do get an ID in the results. If you want to get the sequences coding for this cells you can use this program
 4. bd_get_single_cell will select only one single cell from the fastq files.


# Limitations / differences

This program is totally untested and under heavy development.
This is only the first draft - let's see where this heads to.


