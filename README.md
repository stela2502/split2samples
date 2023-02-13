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


Repeated runs on my desktop did faile after up tpo 38h of calculations.
Therefore I report the seven bridges run time here: 7 hours, 39 minutes.
That is significantly faster than the 38h on my system, but it nevertheless is ~4x slower than my Rhapsody implementation.


## The last run on my desktop:

```
user    2283m4,866s
```
Maximum memory requirement over this time: ~ 22.03 Gb.

More than 38 hours the BD system stopped with this error:
<details>
    <summary>

```
40 (function(){return ((JSON.parse(self[0].contents).max_count));})()
stdout was: ''
stderr was: 'evalmachine.<anonymous>:40
(function(){return ((JSON.parse(self[0].contents).max_count));})()
```

</summary>

```
Running RSEC on Cxcr6
Traceback (most recent call last):
  File "/opt/conda/bin/mist_annotate_molecules.py", line 33, in <module>
    sys.exit(load_entry_point('mist==1.11.1', 'console_scripts', 'mist_annotate_molecules.py')())
  File "src/mist/apps/AnnotateMolecules.py", line 172, in mist.apps.AnnotateMolecules.main
  File "src/mist/apps/utils.py", line 349, in mist.apps.utils.node_timer.node_timer_inner
  File "src/mist/apps/utils.py", line 350, in mist.apps.utils.node_timer.node_timer_inner
  File "src/mist/apps/AnnotateMolecules.py", line 76, in mist.apps.AnnotateMolecules.annotate_molecules
  File "src/mist/apps/AnnotateMolecules.py", line 77, in mist.apps.AnnotateMolecules.annotate_molecules
  File "src/mist/apps/AnnotateMolecules.py", line 78, in mist.apps.AnnotateMolecules.annotate_molecules
  File "src/mist/apps/AnnotateMolecules.py", line 105, in genexpr
  File "/opt/conda/lib/python3.9/collections/__init__.py", line 593, in __init__
    self.update(iterable, **kwds)
  File "/opt/conda/lib/python3.9/collections/__init__.py", line 679, in update
    _count_elements(self, iterable)
  File "src/mist/apps/AnnotateMolecules.py", line 105, in genexpr
  File "/opt/conda/lib/python3.9/csv.py", line 111, in __next__
    row = next(self.reader)
  File "/opt/conda/lib/python3.9/gzip.py", line 313, in read1
    return self._buffer.read1(size)
  File "/opt/conda/lib/python3.9/_compression.py", line 68, in readinto
    data = self.read(len(byte_view))
  File "/opt/conda/lib/python3.9/gzip.py", line 478, in read
    self._read_eof()
  File "/opt/conda/lib/python3.9/gzip.py", line 524, in _read_eof
    raise BadGzipFile("CRC check failed %s != %s" % (hex(crc32),
gzip.BadGzipFile: CRC check failed 0xcb91dc88 != 0x6039edc1
WARNING [job AnnotateMolecules_4] exited with status: 1
ERROR Expecting value: line 1 column 1 (char 0)
script was:
01 "use strict";
02 var inputs = {
03     "AbSeq_UMI": null,
04     "Run_Metadata": {
05         "location": "file:///mnt/data2/tmp/l646vxf_/run_metadata.json",
06         "basename": "run_metadata.json",
07         "nameroot": "run_metadata",
08         "nameext": ".json",
09         "class": "File",
10         "checksum": "sha1$98da781bea38bbfdf942b84dedd3306ecfc1528b",
11         "size": 2177,
12         "http://commonwl.org/cwltool#generation": 0,
13         "path": "/var/lib/cwl/stg27313649-e54e-40d0-9077-92ec21a0c264/run_metadata.json",
14         "dirname": "/var/lib/cwl/stg27313649-e54e-40d0-9077-92ec21a0c264"
15     },
16     "Use_DBEC": null,
17     "Valids": {
18         "location": "file:///mnt/data2/tmp/s_mh_f4k/Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv.gz",
19         "basename": "Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv.gz",
20         "nameroot": "Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv",
21         "nameext": ".gz",
22         "class": "File",
23         "checksum": "sha1$b86533ea9a444b9c629c0ed02b240b79b7d2039f",
24         "size": 257705212,
25         "http://commonwl.org/cwltool#generation": 0,
26         "path": "/var/lib/cwl/stgf51c5461-bd47-462c-8559-cb5a3bd953f1/Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv.gz",
27         "dirname": "/var/lib/cwl/stgf51c5461-bd47-462c-8559-cb5a3bd953f1"
28     }
29 };
30 var self = [];
31 var runtime = {
32     "cores": 1,
33     "ram": 32000,
34     "tmpdirSize": 1024,
35     "outdirSize": 1024,
36     "exitCode": 1,
37     "tmpdir": "/tmp",
38     "outdir": "/nIvmVh"
39 };
40 (function(){return ((JSON.parse(self[0].contents).max_count));})()
stdout was: ''
stderr was: 'evalmachine.<anonymous>:40
(function(){return ((JSON.parse(self[0].contents).max_count));})()
                                        ^

TypeError: Cannot read property 'contents' of undefined
    at evalmachine.<anonymous>:40:41
    at evalmachine.<anonymous>:40:65
    at Script.runInContext (vm.js:130:18)
    at Script.runInNewContext (vm.js:135:17)
    at Object.runInNewContext (vm.js:302:38)
    at Socket.<anonymous> ([eval]:11:57)
    at Socket.emit (events.js:314:20)
    at addChunk (_stream_readable.js:297:12)
    at readableAddChunk (_stream_readable.js:268:11)
    at Socket.Readable.push (_stream_readable.js:213:10)'

Traceback (most recent call last):
  File "/usr/lib/python3/dist-packages/cwltool/sandboxjs.py", line 384, in execjs
    return cast(CWLOutputType, json.loads(stdout))
  File "/usr/lib/python3.10/json/__init__.py", line 346, in loads
    return _default_decoder.decode(s)
  File "/usr/lib/python3.10/json/decoder.py", line 337, in decode
    obj, end = self.raw_decode(s, idx=_w(s, 0).end())
  File "/usr/lib/python3.10/json/decoder.py", line 355, in raw_decode
    raise JSONDecodeError("Expecting value", s, err.value) from None
json.decoder.JSONDecodeError: Expecting value: line 1 column 1 (char 0)

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/usr/lib/python3/dist-packages/cwltool/expression.py", line 393, in do_eval
    return interpolate(
  File "/usr/lib/python3/dist-packages/cwltool/expression.py", line 297, in interpolate
    e = evaluator(
  File "/usr/lib/python3/dist-packages/cwltool/expression.py", line 232, in evaluator
    return execjs(
  File "/usr/lib/python3/dist-packages/cwltool/sandboxjs.py", line 386, in execjs
    raise JavascriptException(
cwltool.sandboxjs.JavascriptException: Expecting value: line 1 column 1 (char 0)
script was:
01 "use strict";
02 var inputs = {
03     "AbSeq_UMI": null,
04     "Run_Metadata": {
05         "location": "file:///mnt/data2/tmp/l646vxf_/run_metadata.json",
06         "basename": "run_metadata.json",
07         "nameroot": "run_metadata",
08         "nameext": ".json",
09         "class": "File",
10         "checksum": "sha1$98da781bea38bbfdf942b84dedd3306ecfc1528b",
11         "size": 2177,
12         "http://commonwl.org/cwltool#generation": 0,
13         "path": "/var/lib/cwl/stg27313649-e54e-40d0-9077-92ec21a0c264/run_metadata.json",
14         "dirname": "/var/lib/cwl/stg27313649-e54e-40d0-9077-92ec21a0c264"
15     },
16     "Use_DBEC": null,
17     "Valids": {
18         "location": "file:///mnt/data2/tmp/s_mh_f4k/Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv.gz",
19         "basename": "Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv.gz",
20         "nameroot": "Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv",
21         "nameext": ".gz",
22         "class": "File",
23         "checksum": "sha1$b86533ea9a444b9c629c0ed02b240b79b7d2039f",
24         "size": 257705212,
25         "http://commonwl.org/cwltool#generation": 0,
26         "path": "/var/lib/cwl/stgf51c5461-bd47-462c-8559-cb5a3bd953f1/Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv.gz",
27         "dirname": "/var/lib/cwl/stgf51c5461-bd47-462c-8559-cb5a3bd953f1"
28     }
29 };
30 var self = [];
31 var runtime = {
32     "cores": 1,
33     "ram": 32000,
34     "tmpdirSize": 1024,
35     "outdirSize": 1024,
36     "exitCode": 1,
37     "tmpdir": "/tmp",
38     "outdir": "/nIvmVh"
39 };
40 (function(){return ((JSON.parse(self[0].contents).max_count));})()
stdout was: ''
stderr was: 'evalmachine.<anonymous>:40
(function(){return ((JSON.parse(self[0].contents).max_count));})()
                                        ^

TypeError: Cannot read property 'contents' of undefined
    at evalmachine.<anonymous>:40:41
    at evalmachine.<anonymous>:40:65
    at Script.runInContext (vm.js:130:18)
    at Script.runInNewContext (vm.js:135:17)
    at Object.runInNewContext (vm.js:302:38)
    at Socket.<anonymous> ([eval]:11:57)
    at Socket.emit (events.js:314:20)
    at addChunk (_stream_readable.js:297:12)
    at readableAddChunk (_stream_readable.js:268:11)
    at Socket.Readable.push (_stream_readable.js:213:10)'

ERROR [job AnnotateMolecules_4] Job error:
('Error collecting output for parameter \'Max_Count\': 1.11.1.cwl:231:25: Expression evaluation error:\n1.11.1.cwl:231:25: Expecting value: line 1 column 1 (char 0)\n1.11.1.cwl:231:25: script was:\n1.11.1.cwl:231:25: 01 "use strict";\n1.11.1.cwl:231:25: 02 var inputs = {\n1.11.1.cwl:231:25: 03     "AbSeq_UMI": null,\n1.11.1.cwl:231:25: 04     "Run_Metadata": {\n1.11.1.cwl:231:25: 05         "location": "file:///mnt/data2/tmp/l646vxf_/run_metadata.json",\n1.11.1.cwl:231:25: 06         "basename": "run_metadata.json",\n1.11.1.cwl:231:25: 07         "nameroot": "run_metadata",\n1.11.1.cwl:231:25: 08         "nameext": ".json",\n1.11.1.cwl:231:25: 09         "class": "File",\n1.11.1.cwl:231:25: 10         "checksum": "sha1$98da781bea38bbfdf942b84dedd3306ecfc1528b",\n1.11.1.cwl:231:25: 11         "size": 2177,\n1.11.1.cwl:231:25: 12         "http://commonwl.org/cwltool#generation": 0,\n1.11.1.cwl:231:25: 13         "path": "/var/lib/cwl/stg27313649-e54e-40d0-9077-92ec21a0c264/run_metadata.json",\n1.11.1.cwl:231:25: 14         "dirname": "/var/lib/cwl/stg27313649-e54e-40d0-9077-92ec21a0c264"\n1.11.1.cwl:231:25: 15     },\n1.11.1.cwl:231:25: 16     "Use_DBEC": null,\n1.11.1.cwl:231:25: 17     "Valids": {\n1.11.1.cwl:231:25: 18         "location": "file:///mnt/data2/tmp/s_mh_f4k/Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv.gz",\n1.11.1.cwl:231:25: 19         "basename": "Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv.gz",\n1.11.1.cwl:231:25: 20         "nameroot": "Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv",\n1.11.1.cwl:231:25: 21         "nameext": ".gz",\n1.11.1.cwl:231:25: 22         "class": "File",\n1.11.1.cwl:231:25: 23         "checksum": "sha1$b86533ea9a444b9c629c0ed02b240b79b7d2039f",\n1.11.1.cwl:231:25: 24         "size": 257705212,\n1.11.1.cwl:231:25: 25         "http://commonwl.org/cwltool#generation": 0,\n1.11.1.cwl:231:25: 26         "path": "/var/lib/cwl/stgf51c5461-bd47-462c-8559-cb5a3bd953f1/Rhapsody_SV_index2_Sorted_Valid_Reads.C.csv.gz",\n1.11.1.cwl:231:25: 27         "dirname": "/var/lib/cwl/stgf51c5461-bd47-462c-8559-cb5a3bd953f1"\n1.11.1.cwl:231:25: 28     }\n1.11.1.cwl:231:25: 29 };\n1.11.1.cwl:231:25: 30 var self = [];\n1.11.1.cwl:231:25: 31 var runtime = {\n1.11.1.cwl:231:25: 32     "cores": 1,\n1.11.1.cwl:231:25: 33     "ram": 32000,\n1.11.1.cwl:231:25: 34     "tmpdirSize": 1024,\n1.11.1.cwl:231:25: 35     "outdirSize": 1024,\n1.11.1.cwl:231:25: 36     "exitCode": 1,\n1.11.1.cwl:231:25: 37     "tmpdir": "/tmp",\n1.11.1.cwl:231:25: 38     "outdir": "/nIvmVh"\n1.11.1.cwl:231:25: 39 };\n1.11.1.cwl:231:25: 40 (function(){return ((JSON.parse(self[0].contents).max_count));})()\n1.11.1.cwl:231:25: stdout was: \'\'\n1.11.1.cwl:231:25: stderr was: \'evalmachine.<anonymous>:40\n1.11.1.cwl:231:25: (function(){return ((JSON.parse(self[0].contents).max_count));})()\n1.11.1.cwl:231:25:                                         ^\n1.11.1.cwl:231:25: \n1.11.1.cwl:231:25: TypeError: Cannot read property \'contents\' of undefined\n1.11.1.cwl:231:25:     at evalmachine.<anonymous>:40:41\n1.11.1.cwl:231:25:     at evalmachine.<anonymous>:40:65\n1.11.1.cwl:231:25:     at Script.runInContext (vm.js:130:18)\n1.11.1.cwl:231:25:     at Script.runInNewContext (vm.js:135:17)\n1.11.1.cwl:231:25:     at Object.runInNewContext (vm.js:302:38)\n1.11.1.cwl:231:25:     at Socket.<anonymous> ([eval]:11:57)\n1.11.1.cwl:231:25:     at Socket.emit (events.js:314:20)\n1.11.1.cwl:231:25:     at addChunk (_stream_readable.js:297:12)\n1.11.1.cwl:231:25:     at readableAddChunk (_stream_readable.js:268:11)\n1.11.1.cwl:231:25:     at Socket.Readable.push (_stream_readable.js:213:10)\'', {})
WARNING [job AnnotateMolecules_4] completed permanentFail
INFO [workflow ] starting step VDJ_Preprocess_Reads_IG
INFO [step VDJ_Preprocess_Reads_IG] start
INFO [workflow VDJ_Preprocess_Reads_IG] start
INFO [workflow VDJ_Preprocess_Reads_IG] starting step VDJ_Trim_Reads_2
INFO [step VDJ_Trim_Reads_2] start
INFO [job VDJ_Trim_Reads_2] /mnt/data2/tmp/c9n6z9hk$ singularity \
    --quiet \
    exec \
    --contain \
    --ipc \
    --cleanenv \
    --userns \
    --home \
    /mnt/data2/tmp/c9n6z9hk:/nIvmVh \
    --bind \
    /mnt/data2/tmp/ntx1seu5:/tmp \
    --pwd \
    /nIvmVh \
    /mnt/data2/RhapsodyTest/S2/bdgenomics_rhapsody:1.11.1.sif \
    VDJ_Trim_Reads.sh
Number of arguments is not 1.
    Usage: bash ./VDJ_Trim_Reads.sh <fastq file>
        fastq file: a fastq file that will undergo trimming
    
INFO [job VDJ_Trim_Reads_2] completed success
INFO [step VDJ_Trim_Reads_2] completed success
INFO [workflow VDJ_Preprocess_Reads_IG] starting step VDJ_num_splits_2
INFO [step VDJ_num_splits_2] start
INFO [step VDJ_num_splits_2] completed success
INFO [workflow VDJ_Preprocess_Reads_IG] starting step VDJ_RSEC_Reads_2
INFO [step VDJ_RSEC_Reads_2] start
INFO [job VDJ_RSEC_Reads_2] /mnt/data2/tmp/7ztzjmai$ singularity \
    --quiet \
    exec \
    --contain \
    --ipc \
    --cleanenv \
    --userns \
    --home \
    /mnt/data2/tmp/7ztzjmai:/nIvmVh \
    --bind \
    /mnt/data2/tmp/4axkzt_m:/tmp \
    --pwd \
    /nIvmVh \
    /mnt/data2/RhapsodyTest/S2/bdgenomics_rhapsody:1.11.1.sif \
    mist_vdj_rsec_reads.py \
    --num-splits \
    46
Running with options:
{'num_splits': 46, 'vdj_valid_reads': None}
Beginning Vdj Rsec Reads at 2023-02-12 08:20:03
No valid reads file for this chain type
Completed Vdj Rsec Reads at 2023-02-12 08:20:03
INFO [job VDJ_RSEC_Reads_2] completed success
INFO [step VDJ_RSEC_Reads_2] completed success
INFO [workflow VDJ_Preprocess_Reads_IG] completed success
INFO [step VDJ_Preprocess_Reads_IG] completed success
INFO [workflow ] starting step VDJ_Assemble_and_Annotate_Contigs_IG
INFO [step VDJ_Assemble_and_Annotate_Contigs_IG] start
INFO [workflow VDJ_Assemble_and_Annotate_Contigs_IG] start
INFO [workflow VDJ_Assemble_and_Annotate_Contigs_IG] starting step VDJ_Assemble_and_Annotate_Contigs_IG_2
WARNING [job step VDJ_Assemble_and_Annotate_Contigs_IG_2] Notice: scattering over empty input in 'RSEC_Reads_Fastq'.  All outputs will be empty.
INFO [step VDJ_Assemble_and_Annotate_Contigs_IG_2] completed success
INFO [workflow VDJ_Assemble_and_Annotate_Contigs_IG] completed success
INFO [step VDJ_Assemble_and_Annotate_Contigs_IG] completed success
INFO [workflow ] starting step VDJ_GatherIGCalls
INFO [step VDJ_GatherIGCalls] start
INFO [workflow VDJ_GatherIGCalls] start
INFO [workflow VDJ_GatherIGCalls] starting step VDJ_GatherCalls
INFO [step VDJ_GatherCalls] start
INFO Using local copy of Singularity image found in /mnt/data2/RhapsodyTest/S2
INFO [job VDJ_GatherCalls] /mnt/data2/tmp/w3305mtb$ singularity \
    --quiet \
    exec \
    --contain \
    --ipc \
    --cleanenv \
    --userns \
    --home \
    /mnt/data2/tmp/w3305mtb:/nIvmVh \
    --bind \
    /mnt/data2/tmp/7o52lagt:/tmp \
    --pwd \
    /nIvmVh \
    /mnt/data2/RhapsodyTest/S2/bdgenomics_rhapsody:1.11.1.sif \
    /bin/sh \
    -c \
    echo "No outputs from PyIR detected in VDJ_GatherCalls"
No outputs from PyIR detected in VDJ_GatherCalls
INFO [job VDJ_GatherCalls] completed success
INFO [step VDJ_GatherCalls] completed success
INFO [workflow VDJ_GatherIGCalls] completed success
INFO [step VDJ_GatherIGCalls] completed success
INFO [workflow ] starting step VDJ_Assemble_and_Annotate_Contigs_TCR
INFO [step VDJ_Assemble_and_Annotate_Contigs_TCR] start
INFO [workflow VDJ_Assemble_and_Annotate_Contigs_TCR] start
INFO [workflow VDJ_Assemble_and_Annotate_Contigs_TCR] starting step VDJ_Assemble_and_Annotate_Contigs_TCR_2
WARNING [job step VDJ_Assemble_and_Annotate_Contigs_TCR_2] Notice: scattering over empty input in 'RSEC_Reads_Fastq'.  All outputs will be empty.
INFO [step VDJ_Assemble_and_Annotate_Contigs_TCR_2] completed success
INFO [workflow VDJ_Assemble_and_Annotate_Contigs_TCR] completed success
INFO [step VDJ_Assemble_and_Annotate_Contigs_TCR] completed success
INFO [workflow ] starting step VDJ_GatherTCRCalls
INFO [step VDJ_GatherTCRCalls] start
INFO [workflow VDJ_GatherTCRCalls] start
INFO [workflow VDJ_GatherTCRCalls] starting step VDJ_GatherCalls_2
INFO [step VDJ_GatherCalls_2] start
INFO [job VDJ_GatherCalls_2] /mnt/data2/tmp/g89zad01$ singularity \
    --quiet \
    exec \
    --contain \
    --ipc \
    --cleanenv \
    --userns \
    --home \
    /mnt/data2/tmp/g89zad01:/nIvmVh \
    --bind \
    /mnt/data2/tmp/wsvn1gyd:/tmp \
    --pwd \
    /nIvmVh \
    /mnt/data2/RhapsodyTest/S2/bdgenomics_rhapsody:1.11.1.sif \
    /bin/sh \
    -c \
    echo "No outputs from PyIR detected in VDJ_GatherCalls"
No outputs from PyIR detected in VDJ_GatherCalls
INFO [job VDJ_GatherCalls_2] completed success
INFO [step VDJ_GatherCalls_2] completed success
INFO [workflow VDJ_GatherTCRCalls] completed success
INFO [step VDJ_GatherTCRCalls] completed success
ERROR Workflow cannot make any more progress.
WARNING Final process status is permanentFail

```
</details>


# Additional Programs

There are several other programs in this package:

 1. split2samples will split the BD Rhapsody fastq files into sample spceific fastq files. This script is older and ~4 times slower in creating just the fastq files when compared to quantifyRhapsody quantifying the data.
 2. demux10x is a small spin off that actually processes 10x single cell data and searches for a set fasta entries.
 3. bd_cell_id_2_seq BD Rhapsody cells do get an ID in the results. If you want to get the sequences coding for this cells you can use this program
 4. bd_get_single_cell will select only one single cell from the fastq files.


# Limitations / differences

This program is totally untested and under heavy development.
This is only the first draft - let's see where this heads to.


