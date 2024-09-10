# 2.2.2

The Cigars are now - hopefully standard conform.
Lots of bugs fixed. Speed of genomic_mapper is still far lower than the fast_mapper, but I have more confidence in it's mapping quality.


# 2.2.1

Using 10x data the index length of 400 bp is too short.
It has now been changed to 500 bp for both indices.

Indices from previouse versions will not work any more!

# 2.2.0

The created SAM files now are accepted by samtools.

# 2.1.2

quantify_gene_mapper now creates sam out files for the matches you wanted reported.
This allowes for a simple mappiing against e.g. the ChrM as the mapper simply indexes a fasta entry on the fly.

It is still unclear if end tagged expression data has the capability to get enough genetic information on a mitochodria
to use this it as a cell index, but this is a first step to answer that.

# 2.1.1

Adding a new index type - istead of indexing 32bp fragments the whole fasta entry is now stored as 2bit encoded u8 Vec. This vector is later on used to find the mapping areas for the reads.

This lead to two new binaries:

check_read_gene_mapper: A tool to check a single reads mapping in detail.
quantify_gene_mapper: Map whichever single cell data using the new index.

# 2.1.0

I finally found fishy reads and had to further improve on the mapping.
Quantify_rhapsody_multi now gives the user access to these new options:

``--min-matches 1`` I one 8bp intitial match (100% identity) and one 32 bp relaxed match of at least 5 tries identify exactly one gene this read will be tagged as coming from that gene. This value is also used to filter if multiple genes are detected, but there the nw-val is more important.
``--highest-humming-val 0.9`` The humming value is used for a quick filtering of really useless reads. It is calculated as absolute difference between all trimers of both the target and the search 32bp fragment. This value is divided by the total length of the comparison. 0.9 is the default value and is very inclusive.
``--highest-nw-val 0.3`` The nw value is a Needleman-Wunsch inspired value. All 32bp fragments that pass the humming test, the initial table of the Needleman Wunsch algorithm is calculated and the final value from that comparison is again divided by the total length of the initial (max) 32bp fragments. In the end both the amount of passing matches as well as the mean nw value of a read will be used to identify the matching gene.
Values above 0.3 will lead to a lot of false postives.

# 2.0.0

This is a huge update. The mapper has become even better. Is is almost as fast as before the huge update, but finds lots more genes now.
And I can not find any read (up to now) that would be a false positive.

The genomic analysis should in theory be working now too, as the index has gotten the most updates this time.

There is now a new binary - te_analysis. It can do almost everything the quantify_rhapsody_multi can do, but it does not index fasta files.
The enw binary needs index folders to work. BUT - it can utilize two indices at the same time!

The expression index is matched using fuzzy matching and the te index is match using exact matches. So using this you can e.g. compare two indices
or index both human and mouse genes or you can index genes and transposable elements (TE) at the same time.

The TE's are also collected per cell - so you can analyze a whole e.g. 10x sample using one run. The first version needed ~1h per 100 mio reads is only TE elements from chr1 are checked.
More info when available.

# 1.2.5

Mapper has significantly improved: both the false positive as well as false negative rate had improved.
The improvement was possible by using a needleman-wunsch inspired algorithm.
The 32 bp matches are now tolerant to not only bp mismatches, but also insertions and deletions.

There is no longer a need to exclude polyA containing reads.

Compared to the mere bp replacement matching we e.g. find almost 10x more reads from the Ighm locus in the test data.
And none of the reads I have seen so far looked like a not Ighm transcript. All I checked were also mapped to the Ighm transcripts using NCBI BLAST (I only checked the strange looking ones).

I have added the Ighm reads that were detected with both settings to this repository.


# 1.2.4

PolyA containing R2 reads are now filtered out if a PolyA streatch of at least 15 A's is detected in the last 30 bp of a R2 read. 

# 1.2.3

Analysis of 10x data should now also be possible.


# 1.2.2

It is now possible to create a genome wide index file.
Matching to an genome wide index is horribly slow. I need to identify the problem.


# 1.2.1

Sample table now also contains a n column with total reads over all samples.
This should be used as a weight to qualify the FractionTotal.

Mapping with smaller kmers_site values (like 16 or even 10) seams to have an overall
positive influence of mapping efficiency. This feels a little fishy.
I will look into that.

# 1.2.0

The mapper has changed completely.
Switched from a BTreeMap to a 2bit u16 initial vector with the 2bit u16 entries as index and later on 2bit u64 to match the following sequences.
Missing sequences are 'replaced' by 'A' == 0.

The index can now also be calculated over the whole genome using the new create_index tool.
In addition there is a small helper that converts a up to u128 integer to Nucleotide sequence (int_2_seq).


# 0.2.0

Support for three bd primer types added:
```
    -v, --version <VERSION>          the version of beads you used v1, v2.96 or v2.384
```

# 0.1.1

Thinking about a two step process to get sample specific fastq files.
First run to collect the sample specific cell IDs which are hidden in R1.
Second run to actually split the fastq files based on the cell ids.

# 0.1.0

The first implementation mainly performed by Rob to get me up to speed.
Now the bioinformatics becomes the focus.
