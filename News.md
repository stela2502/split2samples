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
