target/release/quantify_rhapsody -r testData/cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f testData/cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o testData/BD_results/Rustody_S1 -s mouse  -e testData/2276_20220531_chang_to_rpl36a_amplicons.fasta -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96 > testData/BD_results/Rustody_S1/mapped_reads.txt
target/release/quantify_rhapsody -r testData/cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f testData/cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o testData/BD_results/Rustody_S1 -s mouse  -e testData/addOn.fa -i testData/mapperTest/index -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96 > testData/BD_results/Rustody_S1_index/mapped_reads.txt

Rscript idenitfy_index_failing_reads.R

zgrep -f R1_ids_not_in_index.txt testData/cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -A3 | sed '/--/d' | gzip > testData/failing_in_index_cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz
zgrep -f R2_ids_not_in_index.txt testData/cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -A3 | sed '/--/d' | gzip > testData/failing_in_index_cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz


target/release/quantify_rhapsody_multi -r testData/failing_in_index_cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f testData/failing_in_index_cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o testData/BD_results/Rustody_S1_failed_in_indexing -s mouse  -e testData/2276_20220531_chang_to_rpl36a_amplicons.fasta -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96 
target/release/quantify_rhapsody_multi -r testData/failing_in_index_cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f testData/failing_in_index_cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o testData/BD_results/Rustody_S1_index_failed -s mouse  -e testData/addOn.fa -i testData/mapperTest/index -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96

# there are still indices that can not be mapped using the muti program?!

target/release/quantify_rhapsody -r testData/failing_in_index_cells.1.Rhapsody_SV_index1_S1_R1_001.fastq.gz -f testData/failing_in_index_cells.1.Rhapsody_SV_index1_S1_R2_001.fastq.gz -o testData/BD_results/Rustody_S1_failed_in_indexing -s mouse  -e testData/2276_20220531_chang_to_rpl36a_amplicons.fasta -a testData/MyAbSeqPanel.fasta -m 200 -v v2.96  > testData/BD_results/Rustody_S1_failed_in_indexing/look_into_these_indexes.txt


# damn these are all mapping using the not multi program?! How so??