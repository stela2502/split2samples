library(Seurat)

system( paste( "target/release/quantify_rhapsody ",
	"-r  testData/1e5_mRNA_S1_R1_001.fastq.gz ",
	"-f testData/1e5_mRNA_S1_R2_001.fastq.gz ",
	"-o testData/output_1e5 ",
	"-s mouse  ",
	"-e testData/genes.fasta ",
	"-a testData/MyAbSeqPanel.fasta ",
	"-m 10 -v v1"  ) )


data = Read10X( "testData/output_1e5/BD_Rhapsody_expression")

good = readRDS(file = "Rtest/validatesRes.rds")

if ( all.equal( data, good ) ){
	print ( "Test OK" )
}else {
	print ( "Test FAIL")
}
