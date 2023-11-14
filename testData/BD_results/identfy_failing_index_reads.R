fasta = readLines('testData/BD_results/Rustody_S1/mapped_reads.txt')
index = readLines('testData/BD_results/Rustody_S1_index/mapped_reads.txt') 
fasta[1:10]
fasta[1:100]
fastq
fasta
q()
fasta = readLines('testData/BD_results/Rustody_S1/mapped_reads.txt')
index = readLines('testData/BD_results/Rustody_S1_index/mapped_reads.txt') 
index
fasta = fasta[ grep('A00681', fasta) ]
fasta
index = index[ grep('A00681', index) ]
lsngth(fasta)
length(fasta)
length(index)
intersect(fasta, index)
length( (fasta, index) )
`%nin%` <- Negate(`%in%`)
fasta[which(fasta %nin% index)]
length(fasta[which(fasta %nin% index)] )
writeLines( fasta[which(fasta %nin% index)], "fastq_ids_not_in_index.txt", sep="\n" )
system('head fastq_ids_not_in_index.txt')
q()
