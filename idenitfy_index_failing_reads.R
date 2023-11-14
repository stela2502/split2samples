fasta = readLines('testData/BD_results/Rustody_S1/mapped_reads.txt')
index = readLines('testData/BD_results/Rustody_S1_index/mapped_reads.txt')

fasta = fasta[ grep('A00681', fasta) ]
index = index[ grep('A00681', index) ]

`%nin%` <- Negate(`%in%`)
interesting = fasta[which(fasta %nin% index)]
length( interesting )

R1 = unlist( lapply( interesting[grep ("^R1", interesting )], function(x){ sub(".{3}","", x )}  ))
R2 = unlist( lapply( interesting[grep ("^R2", interesting )], function(x){ sub(".{3}","", x )}  ))

writeLines( R1, "R1_ids_not_in_index.txt", sep="\n" )
writeLines( R2, "R2_ids_not_in_index.txt", sep="\n" )