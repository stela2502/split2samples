
unlink(  "Rhapsody_cell_keys_v2_rust.txt" )
length = 6
i = 0

for (str in scan('Rhapsody_cell_keys_v2.txt', what=character())) {
	str = gsub("[\r\n]", "", str)
	i = i+1
	if ( i == 6 ){
		cat ( paste(sep="", str, "\n"), file="Rhapsody_cell_keys_v2_rust.txt", append = TRUE)
		i = 0
	}
	else {
		cat( paste ( sep ="" ,str, " "), file= "Rhapsody_cell_keys_v2_rust.txt", append=TRUE )
	}
}
