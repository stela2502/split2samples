// tests/genes_mapper.rs

#[cfg(test)]
mod tests {

	use rustody::genes_mapper::GenesMapper;
	use rustody::genes_mapper::NeedlemanWunschAffine;
	use needletail::parse_fastx_file;
	use std::path::Path;
	use rustody::errors::MappingError;

	#[test]
	fn test_create(){
		let mut genes = GenesMapper::new(0);

		let ex = "testData/genes.fasta";
		let mut i=0;
		if Path::new(ex).exists(){
	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }

	    assert_eq!( genes.get_gene_count() , i, "Got the expecetd gene counts out" );
	    assert_eq!( genes.get_max_gene_id(), i, "Got the last id out" );
	}

	#[test]
	fn test_create_offset(){
		let mut genes = GenesMapper::new(100);

		let ex = "testData/genes.fasta";
		let mut i=0;
		if Path::new(ex).exists(){
	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }

	    assert_eq!( genes.get_gene_count() , i, "Got the expecetd gene counts out" );

	    assert_eq!( genes.get_max_gene_id(), 100 +i, "Got the last id out" );
	}

	#[test]
	fn test_indexed_genes(){
		let mut genes = GenesMapper::new(100);

		let ex = "testData/genes.fasta";
		#[allow(unused_variables)] //otherwise the compiler falsly brags about this.
		let mut i = 0;
		if Path::new(ex).exists(){
	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }

	    let indexed = genes.as_indexed_genes();

	    let indexed_names = indexed.get_all_gene_names();

	    for i in 0..5{
	    	assert_eq!( indexed.ids_for_gene_names( &vec![indexed_names[i].to_string()])[0], genes.extern_id_for_gname( &indexed_names[i] ).unwrap(), "the id was as expoected" );
	    }

	}

	#[test]
	fn test_wrap_arount_to_max(){
		
		let mut genes = GenesMapper::new(100);
		let mut nwa = NeedlemanWunschAffine::new();
		
		let ex = "testData/genes.fasta";
		#[allow(unused_variables)] //otherwise the compiler falsly brags about this.
		let mut i = 0;
		if Path::new(ex).exists(){
	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }
	    //seq: &[u8], _cellid:u32, nwa: &mut NeedlemanWunschAffine 
	    let result = genes.get_strict( b"CTGCCCCTCTTTTGTGTTGTCTTTTTTTCTTAGACTATCTGTCCTTTCTCCTTGATTTCTAAACTATGTTATTT", 1, &mut nwa );
	    // result is a Result< Vec<MapperResult>, MappingError >
	    match result {
	    	Ok(mapper_result) => {
	    		// the SAM entry for this match was 
				//..    0       Ube2c   18446744073709551615    38      1I20M1D53M      *       0       0       CTGCCCCTCTTTTGTGTTGTCTTTTTTTCTTAGACTATCTGTCCTTTCTCCTTGATTTCTAAACTATGTTATTT      AAAAAEEEEA</EEEE/<EAAEE//EEE/EEEEEEEAAEAAE<E//EE/EE/<EAAEEEE<//EA//EE<A<EE      NH:i:1  HI:i:1  AS:i:38 nM:i:0.026666667        RE:A:I  li:i:0  BC:Z:AAAAACTC   QT:Z:/EEEEEEE   CR:Z:AACTTCTCCACCTGAGTCAAGGGTCAG        CY:Z:AAAAAEEEE6EEE/EEEEEEEEEEEEE        CB:Z:AACTTCTCCACCTGAGTCAAGGGTCAG-1      UR:Z:AAAAACTC   UZ:Z:/EEEEEEE   UB:Z:AAAAACTC   RG:Z:Sample4:0:1:HN2CKBGX9:1
				// which leads to the assumption that we should actually clip the starting I from the match - that could fix everything!
				//panic!("I got a result {mapper_result:?}")
				//I got a result [MapperResult { gene_id: 440, start: 18446744073709551615, save: true, cigar: Some(Cigar { cigar: "1I20M1D52M", fixed: Some(Na) }), mapq: 38, score: 9, nw: 0.12328767, edit_dist: 0.027027028, gene_name: "Ube2c", db_length: 330 }]
				// cool reproduced the error!
				assert_eq!( &format!("{}",mapper_result[0].get_cigar()), "1I20M1D52M - Some(Na)", "Starting 1I still in the match? {}", mapper_result[0].start() );

				assert_eq!( mapper_result[0].start(), 0, "start integer overflow? {}", mapper_result[0].start() );

	    	},
	    	Err(e) => {
	    		panic!("Insetad of a result I got the error {e:?}")
	    	}
	    }
		
	}

	//TCAAGTGTCTGATGTTCCCTGTGAGCTTGGGTTCAGTGTGAAGAACTGTGGAGCCCAGCCTGCCCTGCACACC
	#[test]
	fn test_bad_cigar_length(){
		
		let mut genes = GenesMapper::new(100);
		let mut nwa = NeedlemanWunschAffine::new();
		//nwa.set_debug(true);
		//genes.debug(Some(true));
	
		let ex = "testData/genes.fasta";
		#[allow(unused_variables)] //otherwise the compiler falsly brags about this.
		let mut i = 0;
		if Path::new(ex).exists(){
	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }
	    //seq: &[u8], _cellid:u32, nwa: &mut NeedlemanWunschAffine 

	    //
	    //
	    //TCAAGTGTCTGATGTTCCCTGTGAGCTTGGGTTCAGTGTGAAGAACTGTGGAGCCCAGCCTGCCCTGCACACC
	    // CTGCCCCTCTTTTGTGTTGTCTTTTTTCTTAGACTATCTGTCCTTTCTCCTTGATTTCTA

	    let result = genes.get_strict( b"TCAAGTGTCTGATGTTCCCTGTGAGCTTGGGTTCAGTGTGAAGAACTGTGGAGCCCAGCCTGCCCTGCACACC", 1, &mut nwa );
	    // result is a Result< Vec<MapperResult>, MappingError >
	    match result {
	    	Ok(mapper_result) => {
	    		// the SAM entry for this match was 
				//..    0       Ube2c   18446744073709551615    38      1I20M1D53M      *       0       0       CTGCCCCTCTTTTGTGTTGTCTTTTTTTCTTAGACTATCTGTCCTTTCTCCTTGATTTCTAAACTATGTTATTT      AAAAAEEEEA</EEEE/<EAAEE//EEE/EEEEEEEAAEAAE<E//EE/EE/<EAAEEEE<//EA//EE<A<EE      NH:i:1  HI:i:1  AS:i:38 nM:i:0.026666667        RE:A:I  li:i:0  BC:Z:AAAAACTC   QT:Z:/EEEEEEE   CR:Z:AACTTCTCCACCTGAGTCAAGGGTCAG        CY:Z:AAAAAEEEE6EEE/EEEEEEEEEEEEE        CB:Z:AACTTCTCCACCTGAGTCAAGGGTCAG-1      UR:Z:AAAAACTC   UZ:Z:/EEEEEEE   UB:Z:AAAAACTC   RG:Z:Sample4:0:1:HN2CKBGX9:1
				// which leads to the assumption that we should actually clip the starting I from the match - that could fix everything!
				//panic!("I got a result {mapper_result:?}")
				//I got a result [MapperResult { gene_id: 440, start: 18446744073709551615, save: true, cigar: Some(Cigar { cigar: "1I20M1D52M", fixed: Some(Na) }), mapq: 38, score: 9, nw: 0.12328767, edit_dist: 0.027027028, gene_name: "Ube2c", db_length: 330 }]
				// cool reproduced the error!
				println!("I got this MapperResult Vec: {mapper_result:?}");

				assert_eq!( &format!("{}",mapper_result[0].get_cigar()), "26M1X2I3M2X3M1X22M15S - Some(StartInsert)", "complicated cigar {}", mapper_result[0].get_cigar() );

				assert_eq!( mapper_result[0].get_cigar().calculate_covered_nucleotides( &format!("{}",mapper_result[0].get_cigar()) ),( 0,0 ), "start - integer overflow? {}", mapper_result[0].start() );

	    	},
	    	Err(e) => {
	    		assert_eq!( e, MappingError::NoMatch, "Actually this makes more sense for this read" );
	    	}
	    }
		
	}




	//NB501227:252:H2GVTBGXF:1:11101:21151:1444 2:N:0:GCTACGCT+AGGCTATA0+64   0       H2-K1   18446744073709551605    38      1M1X61M1X       *       0       0       CCTCTGCCCTGTGAAGTGTCTGATGTTCCCTGTGAGCCTATGGACTCAATGTGAAGAACTGTGG        AAAAA/E<E<<E6AE/EAE/A//6/E<AAEAEEAA/EAEEE/EEAE/A<A/E/EE6/AAE/A//        NH:i:1  HI:i:1  AS:i:38 nM:i:0.03125    RE:A:I  li:i:0  BC:Z:AAGGACGA   QT:Z:EEAE6EEE   CR:Z:TAAGTTCGAAGGATTCAAATGGGACTC        CY:Z:6AAAAEEAA/A/EEEEAEEEAEEEEE6        CB:Z:TAAGTTCGAAGGATTCAAATGGGACTC-1      UR:Z:AAGGACGA   UZ:Z:EEAE6EEE   UB:Z:AAGGACGA   RG:Z:Sample4:0:1:HN2CKBGX9:1
	#[test]
	fn test_wrap_arount_to_max2(){
		
		let mut genes = GenesMapper::new(100);
		let mut nwa = NeedlemanWunschAffine::new();
		
		nwa.set_debug(true);

		let ex = "testData/genes.fasta";
		#[allow(unused_variables)] //otherwise the compiler falsly brags about this.
		let mut i = 0;
		if Path::new(ex).exists(){
	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }
	    genes.debug(Some(true));
	    //seq: &[u8], _cellid:u32, nwa: &mut NeedlemanWunschAffine 
	    let result = genes.get_strict( b"CCTCTGCCCTGTGAAGTGTCTGATGTTCCCTGTGAGCCTATGGACTCAATGTGAAGAACTGTGG", 1, &mut nwa );
	    // result is a Result< Vec<MapperResult>, MappingError >
	    match result {
	    	Ok(mapper_result) => {
	    		// the SAM entry for this match was 
				assert_eq!( &format!("{}",mapper_result[0].get_cigar()), "1M1X51M - Some(Na)", "Starting 1I still in the match? {}", mapper_result[0].start() );

				assert_eq!( mapper_result[0].start(), 0, "start integer overflow? {}", mapper_result[0].start() );

	    	},
	    	Err(e) => {
	    		panic!("Insetad of a result I got the error {e:?}")
	    	}
	    }
		
	}

}
