// tests/genes_mapper.rs

#[cfg(test)]
mod tests {

	use rustody::genes_mapper::{GenesMapper, NeedlemanWunschAffine, CigarEndFix, Cigar};
	use needletail::parse_fastx_file;
	use std::path::Path;
	use std::fs;
	use rustody::errors::MappingError;


	fn test_this_seqence( seq: &[u8], database:String, cigar:Option<Cigar>, gene: Option<String>, start:Option<usize>, err: Option<MappingError> ){
    	let mut genes = GenesMapper::new(100);
		let mut nwa = NeedlemanWunschAffine::new();
		
		nwa.set_debug(true);

		//let database = "testData/genes.fasta";
		#[allow(unused_variables)] //otherwise the compiler falsly brags about this.
		let mut i = 0;
		if Path::new(&database).exists(){
	    	let mut expr_file = parse_fastx_file(database).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), &id, &id, &id, 0 );
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
	    let result = genes.get_strict( seq, 1, &mut nwa );
	    // result is a Result< Vec<MapperResult>, MappingError >
	    
	    let (matches, err_matches) = match result {
	    	Ok(mapper_result) => {
	    		// the SAM entry for this match was
	    		let mut matches = 0;
	    		matches += match gene {
	    			Some(gene_name) => {
	    				assert_eq!( mapper_result[0].get_name(), gene_name, "not the expected gene name {gene_name}?");
	    				1
	    			},
	    			None => {0}
	    		};
	    		matches += match start{
	    			Some(pos) => {
	    				assert_eq!( mapper_result[0].start(), pos, "not the expected start position on the database {pos}?");
	    				1
	    			},
	    			None => {0}
	    		};

	    		matches += match cigar {
	    			Some(cig) => {
	    				assert_eq!( &format!("{}",mapper_result[0].get_cigar()), &format!("{}",cig), "The Cigar is as expected?" );
	    				1
	    			},
	    			None => {0}
	    		};
	    		(matches, 0)
	    	},
	    	Err(e) => {
	    		let err_matches =  match err{
	    			Some(expected_error) => {
	    				assert_eq!( e, expected_error, "Not the expected error {expected_error:?}" );
	    				1
	    			},
	    			None => {
	    				panic!("And error was thrown even so it was not expected: {e:?}");
	    			}
	    		};
	    		(0, err_matches)
	    	}
	    };
	    match err{
			Some(expected_error) => {
				assert_eq!( err_matches, 1, "The expected error was not thrown {expected_error:?}!");
			},
			None => {
				assert!(matches > 0, "This function needs to run at least one test ({matches}) in order to make sense!")
			}
		};

    }

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
		                    genes.add( &seqrec.seq().to_vec(), &id, &id, &id, 0 );
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
		                    genes.add( &seqrec.seq().to_vec(), &id, &id, &id, 0 );
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
		                    genes.add( &seqrec.seq().to_vec(), &id, &id, &id, 0 );
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

		let mut cig = Cigar::new( "21M1I52M1D" );
		cig.fixed = Some(CigarEndFix::Na);
		test_this_seqence( 
			b"CTGCCCCTCTTTTGTGTTGTCTTTTTTTCTTAGACTATCTGTCCTTTCTCCTTGATTTCTAAACTATGTTATTT",
			"testData/genes.fasta".to_string(),
			Some(cig),
			Some( "Ube2c".to_string() ),
			Some(0),
			None
		);
	
	}

	////..    0       Ube2c   18446744073709551615    38      1I20M1D53M      *       0       0       CTGCCCCTCTTTTGTGTTGTCTTTTTTTCTTAGACTATCTGTCCTTTCTCCTTGATTTCTAAACTATGTTATTT      AAAAAEEEEA</EEEE/<EAAEE//EEE/EEEEEEEAAEAAE<E//EE/EE/<EAAEEEE<//EA//EE<A<EE      NH:i:1  HI:i:1  AS:i:38 nM:i:0.026666667        RE:A:I  li:i:0  BC:Z:AAAAACTC   QT:Z:/EEEEEEE   CR:Z:AACTTCTCCACCTGAGTCAAGGGTCAG        CY:Z:AAAAAEEEE6EEE/EEEEEEEEEEEEE        CB:Z:AACTTCTCCACCTGAGTCAAGGGTCAG-1      UR:Z:AAAAACTC   UZ:Z:/EEEEEEE   UB:Z:AAAAACTC   RG:Z:Sample4:0:1:HN2CKBGX9:1
	#[test]
	fn test_bad_cigar_length(){
		
		//let mut cig = Cigar::new( "26M1X2I3M2X3M1X22M15S" );
		//cig.fixed = Some(CigarEndFix::StartInsert);
		test_this_seqence( 
			b"TCAAGTGTCTGATGTTCCCTGTGAGCTTGGGTTCAGTGTGAAGAACTGTGGAGCCCAGCCTGCCCTGCACACC",
			"testData/genes.fasta".to_string(),
			None,
			None,
			None,
			Some(MappingError::NoMatch),
		);
	}




	//NB501227:252:H2GVTBGXF:1:11101:21151:1444 2:N:0:GCTACGCT+AGGCTATA0+64   0       H2-K1   18446744073709551605    38      1M1X61M1X       *       0       0       CCTCTGCCCTGTGAAGTGTCTGATGTTCCCTGTGAGCCTATGGACTCAATGTGAAGAACTGTGG        AAAAA/E<E<<E6AE/EAE/A//6/E<AAEAEEAA/EAEEE/EEAE/A<A/E/EE6/AAE/A//        NH:i:1  HI:i:1  AS:i:38 nM:i:0.03125    RE:A:I  li:i:0  BC:Z:AAGGACGA   QT:Z:EEAE6EEE   CR:Z:TAAGTTCGAAGGATTCAAATGGGACTC        CY:Z:6AAAAEEAA/A/EEEEAEEEAEEEEE6        CB:Z:TAAGTTCGAAGGATTCAAATGGGACTC-1      UR:Z:AAGGACGA   UZ:Z:EEAE6EEE   UB:Z:AAGGACGA   RG:Z:Sample4:0:1:HN2CKBGX9:1
	#[test]
	fn test_wrap_arount_to_max2(){

		let mut cig = Cigar::new( "1M1X51M" );
		cig.fixed = Some(CigarEndFix::Na);
		test_this_seqence( 
			b"CCTCTGCCCTGTGAAGTGTCTGATGTTCCCTGTGAGCCTATGGACTCAATGTGAAGAACTGTGG",
			"testData/genes.fasta".to_string(),
			Some(cig),
			Some( "H2-K1".to_string() ),
			Some(0),
			None
		);
		
	}

	//A00519:856:HWHNJDSXY:2:1131:16712:3662 2:N:0:ATGACGTCGC+ATCCTGACCT      0       chrM    1       39      72M1X12M        *       0       0       CGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAG      FFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFFFFFFFFFFFFFFFFFF      NH:i:1  HI:i:1  AS:i:39 nM:i:0.011764706        RE:A:I  li:i:0  BC:Z:AACTGTCTGCGT       QT:Z:FFFFFFFFFFFF       CR:Z:TACTTGTTCTTACTCG   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:TACTTGTTCTTACTCG-1 UR:Z:AACTGTCTGCGT       UZ:Z:FFFFFFFFFFFF       UB:Z:AACTGTCTGCGT       RG:Z:Sample4:0:1:HN2CKBGX9:1
    //The read sequence is too long 
    
    // CGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAG

    #[test]
	fn test_cigar_sequence_length_mismatch(){

		let mut cig = Cigar::new( "72M1X12M" );
		cig.fixed = Some(CigarEndFix::Na);
		// cigar:Option<Cigar>, gene: Option<String>, start:Option<usize>, err
		test_this_seqence( 
			b"CGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAG",
			"testData/ChrM.fasta.gz".to_string(),
			Some(cig),
			Some( "chrM".to_string() ),
			Some(0),
			None
		);
	}

	#[test]
	fn create_save_and_load_an_index(){

		let mut genes = GenesMapper::new(100);
		let database = "testData/genes.fasta";

		#[allow(unused_variables)] //otherwise the compiler falsly brags about this.
		let mut i = 0;
		if Path::new(&database).exists(){
	    	let mut expr_file = parse_fastx_file(database).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), &id, &id, &id,0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }

		assert_eq!( genes.len(), 466, "the right amount of genes" );
		assert_eq!( genes.depth(), 41716, "the expected amount of mappers");

		let outpath = "testData/output_index_test/genes";

		if fs::metadata(&outpath).is_err() {
        if let Err(err) = fs::create_dir_all(&outpath) {
	            eprintln!("Error creating directory {}: {}", &outpath, err);
	        } else {
	            println!("New output directory created successfully!");
	        }
	    }

		match genes.write_index( "testData/output_index_test/genes"  ){
			Ok(_) => { },
			Err(e) => panic!("Writing of the index failed with the error {e}")
		};

		assert_eq!( Path::new("testData/output_index_test/genes/index.bin").exists(), true, "the index file exists");
		assert_eq!( Path::new("testData/output_index_test/genes/indexed_sequences.fa.gz").exists(), true, "the fasta file exists");

		let mut from_binary = match GenesMapper::load_index( "testData/output_index_test/genes" ){
			Ok(mapper) => mapper,
			Err(e) => {panic!("The loading of the index failed with the error {e}")}
		};

		assert_eq!( genes.len(), from_binary.len(), "from binary the right amount of genes" );
		assert_eq!( genes.depth(), from_binary.depth(), "from binary the expected amount of mappers");


	}
	

}
