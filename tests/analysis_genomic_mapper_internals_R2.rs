// /tests/analysis_gene_mapper_internals.rs

#[cfg(test)]
mod tests {

	use rustody::analysis_genomic_mapper::AnalysisGenomicMapper;
	use rustody::genes_mapper::GenesMapper;
	use rustody::mapping_info::MappingInfo;
	use rustody::genes_mapper::SeqRec;
	use rustody::errors::MappingError;

	use std::path::Path;
	use std::fs;

	use needletail::parse_fastx_file;



	// we need a temp index testData/index/tmp_idx which is ignored by git anyhow

	fn test_this_seqence( seq: &[u8], database:&str, sam_line: Option<&str>, err:Option<MappingError> ){

		let mut results = MappingInfo::new( None, 20.0, 10, None );

		let idx_path = "testData/index/tmp_idx";

		if fs::metadata( idx_path ).is_err() {
	        if let Err(err) = fs::create_dir_all( idx_path ) {
	            eprintln!("Error creating directory {}: {}", idx_path, err);
	        } else {
	            println!("New output directory created successfully!");
	        }
	    }

		let mut idx = GenesMapper::new(0);
    	if Path::new(&database).exists(){

	    	let mut expr_file = parse_fastx_file(database).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    idx.add( &seqrec.seq().to_vec(), id, id, "genetag", 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }

	    }else {
	    	eprintln!("Expression file could not be read - ignoring")
	    }
	    idx.write_index( idx_path );

		// new(_gene_kmers:usize, version:String, specie: String, index:Option<String>, num_threads:usize, exp:&str, _debug: bool  ) -> Self{
		let mut worker = AnalysisGenomicMapper::new( 32, "v1".to_string(), "mouse".to_string(), Some(idx_path.to_string() ), 1, "bd", true);
		let pos = &[0,9, 21,30, 43,52, 52,60 ];

		worker.debug( Some(true) );
		// that contains a cell id for the version of the bd tool

		let r1 = SeqRec::new( b"SomeRead1", b"AGGAGATTAACTGGCCTGCGAGCCTGTTCAGGTAGCGGTGACGACTACATATGCTGCACATTTTTT", b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" );
	    // a read for the seq you wanted analyzed
	    let qual: Vec<u8>  = vec![b'F'; seq.len()];
	    let r2 = SeqRec::new( b"SomeRead2", seq, qual.as_slice() );
	    let data= vec![ (r1, r2 )];

	    let (single_cell_data, sam_strings) = worker.analyze_paralel( &data, &mut results, pos );
	    match err{
	    	Some(e) => {
	    		assert_eq!( single_cell_data.is_empty(), true, "no results in the data object");
	    		assert_eq!( sam_strings.len(), 0, "I go no result for the search" );
	    	},
	    	None=> {
	    		assert_eq!( single_cell_data.is_empty(), false, "there no result in the data object and I expected '{sam_line:?}'");
	    		match sam_line {
	    			Some( sam ) => {
	    				if sam_strings.is_empty(){
			    			panic!("I go no result instead of a sam line!");
			    		}
			    		assert_eq!(sam_strings[0], sam, "We got the expected sam line?");
	    			},
	    			None => {
	    				assert_eq!( sam_strings.len(), 0, "I go no result for the search" );
	    			}
	    		}
	    		
	    	}
	    }
	    
	}

	#[test]
	fn chrM_over_the_edge( ){
		let seq = b"CGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAG";
		let database = "testData/ChrM.fasta.gz";

		let bam_line= "SomeRead20+85\t0\tchrM\t1\t39\t72M1X12M\t*\t0\t0\tCGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.011764706\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );
	}

	#[test]
	fn should_map_not_to_center(){
		let seq = b"CCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCATCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGG";
		let database = "testData/ChrM.fasta.gz";
		let bam_line= "SomeRead2\t0\tchrM\t9731\t40\t90M\t*\t0\t0\tCCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCATCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );
	}

	#[test]
	fn simulate_chrm_read_over_start(){
		let name = "Btla";
		let database = "testData/genes.fasta";
		// start and end of the 'contig'
		//          eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeesssssssssssssssssssssssssss
		let seq = b"ATGTATGTTGTAGCTCCTCAAATAAATTTGTTCCAGCATTAgcactctcacttactaagcATGTTCTA";
		let bam_line= "SomeRead2\t0\t";
		test_this_seqence( seq, database, None, Some(MappingError::NoMatch) );
	}

	#[test]
	fn simulate_chrm_read_over_start2(){
		let name = "Btla";
		let database = "testData/genes.fasta";
		// start and end of the 'contig'
		//          eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeessssssss
		let seq = b"ATGTATGTTGTAGCTCCTCAAATAAATTTGTTCCAGCATTAgcactctc";
		//                               the right gene               the smaller sequence.
		let bam_line= "SomeRead20+41\t0\tBtla\t145\t40\t41M\t*\t0\t0\tATGTATGTTGTAGCTCCTCAAATAAATTTGTTCCAGCATTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );

	}

	#[test]
	fn simulate_obscue_error1(){
		let name = "Btla";
		let database = "testData/ChrM.fasta.gz";
		// start and end of the 'contig'
		//          eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeessssssss
		let seq = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTG";
		//                               the right gene               the smaller sequence.
		let bam_line= "SomeRead2\t0\tchrM\t1\t39\t72M1X17M\t*\t0\t0\tGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.011111111\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );

	}

    #[test]
    fn test_buffer_overflow_issue_in_calculate_cigar(){
        let seq = b"CAAGCAGTTTGCACGTTTGTGATTCTAGAGAGAGAAGACGACGGCGAAGTAGGAGTGG";
        let bam_line= "SomeRead2";
        let database = "testData/Srsf11.fasta.gz";
        test_this_seqence( seq, database, Some(bam_line), None );
    }
    /*#[test]
    fn identify_the_better_database_entry() {
        let name ="Rpl11_int";
        let database = "testData/problematic_match.fasta".to_string();
        let seq = b"GGAGAAAGGCCTGAAGGTGCGGGAGTATGAGTTGCGGAAAAATAACTTCTCGGATACTGGAAACTTTGGTTTTGGAATTCAAGAACACAT";
        let bam_line= "SomeRead20+78\t0\tRpl11\t320\t40\t78M\t*\t0\t0\tGGAGAAAGGCCTGAAGGTGCGGGAGTATGAGTTGCGGAAAAATAACTTCTCGGATACTGGAAACTTTGGTTTTGGAAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
        test_this_seqence( seq, database, Some(bam_line), None );

    }*/

    /*#[test]
    fn idenitfy_sample1() {
    	let database = "testData/problematic_match.fasta.gz".to_string();
    	let seq = b"ATTGTCAAGATGCTACCGTTCAGAGAAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAGAAAA";
    	// sample ids do not create a sam line!
    	test_this_seqence( seq, database, None, None );
    }*/

    // A00681:881:H3MV7DSX7:1:1107:7021:16532 2:N:0:CACAATCCCA+TTGTGGATAT      0       Defb39_int      10595   40      17M18446744073709551274N73M     *       0       0       GAACTAACCAGTACCCCGAGCTCTTGACTCTAGCTGCATATGTATCAAAAGATGGCCTAGTCGGCCATCACTGGAAAGAGAGGCCCATTG      FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFFFFFFF      NH:i:1  HI:i:1  AS:i:40 nM:i:0  RE:A:I  li:i:0  BC:Z:AAGCCTGGAGAC       QT:Z:FFFFFFFFFFFF       CR:Z:GTGTCCTGTCATGACT   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:GTGTCCTGTCATGACT-1 UR:Z:AAGCCTGGAGAC       UZ:Z:FFFFFFFFFFFF       UB:Z:AAGCCTGGAGAC       RG:Z:Sample4:0:1:HN2CKBGX9:1
    // GAACTAACCAGTACCCCGAGCTCTTGACTCTAGCTGCATATGTATCAAAAGATGGCCTAGTCGGCCATCACTGGAAAGAGAGGCCCATTG
    /*#[test]
    fn real_live_splice_site() {
    	let database = "testData/Defb39_human.fasta.gz".to_string();
    	let seq = b"GAACTAACCAGTACCCCGAGCTCTTGACTCTAGCTGCATATGTATCAAAAGATGGCCTAGTCGGCCATCACTGGAAAGAGAGGCCCATTG";
    	// sample ids do not create a sam line!
    	test_this_seqence( seq, database, Some("Text"), None );
    }*/

}
