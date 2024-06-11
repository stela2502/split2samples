// /tests/analysis_gene_mapper_internals.rs

#[cfg(test)]
mod tests {

	use rustody::analysis::AnalysisGeneMapper;
	use rustody::mapping_info::MappingInfo;
	use rustody::genes_mapper::SeqRec;
	use rustody::errors::MappingError;

	fn test_this_seqence( seq: &[u8], database:String, sam_line: Option<&str>, err:Option<MappingError> ){

		let mut results = MappingInfo::new( None, 20.0, 10, None );
		let mut worker = AnalysisGeneMapper::new( 32, "v1".to_string(), Some(database),
    		None, "mouse".to_string(), None, 1, "bd", true);
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
		let database = "testData/ChrM.fasta.gz".to_string();

		let bam_line= "SomeRead20+85\t0\tchrM\t1\t39\t72M1X12M\t*\t0\t0\tCGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.011764706\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );
	}

	#[test]
	fn should_map_not_to_center(){
		let seq = b"CCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCATCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGG";
		let database = "testData/ChrM.fasta.gz".to_string();
		let bam_line= "SomeRead2\t0\tchrM\t9731\t40\t90M\t*\t0\t0\tCCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCATCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );
	}

	#[test]
	fn simulate_chrm_read_over_start(){
		//let name = "Btla";
		let database = "testData/genes.fasta".to_string();
		// start and end of the 'contig'
		//          eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeesssssssssssssssssssssssssss
		let seq = b"ATGTATGTTGTAGCTCCTCAAATAAATTTGTTCCAGCATTAgcactctcacttactaagcATGTTCTA";
		let bam_line= "SomeRead2\t0\t";
		test_this_seqence( seq, database, None, Some(MappingError::NoMatch) );
	}

	#[test]
	fn simulate_chrm_read_over_start2(){
		//let name = "Btla";
		let database = "testData/genes.fasta".to_string();
		// start and end of the 'contig'
		//          eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeessssssss
		let seq = b"ATGTATGTTGTAGCTCCTCAAATAAATTTGTTCCAGCATTAgcactctc";
		//                               the right gene               the smaller sequence.
		let bam_line= "SomeRead20+41\t0\tBtla\t145\t40\t41M\t*\t0\t0\tATGTATGTTGTAGCTCCTCAAATAAATTTGTTCCAGCATTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );

	}

	#[test]
	fn simulate_obscue_error1(){
		//let name = "Btla";
		let database = "testData/ChrM.fasta.gz".to_string();
		// start and end of the 'contig'
		//          eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeessssssss
		let seq = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTG";
		//                               the right gene               the smaller sequence.
		let bam_line= "SomeRead2\t0\tchrM\t1\t39\t72M1X17M\t*\t0\t0\tGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.011111111\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );

	}

    #[test]
    fn identify_the_better_database_entry() {
        //let name ="Rpl11_int";
        let database = "testData/problematic_match.fasta.gz".to_string();
        let seq = b"GGAGAAAGGCCTGAAGGTGCGGGAGTATGAGTTGCGGAAAAATAACTTCTCGGATACTGGAAACTTTGGTTTTGGAATTCAAGAACACAT";
        let bam_line= "SomeRead20+78\t0\tRpl11\t320\t40\t78M\t*\t0\t0\tGGAGAAAGGCCTGAAGGTGCGGGAGTATGAGTTGCGGAAAAATAACTTCTCGGATACTGGAAACTTTGGTTTTGGAAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
        test_this_seqence( seq, database, Some(bam_line), None );
        //test_this_seqence( seq, database, None, Some( MappingError::NoMatch ) );
    }

    #[test]
    fn idenitfy_sample1() {
    	let database = "testData/problematic_match.fasta.gz".to_string();
    	let seq = b"ATTGTCAAGATGCTACCGTTCAGAGAAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAGAAAA";
    	// sample ids do not create a sam line!
    	test_this_seqence( seq, database, None, None );
    }

    // A00681:881:H3MV7DSX7:1:1107:7021:16532 2:N:0:CACAATCCCA+TTGTGGATAT      0       Defb39_int      10595   40      17M18446744073709551274N73M     *       0       0       GAACTAACCAGTACCCCGAGCTCTTGACTCTAGCTGCATATGTATCAAAAGATGGCCTAGTCGGCCATCACTGGAAAGAGAGGCCCATTG      FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFFFFFFF      NH:i:1  HI:i:1  AS:i:40 nM:i:0  RE:A:I  li:i:0  BC:Z:AAGCCTGGAGAC       QT:Z:FFFFFFFFFFFF       CR:Z:GTGTCCTGTCATGACT   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:GTGTCCTGTCATGACT-1 UR:Z:AAGCCTGGAGAC       UZ:Z:FFFFFFFFFFFF       UB:Z:AAGCCTGGAGAC       RG:Z:Sample4:0:1:HN2CKBGX9:1
    // GAACTAACCAGTACCCCGAGCTCTTGACTCTAGCTGCATATGTATCAAAAGATGGCCTAGTCGGCCATCACTGGAAAGAGAGGCCCATTG
    #[test]
    fn real_live_splice_site() {
    	let database = "testData/Defb39_human.fasta.gz".to_string();
    	let seq = b"GAACTAACCAGTACCCCGAGCTCTTGACTCTAGCTGCATATGTATCAAAAGATGGCCTAGTCGGCCATCACTGGAAAGAGAGGCCCATTG";
    	// sample ids do not create a sam line!
    	test_this_seqence( seq, database, None, Some( MappingError::NoMatch ) );
    }

    //insertion/deletion border cases
    #[test]
    fn insert_start() {
    	// start Insert ENSMUST00000005017 Hdgf
    	let seq = b"TTGGTCTCTGGTGTTTTCTCACATCCAGTTTGTAGCCCACTAAGAATACAGGGAAAGTGTCCTCCAGCCTCTTCTAGTGGTTTCTTACA";
    	let database = "testData/Hdgf.fasta.gz".to_string();
		let bam_line= "SomeRead2\t0\tENSMUST00000005017\t269\t39\t53M1D36M\t*\t0\t0\tTTGGTCTCTGGTGTTTTCTCACATCCAGTTTGTAGCCCACTAAGAATACAGGGAAAGTGTCCTCCAGCCTCTTCTAGTGGTTTCTTACA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.022222223\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );
    }
    #[test]
    fn insert_end() {
    	// end Insert ENSMUST00000026565 Ifitm3
    	let seq = b"ACCTTGGTCCTCAGCATCCTGATGGTTGTTATCACCATTGTTAGTGTCATCATCATTGTTCTTAACGCTCAAAAACCTTCACACTTAAT";
    	let database = "testData/Ifitm3.fasta.gz".to_string();
    	let bam_line= "SomeRead2\t0\tENSMUST00000026565\t275\t39\t70M1I18M\t*\t0\t0\tACCTTGGTCCTCAGCATCCTGATGGTTGTTATCACCATTGTTAGTGTCATCATCATTGTTCTTAACGCTCAAAAACCTTCACACTTAAT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.022222223\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
    	test_this_seqence( seq, database, Some(bam_line), None );
    }
    #[test]
    fn deletion_start() {
    	// start deletion ENSMUST00000109641 Sec61g
		let seq = b"ATTATTGTGGTGGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTG";
		let database = "testData/Sec61g.fasta.gz".to_string();
		let bam_line= "SomeRead2\t0\tENSMUST00000109641\t286\t39\t8M1D81M\t*\t0\t0\tATTATTGTGGTGGCTGAGTCCTTCTCATCATGGGACGAGTGAGCCAGAGCGGGGGAAAGGGCATGAAGTAAAGCGTTGCCTGAATGCTG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.022222223\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );
    }
	#[test]
    fn deletion_end() {
		// end deletion ENSMUST00000045897 Ptma
		let seq = b"GGTGTGACCATGTTCATTATAATCTCAAAGGAGAAAAAAAAACCTTGTAAAAAAAAGCAAAAACAACAACAAAAAAACAATCTTATTCC";
		let database = "testData/Ptma.fasta.gz".to_string();
		let bam_line= "SomeRead2\t0\tENSMUST00000045897\t184\t39\t33M1D56M\t*\t0\t0\tGGTGTGACCATGTTCATTATAATCTCAAAGGAGAAAAAAAAACCTTGTAAAAAAAAGCAAAAACAACAACAAAAAAACAATCTTATTCC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.022222223\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, Some(bam_line), None );
	}

	#[test]
	fn what_is_the_problem_here() {
		let seq = b"TCGTGTCTCTAGCTGCATATGTAGCAGAGAATGGCCTAGTCGGCCATCACTGGGAAGAGAGGCCCCTTGGTCTTGCAAACTTTATATGC";
		let database = "testData/Hectd2os.fasta.gz".to_string();
		let bam_line= "Fix Me";
		test_this_seqence( seq, database, Some(bam_line), None );
	}

}
