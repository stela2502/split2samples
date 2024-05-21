// /tests/analysis_gene_mapper_internals.rs

#[cfg(test)]
mod tests {

	use rustody::analysis_genemapper::AnalysisGeneMapper;
	use rustody::mapping_info::MappingInfo;
	use rustody::genes_mapper::SeqRec;
	use rustody::errors::MappingError;

	fn test_this_seqence( seq: &[u8], database:String, sam_line: &str, err:Option<MappingError> ){

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
	    		assert_eq!( sam_strings.len(), 0, "I go no result for the search" );
	    	},
	    	None=> {
	    		assert_eq!(sam_strings[0], sam_line, "We got the expected sam line?");
	    	}
	    }
	    
	}

	#[test]
	fn chrM_over_the_edge( ){
		let seq = b"CGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAG";
		let database = "testData/ChrM.fasta.gz".to_string();

		let bam_line= "SomeRead20+85\t0\tchrM\t1\t39\t72M1X12M\t*\t0\t0\tCGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.011764706\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, bam_line, None );
	}

	#[test]
	fn should_map_not_to_center(){
		let seq = b"CCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCATCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGG";
		let database = "testData/ChrM.fasta.gz".to_string();
		let bam_line= "SomeRead2\t0\tchrM\t9731\t40\t90M\t*\t0\t0\tCCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCATCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, bam_line, None );
	}

	#[test]
	fn simulate_chrM_read_over_start(){
		let name = "Btla";
		let database = "testData/genes.fasta".to_string();
		// start and end of the 'contig'
		//          eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeesssssssssssssssssssssssssss
		let seq = b"ATGTATGTTGTAGCTCCTCAAATAAATTTGTTCCAGCATTAgcactctcacttactaagcATGTTCTA";
		let bam_line= "SomeRead2\t0\t";
		test_this_seqence( seq, database, bam_line, Some(MappingError::NoMatch) );
	}

	#[test]
	fn simulate_chrM_read_over_start2(){
		let name = "Btla";
		let database = "testData/genes.fasta".to_string();
		// start and end of the 'contig'
		//          eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeessssssss
		let seq = b"ATGTATGTTGTAGCTCCTCAAATAAATTTGTTCCAGCATTAgcactctc";
		//                               the right gene               the smaller sequence.
		let bam_line= "SomeRead20+41\t0\tBtla\t145\t40\t41M\t*\t0\t0\tATGTATGTTGTAGCTCCTCAAATAAATTTGTTCCAGCATTA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, bam_line, None );

	}

	#[test]
	fn simulate_obscue_error1(){
		let name = "Btla";
		let database = "testData/ChrM.fasta.gz".to_string();
		// start and end of the 'contig'
		//          eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeessssssss
		let seq = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTG";
		//                               the right gene               the smaller sequence.
		let bam_line= "SomeRead2\t0\tchrM\t1\t39\t72M1X17M\t*\t0\t0\tGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:39\tnM:i:0.011111111\tRE:A:I\tli:i:0\tBC:Z:GCTGCACA\tQT:Z:FFFFFFFF\tCR:Z:AGGAGATTAGCCTGTTCAACTACATAT\tCY:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AGGAGATTAGCCTGTTCAACTACATAT-1\tUR:Z:GCTGCACA\tUZ:Z:FFFFFFFF\tUB:Z:GCTGCACA\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, database, bam_line, None );

	}

}