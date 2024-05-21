// /tests/analysis_gene_mapper_internals.rs

#[cfg(test)]
mod tests {

	use rustody::analysis_genemapper::AnalysisGeneMapper;
	use rustody::mapping_info::MappingInfo;
	use rustody::genes_mapper::SeqRec;
	use rustody::errors::MappingError;

	fn test_this_seqence( seq: &[u8], exp: &str, version:String, specie:String,  sam_line: &str, err:Option<MappingError> ){

		let mut results = MappingInfo::new( None, 20.0, 10, None );
		let mut worker = AnalysisGeneMapper::new( 32, version, Some("testData/ChrM.fasta.gz".to_string()),
    		None, specie, None, 1, exp, true);

		// the pos is not used any more in the function - so any value here is OK
		let pos = &[0,9, 21,30, 43,52, 52,60 ];

		worker.debug( Some(true) );
		let qual: Vec<u8>  = vec![b'F'; seq.len()];
		// that contains a cell id for the version of the bd tool
		let r1 = SeqRec::new( b"SomeRead1", seq, qual.as_slice() );
	    // a read for the seq you wanted analyzed
	    let qual: Vec<u8>  = vec![b'F'; seq.len()];

	    let r2 = SeqRec::new( b"SomeRead2", 
	    	b"CATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGA", 
	    	&[b'F'; 79]
	    	);
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
		let seq = b"AGCTTTAAGGAGCATACTTTACAAGGAT";
		let specie = "mouse".to_string();
		let exp="10x";
		let version = "Single Cell Multiome (ATAC+GEX) v1".to_string();

		let bam_line= "SomeRead2\t0\tchrM\t12115\t40\t79M\t*\t0\t0\tCATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:CTTTACAAGGAT\tQT:Z:FFFFFFFFFFFF\tCR:Z:AGCTTTAAGGAGCATA\tCY:Z:FFFFFFFFFFFFFFFF\tCB:Z:AGCTTTAAGGAGCATA-1\tUR:Z:CTTTACAAGGAT\tUZ:Z:FFFFFFFFFFFF\tUB:Z:CTTTACAAGGAT\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, exp, version, specie, bam_line, None );
	}

	#[test]
	fn chrM_cellid_problem( ){
		let seq = b"GCTTGTTGTTACTTGCCGGAAATCGTAC";
		let specie = "mouse".to_string();
		let exp="10x";
		let version = "Single Cell Multiome (ATAC+GEX) v1".to_string();

		let bam_line= "SomeRead2\t0\tchrM\t12115\t40\t79M\t*\t0\t0\tCATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tNH:i:1\tHI:i:1\tAS:i:40\tnM:i:0\tRE:A:I\tli:i:0\tBC:Z:CGGAAATCGTAC\tQT:Z:FFFFFFFFFFFF\tCR:Z:GCTTGTTGTTACTTGC\tCY:Z:FFFFFFFFFFFFFFFF\tCB:Z:GCTTGTTGTTACTTGC-1\tUR:Z:CGGAAATCGTAC\tUZ:Z:FFFFFFFFFFFF\tUB:Z:CGGAAATCGTAC\tRG:Z:Sample4:0:1:HN2CKBGX9:1";
		test_this_seqence( seq, exp, version, specie, bam_line, None );
	}

}