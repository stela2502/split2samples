//test/needleman_wunsch_affine.rs


#[cfg(test)]
mod tests {
	use rustody::genes_mapper::NeedlemanWunschAffine;
	use rustody::genes_mapper::gene_data::GeneData;
	use rustody::genes_mapper::cigar::Cigar;

	#[test]
	fn test_needleman_wunsch_affine(){
		let mut test = 	NeedlemanWunschAffine::new(10);

		let read = GeneData::new(     b"AAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCATAGGTTTGGTCCTAGC", 
			"read1","gene1", 0  );
		let database = GeneData::new( b"AAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCATAGGTTTGGTCCTAGC", //CTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATC", 
			"database", "database1", 0);

		let nw = test.needleman_wunsch_affine( &read, &database, 1.0 );
		assert_eq!( nw, 0.0, "perfect match has nw of 0");

		let mut cigar= Cigar::new( "" );
		cigar.convert_to_cigar( &test.cigar_vec() );

		assert_eq!( format!("{}",cigar), "50M - None", "get a perfect 50 bp matching result" )


	}

	#[test]
	fn test_needleman_wunsch_affine_gap(){
		let mut test = 	NeedlemanWunschAffine::new(10);

		let read = GeneData::new(     b"AAGCCGGCGTAAAGAGTGAGATCACCCCCATAGGTTTGGTCCTAGC", 
			"read1","gene1", 0  );
		let database = GeneData::new( b"AAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCATAGGTTTGGTCCTAGC", //CTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATC", 
			"database", "database1", 0);

		let nw = test.needleman_wunsch_affine( &read, &database, 1.0 );
		assert!( nw < 0.5, "this NOT perfect match has nw of less than 0.5");

		let mut cigar= Cigar::new( "" );
		cigar.convert_to_cigar( &test.cigar_vec() );

		assert_eq!( format!("{}",cigar), "18M4D28M - None", "get a perfect 50 bp matching result" )


	}


	#[test]
	fn test_needleman_wunsch_affine_insert(){
		let mut test = 	NeedlemanWunschAffine::new(10);

		let read = GeneData::new(     b"AAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCATAGGTTTGGTCCTAGC", 
			"read1","gene1", 0  );
		let database = GeneData::new( b"AAGCCGGCGTAAAGAGTGAGATCACCCCCATAGGTTTGGTCCTAGC", //CTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATC", 
			"database", "database1", 0);

		let nw = test.needleman_wunsch_affine( &read, &database, 1.0 );
		assert!( nw < 0.5, "this NOT perfect match has nw of less than 0.5");

		let mut cigar= Cigar::new( "" );
		cigar.convert_to_cigar( &test.cigar_vec() );

		assert_eq!( format!("{}",cigar), "18M4I28M - None", "get a perfect 50 bp matching result" )


	}


	#[test]
	fn test_needleman_wunsch_affine_large_gap(){
		let mut test = 	NeedlemanWunschAffine::new(10);

		let read = GeneData::new(     b"AAAACGCTTAGCCTAGCATTTCGTGCCAGCCACCGC", 
			"read1","gene1", 0  );
		let database = GeneData::new( b"AAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGC",
			"database", "database1", 0);

		let nw = test.needleman_wunsch_affine( &read, &database, 1.0 );
		assert!( nw < 5.0 ,"mismatch match has nw of less tha  0.5 0");

		let mut cigar= Cigar::new( "" );
		cigar.convert_to_cigar( &test.cigar_vec() );

		assert_eq!( format!("{}",cigar), "17M84D19M - None", "get a perfect 50 bp matching result" )


	}



}
