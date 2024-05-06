// test/creating_cigars.rs



#[cfg(test)]
mod tests {
	use rustody::genes_mapper::Cigar;
	use rustody::genes_mapper::gene_data::GeneData;
	//use rustody::traits::BinaryMatcher;
	use rustody::genes_mapper::NeedlemanWunschAffine;

	#[test]
	fn test_full_match(){
		let seq1 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT";
		let seq2 = &seq1.clone();
		let gd1 = GeneData::new( seq1, "read", "chrM", 0 );
		let gd2 = GeneData::new( seq2, "database", "chrM", 0 );
		let mut cigar = Cigar::new("");
		let mut nwa = NeedlemanWunschAffine::new(91);

		let _nw = nwa.needleman_wunsch_affine( &gd1, &gd2, 0.4);
		cigar.convert_to_cigar( &nwa.cigar_vec() );
		cigar.clean_up_cigar(&gd1, &gd2);

		assert_eq!( &cigar.to_string(), "90M", "I expected 90M as the sequences are the same and they are 90 bp long" );
	}

	#[test]
	fn test_deletion_match(){
		let seq1 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATATTATTTTACTTCTACTCACAAT";
		let seq2 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT";
		let gd1 = GeneData::new( seq1, "read", "chrM", 0 );
		let gd2 = GeneData::new( seq2, "database", "chrM", 0 );
		let mut cigar = Cigar::new("");
		let mut nwa = NeedlemanWunschAffine::new(91);

		let _nw = nwa.needleman_wunsch_affine( &gd1, &gd2, 0.4);
		cigar.convert_to_cigar( &nwa.cigar_vec() );
		cigar.clean_up_cigar(&gd1, &gd2);

		assert_eq!( &cigar.to_string(), "34M34D22M", "I expected 34M34D22M as I manually deleted 34 bp from the read" );
	}

	#[test]
	fn test_insertion_match(){
		let seq2 =                                   b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATATTATTTTACTTCTACTCACAAT";
		let seq1 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT";
		let gd1 = GeneData::new( seq1, "read", "chrM", 0 );
		let gd2 = GeneData::new( seq2, "database", "chrM", 0 );
		let mut cigar = Cigar::new("");
		let mut nwa = NeedlemanWunschAffine::new(91);

		let _nw = nwa.needleman_wunsch_affine( &gd1, &gd2, 0.4 );
		cigar.convert_to_cigar( &nwa.cigar_vec() );
		cigar.clean_up_cigar(&gd1, &gd2);

		assert_eq!( &cigar.to_string(), "34M34I22M", "I expected 34M34I22M as I manually deleted 34 bp from the database" );
	}
}