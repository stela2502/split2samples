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
		let gd1 = GeneData::new( seq1, "read", "read", "chrM", 0 );
		let gd2 = GeneData::new( seq2, "database", "database", "chrM", 0 );
		let mut cigar = Cigar::new("");
		let mut nwa = NeedlemanWunschAffine::new();
		nwa.set_debug(true);

		let _nw = nwa.needleman_wunsch_affine( &gd1, &gd2, 0.4);
		cigar.convert_to_cigar( &nwa.cigar_vec() );
		cigar.clean_up_cigar(&gd1, &gd2);

		assert_eq!( &cigar.to_string(), "90M", "I expected 90M as the sequences are the same and they are 90 bp long" );
	}

	#[test]
	fn test_deletion_match(){

		//
		//
		//ACACCTAATCGGAGGAGCTACTCTAGTATTAATA----------------------------------TTATTTTACTTCTACTCACAAT
		//ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT

		let seq1 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATATTATTTTACTTCTACTCACAAT";
		let seq2 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT";
		let gd1 = GeneData::new( seq1, "read", "read", "chrM", 0 );
		let gd2 = GeneData::new( seq2, "database", "database", "chrM", 0 );
		let mut cigar = Cigar::new("");
		let mut nwa = NeedlemanWunschAffine::new();
		nwa.set_debug(true);

		let _nw = nwa.needleman_wunsch_affine( &gd1, &gd2, 0.4);
		cigar.convert_to_cigar( &nwa.cigar_vec() );
		cigar.clean_up_cigar(&gd1, &gd2);
		//assert_eq!( &cigar.to_string(), "67D23M", "I expected 34M34D22M as I manually deleted 34 bp from the read" );
		assert_eq!( &cigar.to_string(), "32M34D24M", "I expected 34M34D22M as I manually deleted 34 bp from the read" );
	}

	#[test]
	fn test_insertion_match(){
		//
		//
		//ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT
		//ACACCTAATCGGAGGAGCTACTCTAGTATTAATA--------------------------------- TTATTTTACTTCTACTCACAAT

		let seq1 = b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATAAATATTAGCCCACCAACAGCTACCATTACATTTATTATTTTACTTCTACTCACAAT";
		let seq2 =                                   b"ACACCTAATCGGAGGAGCTACTCTAGTATTAATATTATTTTACTTCTACTCACAAT";

		let gd1 = GeneData::new( seq1, "read", "read", "chrM", 0 );
		let gd2 = GeneData::new( seq2, "database", "database", "chrM", 0 );
		let mut cigar = Cigar::new("");
		let mut nwa = NeedlemanWunschAffine::new();
		nwa.set_debug(true);
		
		let _nw = nwa.needleman_wunsch_affine( &gd1, &gd2, 0.4 );
		cigar.convert_to_cigar( &nwa.cigar_vec() );
		println!("The cigar before being cleaned up in any weay: {cigar}");
		cigar.clean_up_cigar(&gd1, &gd2);
		//assert_eq!( &cigar.to_string(), "67I23M", "I expected 34M34I22M as I manually deleted 34 bp from the database" );
		assert_eq!( &cigar.to_string(), "32M34I24M", "I expected 34M34I22M as I manually deleted 34 bp from the database" );
	}
}