// tests/genes.rs

#[cfg(test)]
mod tests {
	use rustody::gene::Gene;

	#[test]
	fn test_new(){
		/*
		pub struct Gene{
			pub chrom:String, // the cromosome id to look for the sequence
			pub start:usize, // the start position for this entry
			pub end:usize, // the end position for this entry
			exons:Vec<[usize;2]>, // a vector of start and end positions
			sense_strand:bool, // sense_strand in the genome true 1->n; false n <- 1
			pub name:String, // the gene symbol
			pub ids:Vec<String>, // e.g. ENSMBL ID and other entries like family name or class 
		}
		*/
		//new(chrom:String, start_s:String, end_s:String, sense_strand_s:String, name:String, ids:Vec<String> )
		let gene = Gene::new( "chr1", "1", "100", 
			"+","testTranscript1", "testGene1", vec!["testFam1".to_string()] );
		assert_eq!( gene.chrom, "chr1", "chr" );
		assert_eq!( gene.start, 1, "start" );
		assert_eq!( gene.end, 100, "end" );
		//assert_eq!( gene.sense_strand, true, "sense" );
		assert_eq!( gene.name, "testGene1", "name" );
		assert_eq!( gene.transcript, "testTranscript1", "name" );
		assert_eq!( gene.ids[0], "testFam1", "family name" );
	}
	#[test]
	fn test_exons(){
		let mut gene = Gene::new( "chr1", "1", "100", 
			"+", "testTranscript1", "testGene1", vec!["testFam1".to_string()] );
		
		gene.add_exon( "1", "20" );
		gene.add_exon( "35", "60" );

		assert_eq!( gene.exon_length(), 2, "exons have been stored");
	}

	#[test]
	fn test_to_mrna(){
		let mut gene = Gene::new( "chr1", "1", "20", 
			"+", "testTranscript1", "testGene1", vec!["testFam1".to_string()] );
		
		gene.add_exon( "1", "5" );
		gene.add_exon( "15", "20" );

		if let Some(seq) = gene.to_mrna( b"ACCCCGAAAAAAAAAATTTT", 100 ) {
			assert_eq!( seq, b"ACCCCAATTTT", "Correct mrna? {seq:?}" );
		}else {
			panic!("to_mrna has not returned a sequence!")
		}

		if let Some(seq) = gene.to_mrna( b"ACCCCGAAAAAAAAAATTT", 100 ) {
			panic!("The sequences is one bp too short - that should return None - not {seq:?}!")
		}else {
			assert!(true); // to show we did good.
		}
		
	}

	#[test]
	fn test_to_mrna_rev(){
		let mut gene = Gene::new( "chr1", "1", "20", 
			"-", "testTranscript1", "testGene1", vec!["testFam1".to_string()] );
		
		gene.add_exon( "1", "5" );
		gene.add_exon( "15", "20" );

		if let Some(seq) = gene.to_mrna( b"ACCCCGAAAAAAAAAATTTT", 100 ) {
			assert_eq!( seq, b"AAAATTGGGGT", "Correct mrna? {seq:?}" );
		}else {
			panic!("to_mrna has not returned a sequence!")
		}

		if let Some(seq) = gene.to_mrna( b"ACCCCGAAAAAAAAAATTT", 100 ) {
			panic!("The sequences is one bp too short - that should return None - not {seq:?}!")
		}else {
			assert!(true); // to show we did good.
		}
		
	}

	#[test]
	fn test_to_mrna2(){
		//          ********************--------------------********************--------------------********************--------------------*********************
		//          AGCTAAATGCAGACTGTGGC                    GTTTGACTTGGTATCTTTTC                    GATTGTGCTGAGACTGGAGT                    TCTAAGAGGATGCTGGCTAT
		//          012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
		let seq = b"AGCTAAATGCAGACTGTGGCAGAATGGAATTTCATGAGAGGTTTGACTTGGTATCTTTTCTAGGCTAAAACTACAAAGAAGATTGTGCTGAGACTGGAGTGCGTTGAGCCCAACTGCAGATCTAAGAGGATGCTGGCTATTAAGAGATGCAAGCATTTTGAATTGGGAGGCGACAAGAAGAGAAAGGGCCAAGTGATCCAGTTCTAAGCAGATTTTGTTATGAAGACAATAAAATCTTGACCTTTC";

		let mut gene = Gene::new( "chr1", "1", "200", 
			"+", "testTranscript1", "testGene1", vec!["testFam1".to_string()] );

		gene.add_exon( "1", "20" );
		gene.add_exon( "41", "60" );
		gene.add_exon( "81", "100" );
		gene.add_exon( "121", "140" );
		


		if let Some(mrna) = gene.to_mrna( seq, 100 ) {
			assert_eq!(mrna.len(), 80, "I had expected a 60bp sequence");
			assert_eq!( std::str::from_utf8(&mrna).expect("Invalid UTF-8 sequence"), "AGCTAAATGCAGACTGTGGCGTTTGACTTGGTATCTTTTCGATTGTGCTGAGACTGGAGTTCTAAGAGGATGCTGGCTAT", "spliced correctly");
		}else {
			panic!("I did not get an mRNA at all.")
		}
	}

	#[test]
	fn test_to_mrna2_sub(){
		//          ********************--------------------********************--------------------********************--------------------*********************
		//          AGCTAAATGCAGACTGTGGC                    GTTTGACTTGGTATCTTTTC                    GATTGTGCTGAGACTGGAGT                    TCTAAGAGGATGCTGGCTAT
		//          012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
		let seq = b"AGCTAAATGCAGACTGTGGCAGAATGGAATTTCATGAGAGGTTTGACTTGGTATCTTTTCTAGGCTAAAACTACAAAGAAGATTGTGCTGAGACTGGAGTGCGTTGAGCCCAACTGCAGATCTAAGAGGATGCTGGCTATTAAGAGATGCAAGCATTTTGAATTGGGAGGCGACAAGAAGAGAAAGGGCCAAGTGATCCAGTTCTAAGCAGATTTTGTTATGAAGACAATAAAATCTTGACCTTTC";

		let mut gene = Gene::new( "chr1", "1", "200", 
			"+", "testTranscript1", "testGene1", vec!["testFam1".to_string()] );

		gene.add_exon( "1", "20" );
		gene.add_exon( "41", "60" );
		gene.add_exon( "81", "100" );
		gene.add_exon( "121", "140" );
		


		if let Some(mrna) = gene.to_mrna( seq, 30 ) {
			assert_eq!(mrna.len(), 30, "I had expected a 30bp sequence");
			assert_eq!( std::str::from_utf8(&mrna).expect("Invalid UTF-8 sequence"), "AGACTGGAGTTCTAAGAGGATGCTGGCTAT", "spliced correctly - and subset");
		}else {
			panic!("I did not get an mRNA at all.")
		}
	}


	#[test]
	fn test_to_mrna2_rev(){
		//          ********************--------------------********************--------------------********************--------------------*********************
		//          AGCTAAATGCAGACTGTGGC                    GTTTGACTTGGTATCTTTTC                    GATTGTGCTGAGACTGGAGT                    TCTAAGAGGATGCTGGCTAT
		//          012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
		let seq = b"AGCTAAATGCAGACTGTGGCAGAATGGAATTTCATGAGAGGTTTGACTTGGTATCTTTTCTAGGCTAAAACTACAAAGAAGATTGTGCTGAGACTGGAGTGCGTTGAGCCCAACTGCAGATCTAAGAGGATGCTGGCTATTAAGAGATGCAAGCATTTTGAATTGGGAGGCGACAAGAAGAGAAAGGGCCAAGTGATCCAGTTCTAAGCAGATTTTGTTATGAAGACAATAAAATCTTGACCTTTC";

		let mut gene = Gene::new( "chr1", "1", "200", 
			"-", "testTranscript1", "testGene1", vec!["testFam1".to_string()] );

		gene.add_exon( "1", "20" );
		gene.add_exon( "41", "60" );
		gene.add_exon( "81", "100" );
		gene.add_exon( "121", "140" );
		

		if let Some(mrna) = gene.to_mrna( seq, 100 ) {
			assert_eq!(mrna.len(), 80, "I had expected a 60bp sequence");
			// manually reversed the sequence!
			assert_eq!( std::str::from_utf8(&mrna).expect("Invalid UTF-8 sequence"), "ATAGCCAGCATCCTCTTAGAACTCCAGTCTCAGCACAATCGAAAAGATACCAAGTCAAACGCCACAGTCTGCATTTAGCT", "spliced correctly");
		}else {
			panic!("I did not get an mRNA at all.")
		}
	}

	#[test]
	fn test_to_mrna2_rev_sub(){
		//          ********************--------------------********************--------------------********************--------------------*********************
		//          AGCTAAATGCAGACTGTGGC                    GTTTGACTTGGTATCTTTTC                    GATTGTGCTGAGACTGGAGT                    TCTAAGAGGATGCTGGCTAT
		//          012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
		let seq = b"AGCTAAATGCAGACTGTGGCAGAATGGAATTTCATGAGAGGTTTGACTTGGTATCTTTTCTAGGCTAAAACTACAAAGAAGATTGTGCTGAGACTGGAGTGCGTTGAGCCCAACTGCAGATCTAAGAGGATGCTGGCTATTAAGAGATGCAAGCATTTTGAATTGGGAGGCGACAAGAAGAGAAAGGGCCAAGTGATCCAGTTCTAAGCAGATTTTGTTATGAAGACAATAAAATCTTGACCTTTC";

		let mut gene = Gene::new( "chr1", "1", "200", 
			"-", "testTranscript1", "testGene1", vec!["testFam1".to_string()] );

		gene.add_exon( "1", "20" );
		gene.add_exon( "41", "60" );
		gene.add_exon( "81", "100" );
		gene.add_exon( "121", "140" );
		

		if let Some(mrna) = gene.to_mrna( seq, 30 ) {
			assert_eq!(mrna.len(), 30, "I had expected a 60bp sequence");
			// manually reversed the sequence!
			assert_eq!( std::str::from_utf8(&mrna).expect("Invalid UTF-8 sequence"), "CAAGTCAAACGCCACAGTCTGCATTTAGCT", "spliced correctly");
		}else {
			panic!("I did not get an mRNA at all.")
		}
	}
}