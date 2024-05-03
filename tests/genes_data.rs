// /tests/genes_data.rs


#[cfg(test)]
mod tests {

	use rustody::genes_mapper::gene_data::GeneData;
	use rustody::traits::BinaryMatcher;
	use rustody::int_to_str::IntToStr;
	use rustody::genes_mapper::Cigar;
	use rustody::genes_mapper::NeedlemanWunschAffine;

	#[test]
	fn test_encode_decode(){
		let seq = b"AACCTTGGGT";
		let encoded = GeneData::encode( seq );

		let tool=IntToStr::new( seq.to_vec(), 32);

		assert_eq!( encoded, tool.u8_encoded, "encoded as expected {encoded:?}?" );

		assert_eq!( encoded, vec![0b01010000, 0b10101111, 0b00001110], "encoded as expected {encoded:?}?" );
		
		if let Some(decoded) = GeneData::decode( &encoded, 10 ){
			assert_eq!( decoded, seq, "encoded as expected {seq:?}?" );
		}else {
			panic!("decode has failed!");
		}
		

		if let Some(decoded_longer) = GeneData::decode( &encoded, 12 ){
			assert_eq!( decoded_longer, b"AACCTTGGGTAA", "encoded as expected {seq:?}?" );
		}else {
			panic!("decode has failed #2!");
		}

		if let Some(_decoded_longer) = GeneData::decode( &encoded, 15 ){
			panic!("decode should have failed - too little data!");
		}else {
			assert!(true) // just to check this as OK.
		}
		
	}

	#[test]
	fn test_di_nuc_tab(){
		//          0000111122
		let seq = b"AACCTTGGGT";
		println!("AACCTTGGGT");
		let obj = GeneData::from_bytes( seq );

		let di = obj.di_nuc_tab();

		let mut exp=vec![0;16];
		exp[ 0b0 ] = 1; // AA
		exp[ 0b100 ] = 1; // CA
		exp[ 0b101 ] = 1; // CC
		exp[ 0b111 ] = 1; // CT
		exp[ 0b1111 ] = 1; // TT
		exp[ 0b1011 ] = 1; // TG
		exp[ 0b1010] = 2; // GG
		exp[ 0b1110 ] = 1; // GT

		assert_eq!( di, exp, "di_nuc_tab as expected {di:?}?" );
	}

	#[test]
	fn test_key_at_position(){
		//          0000111122
		let seq = b"CGTGTGTGGCGTCGTG";
		println!("CGTGTGTGGCGTCGTG");
		let obj = GeneData::from_bytes( seq );
		assert_eq!( obj.as_dna_string(), "CGTGTGTGGCGTCGTG".to_string(), "sequence correct");

		fn key_to_string( key:&u16 ) -> String{
	        let mut data = String::new();
	        for i in 0..8 {
	            let ch = match (key >> (i * 2)) & 0b11 {
	                0b00 => "A",
	                0b01 => "C",
	                0b10 => "G",
	                0b11 => "T",
	                _ => "N",
	            };
	            data += ch;
	        }
	        data
		}

		let kmers = vec![ "CGTGTGTG", "GTGTGTGG", "TGTGTGGC", "GTGTGGCG",
			"TGTGGCGT", "GTGGCGTC", "TGGCGTCG", "GGCGTCGT" ];

		for (id, kmer) in kmers.iter().enumerate() {
			match obj.key_at_position(id){
				Some( (key, start) ) => {
					assert_eq!( key_to_string(&key), kmer.to_string(), "key {}",kmer);
					assert_eq!( start, id, "start key 0");
				},
				None=> {
					panic!("iteration {id} - no kmer/key available!?");
				}
			}
		}
	}

	#[test]
	fn test_get_nucleotide_2bit(){
		let seq = b"CATGTATGGCGTCGTG";
		let obj = GeneData::from_bytes( seq );
		println!("testing this DNA seq: CATGTATGGCGTCGTG");
		//               0123 4567 8901 2345
		//               CGTG TGTG GCGT CGTG
		for (i, exp) in "CATGTATGGCGTCGTG".chars().enumerate(){
			let ch = match obj.get_nucleotide_2bit(i) {
	                Some(0b00) => 'A',
	                Some(0b01) => 'C',
	                Some(0b10) => 'G',
	                Some(0b11) => 'T',
	                Some(4_u8..=u8::MAX) => panic!("All those values must not be returned here!"),
	                None => 'N',
	    	};
	    	assert_eq!( ch, exp , "sequence at position {i} is expecetedd to be {}", exp );
		}
		
	}

	#[test]
	fn test_tri_nuc_tab(){
		//          0000111122
		let seq = b"AACCTTGGGT";
		println!("AACCTTGGGT");
		let obj = GeneData::from_bytes( seq );

		let di = obj.tri_nuc_tab();

		let mut exp=vec![0;64];

		exp[0b10000] = 1;  // AAC
		exp[0b10100] = 1;  // ACC
		exp[0b110101] = 1;  // CCT
		exp[0b111101] = 1;  // CTT
		exp[0b101111] = 1;  // TTG
		exp[0b101011] = 1;  // TGG
		exp[0b101010] = 1;  // GGG
		exp[0b111010] = 1;  // GGT

		assert_eq!( di, exp, "di_nuc_tab as expected {di:?}?" );
	}

	#[test]
	fn test_tri_nuc_tab2(){
		//          0000111122
		let seq = b"ACTGTGCA";
		println!("ACTGTGCA");
		let obj = GeneData::from_bytes( seq );

		let di = obj.tri_nuc_tab();

		let mut exp=vec![0;64];
		exp[0b110100] = 1; // ACT
		exp[0b101101] = 1; // CTG
		exp[0b111011] = 1; // TGT
		exp[0b101110] = 1; // GTG
		exp[0b11011] = 1; // TGC
		exp[0b110] = 1; // GCA

		assert_eq!( di, exp, "di_nuc_tab as expected {di:?}?" );
	}

	#[test]
	fn test_slice(){
		// to bp each
		//          0000111122223333444455556666777788889999
		let seq = b"AAACAAACCCCACCCACCCATTTGTTTGTTTGGGGTGGGT";
		let obj = GeneData::from_bytes( seq );
		if let Some(a_s) = obj.slice(0,10){
			if let Some(content) = a_s.to_bytes(10){
				//println!("I got something here: {a_content:?}");
				assert!( content.len() == 10, "exactly 10 bp ({})", content.len());
				assert_eq!( content, b"AAACAAACCC", "And all are A's");
			}else {
				panic!("to_bytes did not return any value");
			}
		}else {
			panic!("to_bytes did not return any value");
		}
		
		
		println!("Next");

		if let Some(b_s) = obj.slice(5,10){
			if let Some(content) = b_s.to_bytes(10){
				//println!("I got something here: {a_content:?}");
				assert!( content.len() == 10, "exactly 10 bp ({})", content.len());
				assert_eq!( content, b"AACCCCACCC", "5 A's and 5 C's");
			}else {
				panic!("to_bytes did not return any value");
			}
		}else {
			panic!("to_bytes did not return any value");
		}
		

		println!("Next");

		if let Some(b_s) = obj.slice(20,10){
			if let Some(content) = b_s.to_bytes(10){
				//println!("I got something here: {a_content:?}");
				assert!( content.len() == 10, "exactly 10 bp ({})", content.len());
				assert_eq!( content, b"TTTGTTTGTT", "10 T's");
				assert_eq!( b_s.get_start(), 20, "Start set to 20 ({})", b_s.get_start() );
			}else {
				panic!("to_bytes did not return any value");
			}
		}else {
			panic!("to_bytes did not return any value");
		}
		

	}


	#[test]
	fn test_next(){
		let seq = b"TGGTATCTTTTACTTA";
		let mut obj = GeneData::from_bytes( seq );
		let mut id = 0;

		/*
		A= 0b00; C=0b01, G=0b10 T= 0b11

		TGGTATCTTTTACTTA -> TGGT, ATCT, TTTA, CTTA -> 
		     T G G T       T C T A     A T T T     A T T C  
		[ "0b11101011", "0b11011100", "0b111111", "0b111101" ]
		[ "0b11101011", "0b11011100", "0b111111", "0b111101"]

		
		// encoded would be that
		//                                           1111         11111111
		//                   33221100    77665544    11009988     55443322  
		//               [ 0b11101011, 0b11011100, 0b00111111,  0b00111101]
	        1, 2, 3, 4
	        5, 6, 7, 8?  
	    */ 

	    /*
		e.g. the 9th slice (next call):
          111111    1111
          443322    11009988    77
		0b111101  0b00111111, 0b11 -> 0b1111010011111111
	    */
		let exp: Vec::<u16> = vec![ 0b1101110011101011, 0b0011110100111111, 
			0b1111011100111010, 0b1111110111001110, 0b1111111101110011, 0b0011111111011100,
			0b0100111111110111, 0b1101001111111101, 0b1111010011111111
		];

		let encoded: Vec::<u8> = vec![ 0b11101011, 0b11011100, 0b111111, 0b111101];
		assert_eq!( obj.get_encoded(), &encoded, "u8_encoded \n{:?}\nwas not the expected {:?}", 
			obj.get_encoded().iter().map(|&x|  format!("{:b}", x)).collect::<Vec<_>>(), 
			encoded.iter().map(|&x|  format!("{:b}", x)).collect::<Vec<_>>() );

		while let Some(slice) = obj.next() {
			if id == 10 {
				panic!("at pos {id} I got the slice {slice:?} {:b} vs {:b}", slice.0, exp[id]);
			}
			assert_eq!( slice.0, exp[id], "position {id} failed {:b} vs {:b}", slice.0, exp[id] );
			id +=1;
		}
		assert_eq!( id, 9, "9 possible iterations 8bp kmers for 16bp?");
	}

	#[test]
	fn test_needleman_wunsch(){
		// to bp each
		let seq = b"TGGTATCTTTTACTTACCTGCTTGAATACTTG";
		let seq2 = b"TGGTATACTTACCTTTTCTGCTTGAATACTTG";
		let obj = GeneData::from_bytes( seq );
		let obj2 = GeneData::from_bytes( seq2 );

		let val =  obj.needleman_wunsch( &obj2, 0.6, None);
		assert!( ((val > 0.40) & ( val < 0.41 ) ), "Val close enough at 0.40625 ({val})" );
		
	}

	#[test]
	fn test_cigar(){
		//                          ***
		//           012345678901234567890123456789012345678901234
		//           TGGTATCTTTTACTTACCTGCTTGAATACTTG
		//           TGGTATCTTTTACTT   TGCTTGAATACTTG
		let seq  = b"TGGTATCTTTTACTTACCTGCTTGAATACTTG";
		let seq2 = b"TGGTATCTTTTACTTTGCTTGAATACTTG";
		let obj  = GeneData::from_bytes( seq  );
		let obj2 = GeneData::from_bytes( seq2 );
		let mut cigar = Cigar::new("");
		let _val =  obj.needleman_wunsch( &obj2, 0.6, Some(&mut cigar) );

		//assert_eq!( cigar, "15M1X3D12M", "Cigar string was created correctly!" );

		//let mut cig = Cigar::new( "15M1X3D12M" );

		//cig.clean_up_cigar( &obj, &obj2 );


		assert_eq!( cigar.cigar, "15M3D14M", "Cigar string was created correctly!" );


	}

	#[test]
	fn test_needleman_wunsch_cigar(){
		// to bp each                   ***
		//           012345678901234567890123456789012345678901234
		//           TTCATATCGACAATTAGGGTTTACGACCTCGATGTTGGATCAGGT
		//           TTCATATCGACAATTAGGG   ACGACCTCGATGTTGGATCAGGT
		let seq =  b"TTCATATCGACAATTAGGGTTTACGACCTCGATGTTGGATCAGGT";
		let seq2 = b"TTCATATCGACAATTAGGGACGACCTCGATGTTGGATCAGGT"; //missing three T-s in the middle.
		let obj = GeneData::from_bytes( seq );
		let obj2 = GeneData::from_bytes( seq2 );
		let mut cigar = Cigar::new("");
		let _val =  obj.needleman_wunsch( &obj2, 0.6, Some(&mut cigar) );

		assert_eq!( cigar.cigar, "19M3D23M", "Cigar string was created correctly!" );
		
	}

	#[test]
	fn test_needleman_wunsch_cigar_insertion_and_first_mismatch(){
		// to bp each                      *** insertion
		//           0123456789012345678901   2345678901234567890123456789
		//           01234567890123456789012345678901234567890123456789
		//           GTCATATCGACAATTAGGGTTT   ACGACCTCGATGTTGGATCAGGA 22 + 22 bp
		//           TTCATATCGACAATTAGGGTTTGTGACGACCTCGATGTTGGATCAGGA
		let seq =  b"GTCATATCGACAATTAGGGTTTACGACCTCGATGTTGGATCAGGT";
		let seq2 = b"TTCATATCGACAATTAGGGTTTGTGACGACCTCGATGTTGGATCAGGT"; //missing three T-s in the middle.
		let obj = GeneData::from_bytes( seq );
		let obj2 = GeneData::from_bytes( seq2 );
		let mut cigar = Cigar::new("");
		let mut nwa = NeedlemanWunschAffine::new(40);
		let _nw = &nwa.needleman_wunsch_affine( &obj2, &obj, 0.4 );
		cigar.convert_to_cigar( &nwa.cigar_vec() );
		cigar.clean_up_cigar(&obj2, &obj);
		//let _val =  obj.needleman_wunsch( &obj2, 0.6, Some(&mut cigar) );

		//println!("{}", &nwa.to_string( &obj2, &obj ) );

		assert_eq!( cigar.cigar, "1X21M3I23M", "Cigar string was created correctly!" );
		
	}

	#[test]
	fn test_needleman_wunsch_cigar_insertion_and_last_mismatch(){
		// to bp each                      *** insertion
		//           01234567890123456789012345678901234567890123456789
		//           TTCATATCGACAATTAGGGTTT   ACGACCTCGATGTTGGATCAGGA 22 + 23 bp
		//           TTCATATCGACAATTAGGGTTTGTGACGACCTCGATGTTGGATCAGGA
		let seq =  b"TTCATATCGACAATTAGGGTTTACGACCTCGATGTTGGATCAGGC";
		let seq2 = b"TTCATATCGACAATTAGGGTTTGTGACGACCTCGATGTTGGATCAGGT"; //missing three T-s in the middle.
		let obj = GeneData::from_bytes( seq );
		let obj2 = GeneData::from_bytes( seq2 );
		let mut cigar = Cigar::new("");
		let _val =  obj.needleman_wunsch( &obj2, 0.6, Some(&mut cigar) );

		assert_eq!( cigar.cigar, "22M3I22M1X", "Cigar string was created correctly!" );
		
	}

	#[test]
	fn test_needleman_wunsch_cigar_insertion(){
		// to bp each                      *** insertion
		//           01234567890123456789012345678901234567890123456789
		//           TTCATATCGACAATTAGGGTTT   ACGACCTCGATGTTGGATCAGGA 22 + 23 bp
		//           TTCATATCGACAATTAGGGTTTGTGACGACCTCGATGTTGGATCAGGA
		let seq =  b"TTCATATCGACAATTAGGGTTTACGACCTCGATGTTGGATCAGGT";
		let seq2 = b"TTCATATCGACAATTAGGGTTTGTGACGACCTCGATGTTGGATCAGGT"; //missing three T-s in the middle.
		let obj = GeneData::from_bytes( seq );
		let obj2 = GeneData::from_bytes( seq2 );
		let mut cigar = Cigar::new("");
		let _val =  obj.needleman_wunsch( &obj2, 0.6, Some(&mut cigar) );

		assert_eq!( cigar.cigar, "22M3I23M", "Cigar string was created correctly!" );
		
	}

	#[test]
	fn test_needleman_wunsch_cigar_insertion_last_mismatch(){
		// to bp each                      *** insertion
		//           01234567890123456789012345678901234567890123456789
		//           TTCATATCGACAATTAGGGTTT   ACGACCTCGATGTTGGATCAGGA 23 + 24 bp
		//           TTCATATCGACAATTAGGGTTTGTGACGACCTCGATGTTGGATCAGGA
		let seq =  b"TTCATATCGACAATTAGGGTTTACGACCTCGATGTTGGATCAGGG";
		let seq2 = b"TTCATATCGACAATTAGGGTTTGTGACGACCTCGATGTTGGATCAGGT"; //missing three T-s in the middle.
		let obj = GeneData::from_bytes( seq );
		let obj2 = GeneData::from_bytes( seq2 );
		let mut cigar = Cigar::new("");
		let _val =  obj.needleman_wunsch( &obj2, 0.6, Some(&mut cigar) );

		assert_eq!( cigar.cigar, "22M3I22M1X", "Cigar string was created correctly!" );
		
	}

	#[test]
	fn test_key_at() {
		let seq = b"TGGTATCTTTTACTTACCTGCTTGAATACTTG";
		let obj = GeneData::from_bytes( seq );
		let first = obj.key_at(0);
		assert_eq!(first, 56555_u16,  "key was expecetd {first}" );
	}

	
	#[test]
	fn test_slicing_2(){
		let seq = b"GTTCTACATTCTTCATGGCTACTGGATTCCATGGACTCCATGTAATTATTGGATCAACATTCCTTATTGTTTGCCTACTACGACAACTAT";
		let obj = GeneData::from_bytes( seq );
		assert_eq!( obj.get_start(), 0 ,"Start is 0");

		let slice = if let Some(slice) = obj.slice(15,20) {
			assert_eq!( slice.get_start(), 15 ,"Start is 15");
			assert_eq!( slice.len(), 20 ,"Length is 20");
			slice
		}else {
			panic!("I could not get the slice!");
		};

		if let Some(slice_2) = slice.slice(5,10){
			assert_eq!( slice_2.get_start(), 20 ,"Start is 20");
			assert_eq!( slice_2.len(), 10 ,"Start is 0");
		}
	}

	#[test]
	fn test_real_live_issue1(){
		let seq = b"CGCCATCTTCAGCAAACCCTAAAAAGGTATTAAAGTAAGCAAAAGAATCAAACATAAAAACGTTAGGTCAAGGTGTAGCCAATGAAATGG";
		let seq2 = b"CGCCATCTTCAGCAAACCCTAAAAAGGTATTAAAGTAAGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGACACAAAACGTTTGGGGG";
		let obj = GeneData::from_bytes( seq );
		let obj2 = GeneData::from_bytes( seq2 );
		let mut cigar = Cigar::new("");
		let _val =  obj.needleman_wunsch( &obj2, 0.6, Some(&mut cigar) );

		assert_eq!( cigar.cigar, "44M46S", "Cigar string was created correctly!" );
	}




}
