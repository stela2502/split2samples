// /tests/genes_data.rs


#[cfg(test)]
mod tests {

	use rustody::genes_mapper::gene_data::GeneData;
	use rustody::traits::BinaryMatcher;
	use rustody::int_to_str::IntToStr;

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

		if let Some(b_s) = obj.slice(5,15){
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

		if let Some(b_s) = obj.slice(20,30){
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
	fn test_needleman_wunsch(){
		// to bp each
		let seq = b"TGGTATCTTTTACTTACCTGCTTGAATACTTG";
		let seq2 = b"TGGTATACTTACCTTTTCTGCTTGAATACTTG";
		let obj = GeneData::from_bytes( seq );
		let obj2 = GeneData::from_bytes( seq2 );

		let val =  obj.needleman_wunsch( &obj2, 0.6);
		assert!( ((val > 0.40) & ( val < 0.41 ) ), "Val close enough at 0.40625 ({val})" );
		
	}

	#[test]
	fn test_key_at() {
		let seq = b"TGGTATCTTTTACTTACCTGCTTGAATACTTG";
		let obj = GeneData::from_bytes( seq );
		let first = obj.key_at(0);
		assert_eq!(first, 56555_u16,  "key was expecetd {first}" );
	}

}
