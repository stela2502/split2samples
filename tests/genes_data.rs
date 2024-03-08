// /tests/genes_data.rs


#[cfg(test)]
mod tests {

	use rustody::fast_mapper::genes_data::GenesData;
	use rustody::traits::BinaryMatcher;

	#[test]
	fn test_encode_decode(){
		let seq = b"AACCTTGGGT";
		let encoded = GenesData::encode( seq );

		assert_eq!( encoded, vec![0b00000101, 0b11111010, 0b10110000], "encoded as expected {encoded:?}?" );

		if let Some(decoded) = GenesData::decode( &encoded, 10 ){
			assert_eq!( decoded, seq, "encoded as expected {seq:?}?" );
		}else {
			panic!("decode has failed!");
		}
		

		if let Some(decoded_longer) = GenesData::decode( &encoded, 12 ){
			assert_eq!( decoded_longer, b"AACCTTGGGTAA", "encoded as expected {seq:?}?" );
		}else {
			panic!("decode has failed #2!");
		}

		if let Some(_decoded_longer) = GenesData::decode( &encoded, 15 ){
			panic!("decode should have failed - too little data!");
		}else {
			assert!(true) // just to check this as OK.
		}
		
	}

	#[test]
	fn test_di_nuc_tab(){
		let seq = b"AACCTTGGGT";
		println!("AACCTTGGGT");
		let obj = GenesData::from_bytes( seq );

		let di = obj.di_nuc_tab();

		let mut exp=vec![0;16];
		exp[ 0b0 ] = 1; // AA
		exp[ 0b1 ] = 1; // AC
		exp[ 0b101 ] = 1; // CC
		exp[ 0b111 ] = 1; // CT
		exp[ 0b1111 ] = 1; // TT
		exp[ 0b1110 ] = 1; // TG
		exp[ 0b1010] = 2; // GG
		exp[ 0b1011 ] = 1; // GT

		assert_eq!( di, exp, "di_nuc_tab as expected {di:?}?" );
	}


	#[test]
	fn test_tri_nuc_tab(){
		let seq = b"AACCTTGGGT";
		println!("AACCTTGGGT");
		let obj = GenesData::from_bytes( seq );

		let di = obj.tri_nuc_tab();

		let mut exp=vec![0;64];
		exp[ 0b1 ] = 1; // AAC
		exp[ 0b101 ] = 1; // ACC
		exp[ 0b10111 ] = 1; // CCT
		exp[ 0b11111 ] = 1; // CTT
		exp[ 0b111110 ] = 1; // TTG
		exp[ 0b111010 ] = 1; // TGG
		exp[ 0b101010 ] = 1; // GGG
		exp[ 0b101011] = 2; // GGT
		exp[ 0b101011 ] = 1; // GGT

		assert_eq!( di, exp, "di_nuc_tab as expected {di:?}?" );
	}

	#[test]
	fn test_slice(){
		// to bp each
		let seq = b"AAAAAAAAAACCCCCCCCCCTTTTTTTTTTGGGGGGGGGG";
		let obj = GenesData::from_bytes( seq );
		let a_s = obj.slice(0,10);
		if let Some(content) = a_s.to_bytes(10){
			//println!("I got something here: {a_content:?}");
			assert!( content.len() == 10, "exactly 10 bp ({})", content.len());
			assert_eq!( content, b"AAAAAAAAAA", "And all are A's");
		}else {
			panic!("to_bytes did not return any value");
		}
		

		let b_s = obj.slice(5,15);
		if let Some(content) = b_s.to_bytes(10){
			//println!("I got something here: {a_content:?}");
			assert!( content.len() == 10, "exactly 10 bp ({})", content.len());
			assert_eq!( content, b"AAAAACCCCC", "5 A's and 5 C's");
		}else {
			panic!("to_bytes did not return any value");
		}

		let b_s = obj.slice(20,30);
		if let Some(content) = b_s.to_bytes(10){
			//println!("I got something here: {a_content:?}");
			assert!( content.len() == 10, "exactly 10 bp ({})", content.len());
			assert_eq!( content, b"TTTTTTTTTT", "10 T's");
			assert_eq!( b_s.get_start(), 20, "Start set to 20 ({})", b_s.get_start() );
		}else {
			panic!("to_bytes did not return any value");
		}

	}

	#[test]
	fn test_needleman_wunsch(){
		// to bp each
		let seq = b"TGGTATCTTTTACTTACCTGCTTGAATACTTG";
		let seq2 = b"TGGTATACTTACCTTTTCTGCTTGAATACTTG";
		let obj = GenesData::from_bytes( seq );
		let obj2 = GenesData::from_bytes( seq2 );

		let val =  obj.needleman_wunsch( &obj2, 0.6);
		assert!( ((val > 0.40) & ( val < 0.41 ) ), "Val close enough at 0.40625 ({val})" );
		
	}

}
