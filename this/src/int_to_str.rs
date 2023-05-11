// logics copied from https://github.com/COMBINE-lab/kmers/
pub type Base = u8;
pub const A: Base = 0;
pub const C: Base = 1;
pub const G: Base = 2;
pub const T: Base = 3;

struct IntToStr {
	not_needed: usize,
}

/// Here I have my accessorie functions that more or less any of the classes would have use of.
/// I possibly learn a better way to have them...
impl IntToStr {

	pub fn new() -> Self{
		Self{
			not_needed:1
		}
	}

	/// Encode a UTF8 encoded string array [u8] as binary nucleotides
	/// This returns a binary coded array of u8's - in case you have more nucs to encode...
	pub fn encode2bit_u8( array:Vec::<u8> ) -> Vec::<u8> {
		// 4 of the array u8 fit into one result u8
		let mut ret = Vec::<u8>::with_capacity( array.len()/4+1);

		for i in 0..array.len(){
			if i % 4 == 0{ // need to push a new entry	
				let mut val:u8=0;
				ret.push(val);
			}
			ret[i/4] <<=2;
			ret[i/4] |= IntToStr::encode_binary(array[i]);
		}
		ret
	}

	pub fn encode_binary(c: u8) -> Base {
    	// might have to play some tricks for lookup in a const
    	// array at some point
    	match c {
        	b'A' | b'a' => A,
        	b'C' | b'c' => C,
	        b'G' | b'g' => G,
        	b'T' | b't' => T,
	        _ => panic!("cannot decode {c} into 2 bit encoding"),
    	}
	}

	pub fn u64_to_str (kmer_size:usize, km:&u64,  data:&mut String ){
	    let mut i = 0;
	    for u8_rep in km.to_le_bytes(){
	    	IntToStr::u8_to_str( kmer_size, &u8_rep, data );
	    	i +=4;
	        if i >= kmer_size{
	            break;
	        }
	    //println!( "({})Can I get that {u8_rep} to u2 ({i})? {:?}, {:?}, {:?}, {:?}, {:?}", self.seq_len, u8_rep & 0x03,  u8_rep & 0x0C, u8_rep & 0x30 , u8_rep & 0xC0, data );
	    }
	}

	pub fn decode_vec (kmer_size:usize, km:Vec::<u8>, data:&mut String ){
		let mut i = 0;
		for u8_4bp in km{
			i += 4;
			if i >= kmer_size {
				println!("decoding {} bits of this number: {:b}", kmer_size - (i-4), u8_4bp);
				IntToStr::u8_to_str( kmer_size - i +4, &u8_4bp, data );
				break;
			}else {
				//println!("decoding 4 bits of this number: {:b}", u8_4bp);
				IntToStr::u8_to_str( 4, &u8_4bp, data );
			}
		}

	}

	pub fn u8_to_str (kmer_size:usize, u8_rep:&u8,  data:&mut String ){

		let mut loc:u8 = *u8_rep;

		for _i in 0..kmer_size{
			match loc & 0x03{
	            0 => *data +="A",
	            1 => *data +="C",
	            2 => *data +="G",
	            3 => *data +="T",
	            _ => *data +="N",
	       	};
	       	loc >>= 2;
		}
	}

}

#[cfg(test)]
mod tests {

    use crate::int_to_str::IntToStr;

     #[test]
    fn check_conversion_4bp() {

    	let seq = b"AGGC";
    	let binary = IntToStr::encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  1 ); 

    	//                          G G C 
    	assert_eq!( binary, vec![ 0b101001 ]);

    	let mut decoded:String = "".to_string();

		IntToStr::decode_vec(15, binary, &mut decoded );
		assert_eq!( decoded, "AGGC" );      
                            
    }

     #[test]
    fn check_conversion_15bp() {

    	let seq = b"AGGCTTGATAGCGAG";
    	let binary = IntToStr::encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  4 );

    	panic!("{:b} {:b} {:b} {:b}", binary[0], binary[1], binary[2], binary[3] );
		//												  A G C A
		//												0b00100100  
    	//                          G G C     T T G A     T A G C     G A G
    	assert_eq!( binary, vec![ !0b101001, !0b11111000, !0b11001001, !0b100010 ]);

    	let mut decoded:String = "".to_string();

		IntToStr::decode_vec(15, binary, &mut decoded );
		assert_eq!( decoded, "AGGCTTGATAGCGAG" );
    }

     #[test]
    fn check_conversion_1bp() {

    	let seq = b"C";
    	let binary = IntToStr::encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  1 ); 
		//												  A G C A
		//												0b00100100  
    	//                          G G C     T T G A     T A G C     G A G
    	assert_eq!( binary, vec![ 0b1 ]);

    	let mut decoded:String = "".to_string();

		IntToStr::decode_vec(1, binary, &mut decoded );
		assert_eq!( decoded, "C" );
    }

    #[test]
    fn check_conversion_oneA() {

    	let seq = b"A";
    	let binary = IntToStr::encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  1 ); 
		//												  A G C A
		//												0b00100100  
    	//                          G G C     T T G A     T A G C     G A G
    	assert_eq!( binary, vec![ 0b0 ]);

    	let mut decoded:String = "".to_string();

		IntToStr::decode_vec(1, binary, &mut decoded );
		assert_eq!( decoded, "A" );
    }


         #[test]
    fn check_conversion_4_from_15bp() {

    	let seq = b"AGGCCTGTATGA";
    	let binary = IntToStr::encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  3 ); 

    	//                          G G C 
    	assert_eq!( binary[0], 0b101001  );

    	let mut decoded:String = "".to_string();

		IntToStr::decode_vec(4, binary.clone(), &mut decoded );
		assert_eq!( decoded, "AGGC" );      
        decoded.clear();
        IntToStr::decode_vec(3, binary, &mut decoded );
		assert_eq!( decoded, "AGG" );                
    }

}