// logics copied from https://github.com/COMBINE-lab/kmers/
pub type Base = u8;
pub const A: Base = 0;
pub const C: Base = 1;
pub const G: Base = 2;
pub const T: Base = 3;

#[derive(Debug,PartialEq)]
pub struct IntToStr {
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

	/// needed for the secundary mappings
    /// takes the UTF8 encoded sequence and encodes the first 32 into a u64 
	pub fn into_u64(&self,  data:Vec::<u8> ) -> u64{

		let vec = self.encode2bit_u8_sized( data.to_vec(), 32 );
		let a: u64 = *vec.get(0).unwrap_or(&0_u8) as u64;
   		let b: u64 = *vec.get(1).unwrap_or(&0_u8) as u64;
    	let c: u64 = *vec.get(2).unwrap_or(&0_u8) as u64;
    	let d: u64 = *vec.get(3).unwrap_or(&0_u8) as u64;
    	let e: u64 = *vec.get(4).unwrap_or(&0_u8) as u64;
   		let f: u64 = *vec.get(5).unwrap_or(&0_u8) as u64;
    	let g: u64 = *vec.get(6).unwrap_or(&0_u8) as u64;
    	let h: u64 = *vec.get(7).unwrap_or(&0_u8) as u64;
    	let concatenated: u64 = (a << 56) | (b << 48) | (c << 40) | (d << 32) | (e << 24) | (f << 16) | (g << 8) | h as u64;
    	//let concatenated: u64 = (a << 24) | (b << 16) | (c << 8) | d;
    	concatenated

	}

	// this is needed for the initial mappers
	// takes the utf8 encoded seqences and converts them into one u16
	pub fn into_u16(&self,  data:Vec::<u8> ) -> u16{

		let vec = self.encode2bit_u8_sized( data.to_vec(), 16 );
		if vec.len() < 2 {
			panic!("Not enough data to create one (leading) u64 from this data");
		}
		let a: u16 = u16::from(vec[0]);
   		let b: u16 = u16::from(vec[1]);
    	let concatenated: u16 =  (a << 8) | b;
    	concatenated

	}

	pub fn encode2bit_u8_sized( &self, array:Vec::<u8>, kmer_size:usize ) -> Vec::<u8> {
		let mut target = kmer_size;
		let remainder = target % 4;
	    if remainder != 0{
	    	target +=  4 - remainder;
	    }

		let mut ret = Vec::<u8>::with_capacity( target/4);
		let mut i:usize;
		for id in 0..(target/4) {
		//for i in 0..array.len(){
			//let val:u8=0;
			ret.push(0_u8);
			for add in (0..4).rev(){
				ret[id] <<=2;
				i = id*4 + add;

				if array.len()> i {
					//println!("Trying to push to result {id} / {add} the entry {i} ({})",array[i] );
					ret[id] |= self.encode_binary(array[i]);
				}
				else {
					//println!("Pushing an A instead");
					ret[id] |= self.encode_binary(b'A');
				}
			}
		}
		ret
	}

	/// Encode a UTF8 encoded string array [u8] as binary nucleotides
	/// This returns a binary coded array of u8's - in case you have more nucs to encode...
	pub fn encode2bit_u8( &self, array:Vec::<u8> ) -> Vec::<u8> {
		// 4 of the array u8 fit into one result u8

		let mut target = array.len();
		let remainder = target % 4;
	    if remainder != 0{
	    	target +=  4 - remainder;
	    }

		let mut ret = Vec::<u8>::with_capacity( target/4);
		let mut i:usize;
		for id in 0..(target/4) {
		//for i in 0..array.len(){
			//let val:u8=0;
			ret.push(0_u8);
			for add in (0..4).rev(){
				ret[id] <<=2;
				i = id*4 + add;

				if array.len()> i {
					//println!("Trying to push to result {id} / {add} the entry {i} ({})",array[i] );
					ret[id] |= self.encode_binary(array[i]);
				}
				else {
					//println!("Pushing an A instead");
					ret[id] |= self.encode_binary(b'A');
				}
			}
		}
		ret
	}
	// pub fn encode_positions(len: usize) -> Vec::<usize> {

	//     let mut chunk_index:usize;  // Determine the chunk index (0, 1, 2, ...)
	//     let mut remainder = len % 4;    // Determine the remainder within the chunk (0, 1, 2, 3)
	    
	//     let mut target = len;
	//     if remainder != 0{
	//     	target +=  4 - remainder;
	//     }
	// 	println! ("encode_positions for len {len} has the target at {target}");
	//     // Calculate the encoded value based on the chunk index and remainder
	//     let mut ret = Vec::<usize>::with_capacity(target);
	//     for index in 0..target{
	//     	chunk_index = index / 4;
	//     	remainder = index % 4;
	//     	ret.push( (3 - remainder) + chunk_index * 4 );
	//     	println! ("encode_positions at {index} got {}", (3 - remainder) + chunk_index * 4 );
	//     }

	//     ret
	// }



	pub fn encode_binary( &self, c: u8) -> Base {
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

	pub fn u64_to_str ( &self, kmer_size:usize, km:&u64,  data:&mut String ){
	    let final_whole = kmer_size / 4;  // Determine the chunk index (0, 1, 2, ...)
	    let remainder   = kmer_size % 4;    // Determine the remainder within the chunk (0, 1, 2, 3)
	    
	    let array = km.to_le_bytes();
	    for i in 0..final_whole{
	    	self.u8_to_str( 4, &array[i] , data );
	    }
	    self.u8_to_str( remainder, &array[final_whole] , data );
	}

	pub fn decode_vec ( &self, kmer_size:usize, km:Vec::<u8>, data:&mut String ){
		let mut i = 0;
		for u8_4bp in km{
			i += 4;
			if i >= kmer_size {
				println!("decoding {} bits of this number: {:b}", kmer_size - (i-4), u8_4bp);
				self.u8_to_str( kmer_size - i +4, &u8_4bp, data );
				break;
			}else {
				//println!("decoding 4 bits of this number: {:b}", u8_4bp);
				self.u8_to_str( 4, &u8_4bp, data );
			}
		}

	}

	pub fn u8_to_str ( &self, kmer_size:usize, u8_rep:&u8,  data:&mut String ){

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
    fn test_u64_to_str(){

        let num:u64 = 15561;
        //println!("This is the number I want to convert to a Sequence: {num:b}");
        // 11110011001001
        // T T A T A C G 
        let tool = IntToStr::new();
        let mut data:String = "".to_string();
        tool.u64_to_str( 7, &num, &mut data );
        assert_eq!( data, "CGATATT".to_string() )
    } 

     #[test]
    fn check_conversion_4bp() {

    	let seq = b"AGGC";
    	//		   C G G A
    	//         01101000
    	let tool = IntToStr::new();
    	let binary = tool.encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  1 ); 
    	//panic!("{:b}", binary[0] );
    	//                          G G C 
    	assert_eq!( binary, vec![ 0b1101000 ]);

    	let mut decoded:String = "".to_string();

		tool.decode_vec(15, binary, &mut decoded );
		assert_eq!( decoded, "AGGC" );      
                            
    }



     #[test]
    fn check_conversion_15bp() {
        //          0000111122223333   
    	let seq = b"AGGCTTGATAGCGAG";
    	let tool = IntToStr::new();
    	let binary = tool.encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  4 );

    	//panic!("{:b} {:b} {:b} {:b}", binary[0], binary[1], binary[2], binary[3] );
		//												  A G C A
		//												0b00100100  
    	//                          G G C     T T G A     T A G C     G A G

    	//panic!("the binaries I get: {:b} {:b} {:b} {:b} ", binary[0], binary[1], binary[2], binary[3]);
    	assert_eq!( binary, vec![ 0b1101000, 0b101111, 0b1100011, 0b100010 ]);

    	let mut decoded:String = "".to_string();

		tool.decode_vec(15, binary, &mut decoded );
		assert_eq!( decoded, "AGGCTTGATAGCGAG" );
    }

     #[test]
    fn check_conversion_1bp() {

    	let seq = b"C";
    	let tool = IntToStr::new();
    	let binary = tool.encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  1 ); 
		//												  A G C A
		//												0b00100100  
    	//                          G G C     T T G A     T A G C     G A G
    	assert_eq!( binary, vec![ 0b1 ]);

    	let mut decoded:String = "".to_string();

		tool.decode_vec(1, binary, &mut decoded );
		assert_eq!( decoded, "C" );
    }

    #[test]
    fn check_conversion_oneA() {

    	let seq = b"A";
    	let tool = IntToStr::new();
    	let binary = tool.encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  1 ); 
		//												  A G C A
		//												0b00100100  
    	//                          G G C     T T G A     T A G C     G A G
    	assert_eq!( binary, vec![ 0b0 ]);

    	let mut decoded:String = "".to_string();

		tool.decode_vec(1, binary, &mut decoded );
		assert_eq!( decoded, "A" );
    }


    #[test]
    fn check_conversion_4_from_15bp() {

    	let seq = b"AGGCCTGTATGA";
    	let tool = IntToStr::new();
    	let binary = tool.encode2bit_u8( seq.to_vec() );

    	assert_eq!( binary.len(),  3 ); 

    	//                          G G C 
    	assert_eq!( binary[0], 0b1101000  );

    	let mut decoded:String = "".to_string();

		tool.decode_vec(4, binary.clone(), &mut decoded );
		assert_eq!( decoded, "AGGC" );      
        decoded.clear();
        tool.decode_vec(3, binary, &mut decoded );
		assert_eq!( decoded, "AGG" );                
    }

}