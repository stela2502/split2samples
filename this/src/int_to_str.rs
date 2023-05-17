// logics copied from https://github.com/COMBINE-lab/kmers/
pub type Base = u8;
pub const A: Base = 0;
pub const C: Base = 1;
pub const G: Base = 2;
pub const T: Base = 3;

#[derive(Debug,PartialEq)]
pub struct IntToStr {
	long_term_storage: Vec::<u8>, //this will never be shifted nor poped
	storage: Vec::<u8>,  // the initial utf8 encoded data
	u8_encoded: Vec::<u8>, // the 2bit encoded array (4 times compressed)
	lost:usize, // how many times have I lost 4bp?
	pub shifted:usize, // how many times have we shifted our initial sequence?
	kmer_size:usize,
}

// Implement the Index trait for MyClass
use std::ops::Index;

impl Index<usize> for IntToStr {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        &self.u8_encoded.get(index).unwrap_or(&0_u8)
    }
}


impl PartialEq<Vec<u8>> for IntToStr {
    fn eq(&self, other: &Vec<u8>) -> bool {
        &self.u8_encoded == other
    }
}




/// Here I have my accessorie functions that more or less any of the classes would have use of.
/// I possibly learn a better way to have them...
impl IntToStr {

	pub fn new(seq:Vec::<u8>, kmer_size:usize) -> Self{
		// 4 of the array u8 fit into one result u8
		let storage:Vec::<u8> = seq.iter().copied().collect();
		let long_term_storage = seq.iter().copied().collect();

		let u8_encoded = Vec::<u8>::with_capacity( storage.len()/4 +1);
		let lost = 0;
		let shifted = 0;
		let mut ret = Self{
            long_term_storage,
			storage,
			u8_encoded,
			lost,
			shifted,
			kmer_size,
		};

		ret.regenerate();
		return ret
	}

    // pub fn iter(&mut self) -> SeqIterator {
    //     SeqIterator {
    //         my_class: self,
    //     }
    // }

    pub fn next(&mut self) -> Option<(u16, u64)> {

        if self.u8_encoded.len() < 2+ self.kmer_size/4 +1{
            if self.shifted == 3{
                return None
            }else {
                self.shift();
                return self.next()
            }
        }
        let short = self.into_u16().clone();

        self.drop_n(2);
        
        let long = self.into_u64_nbp( self.kmer_size ).clone();
        println!("I will return {short} and {long}");
        Some(( short, long ))
    }

	pub fn len(&self) -> usize{
		self.u8_encoded.len()
	}

	/// this is needed for the initial mappers
	/// takes the self encoded seqences and converts the last entries into one u16
	pub fn into_u16(&self ) -> u16{

		let mut ret: u16 = 0;

    	for i in self.len()-1..=self.len()-2 {
	        ret = (ret << 8) | self.u8_encoded.get(i).copied().unwrap_or(0) as u16;
    	}

    	ret
	}	
	/// this is needed for completing the set
	pub fn into_u32(&self ) -> u32{

		let mut ret: u32 = 0;
    	for i in self.len()-1..=self.len()-4 {
    		ret = (ret << 8) | self.u8_encoded.get(i).copied().unwrap_or(0) as u32;
    	}

    	ret
	}
	/// needed for the secundary mappings
    /// takes the UTF8 encoded sequence and encodes the first 32 into a u64 
	pub fn into_u64(&self ) -> u64{

		let mut ret: u64 = 0;
    	for i in self.len()-1..=self.len()-8{
	        ret = (ret << 8) | self.u8_encoded.get( i ).copied().unwrap_or(0) as u64;
    	}

    	ret
	}

    /// needed for the secundary mappings
    /// takes the UTF8 encoded sequence and encodes the first kmer_len into a u64 
    pub fn into_u64_nbp(&self, kmer_size:usize ) -> u64{

        let mut ret: u64 = 0;
        let target = kmer_size / 4;
        let mut residual = kmer_size % 4;

        for i in self.len()-1..=self.len()-target {
            ret = (ret << 8) | self.u8_encoded.get( i ).copied().unwrap_or(0) as u64;
        }
        // if self.u8_encoded.len() > target{
        //     while residual != 0{
                
        //         residual -=1;
        //     }
        //     // push single byte pairs of the next 2bit u8 into the u64
        //     let mut last
        // }

        ret
    }

	/// this is needed for completing the set
	pub fn into_u128(&self ) -> u128{

		let mut ret: u128 = 0;
    	for i in self.len()-1..=self.len()-16 {
	        ret = (ret << 8) | self.u8_encoded.get(i).copied().unwrap_or(0) as u128;
    	}

    	ret
	}


	pub fn enc_2bit( c: u8) -> Base {
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

	/// drop the last n u8 2 bit encoded sequences (n=1 => 4bp dropped from encoded)
	/// this can be reset using the self.reset() function.
	pub fn drop_n (&mut self, n:usize ){
		for _i in 0..n{
			let _a =self.u8_encoded.pop().unwrap_or(0);
			//println!("drop_n {i} has dropped {a:b}");
			//self.print();
		}
	}

	pub fn shift(&mut self ){
		let tmp = self.shifted +1;
		self.storage.remove(0);
		self.regenerate();
		self.shifted = tmp;
	}

	/// Convert any integer value obtained from any of the self.into_u16.to_le_bytes() like functions
	/// and adds them at the end of the mutable data string.
	pub fn enc_vec_to_string (&self, array:Vec::<u8>, bases:usize, data :&mut String ){
		let mut i = 0;
		let final_whole = bases / 4;  // Determine the chunk index (0, 1, 2, ...)
	    let remainder   = bases % 4;    // Determine the remainder within the chunk (0, 1, 2, 3)
	    
		for u8_4bp in array.iter().rev(){
			i += 4;
			if i >= bases {
				//println!("decoding {} bits of this number: {:b}", kmer_size - (i-4), u8_4bp);
				self.u8_to_str( bases - i +4, &u8_4bp, data );
				break;
			}else {
				//println!("decoding 4 bits of this number: {:b}", u8_4bp);
				self.u8_to_str( 4, &u8_4bp, data );
			}
		}
	}

	/// regenerates the complete object with a new Vec::<u8> utf8 encoded
	pub fn from_vec_u8( &mut self, array:Vec::<u8> ) {
		// 4 of the array u8 fit into one result u8
		self.storage = array.iter().copied().collect();
		self.long_term_storage = array.iter().copied().collect();
		self.regenerate();
	}

	/// regenerate from the long_term_storage
	pub fn deep_refresh(&mut self ){
		self.storage = self.long_term_storage.iter().copied().collect();
		self.regenerate();
	}

	/// regenerate the encoded from the storage
	/// this can be used to re-gain lost sequences
	pub fn regenerate (&mut self ){

		let mut target = self.storage.len();
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

				if self.storage.len()> i {
					//println!("Trying to push to result {id} / {add} the entry {i} ({})",array[i] );
					ret[id] |= self.encode_binary(self.storage[i]);
				}
				else {
					//println!("Pushing an A instead");
					ret[id] |= self.encode_binary(b'A');
				}
			}
		}
		self.u8_encoded = ret.iter().rev().copied().collect();
		self.lost = 0;
		self.shifted = 0;
	}


	/// <unsigned_int>_to_str functions do exactly that.
	/// they regenerate the initial utf8 encoded String from a IntToStr encoded integer.
	/// u64_to_str - convert a IntToStr encoded u64 to utf8 Vec::<u8> with a length of up to 32bp
	pub fn u64_to_str ( &self, kmer_size:usize, km:&u64,  data:&mut String ){
	    let final_whole = kmer_size / 4;  // Determine the chunk index (0, 1, 2, ...)
	    let remainder   = kmer_size % 4;    // Determine the remainder within the chunk (0, 1, 2, 3)
	    
	    println!("The total {final_whole}");

	    let array = km.to_le_bytes();

	    let mut i = 0;
		for u8_4bp in array.iter().rev(){
			i += 4;
			if i >= kmer_size {
				//println!("decoding {} bits of this number: {:b}", kmer_size - (i-4), u8_4bp);
				self.u8_to_str( kmer_size - i +4, &u8_4bp, data );
				break;
			}else {
				//println!("decoding 4 bits of this number: {:b}", u8_4bp);
				self.u8_to_str( 4, &u8_4bp, data );
			}
		}
	}




	/// Convert any integer value obtained from any of the self.into_u16.to_le_bytes() like functions
	/// and adds them at the end of the mutable data string.
	pub fn to_string (&self, bases:usize, data:&mut String ){
		let mut i = 0;
		let final_whole = bases / 4;  // Determine the chunk index (0, 1, 2, ...)
	    let remainder   = bases % 4;    // Determine the remainder within the chunk (0, 1, 2, 3)
	    
		for u8_4bp in self.u8_encoded.iter().rev(){
			i += 4;
			if i >= bases {
				//println!("decoding {} bits of this number: {:b}", kmer_size - (i-4), u8_4bp);
				self.u8_to_str( bases - i +4, &u8_4bp, data );
				break;
			}else {
				//println!("decoding 4 bits of this number: {:b}", u8_4bp);
				self.u8_to_str( 4, &u8_4bp,  data );
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

	pub fn u8_array_to_str ( &self, kmer_size:usize, u8_rep:Vec::<u8>,  data:&mut String ){

		let mut i = 0;
		let final_whole = kmer_size / 4;  // Determine the chunk index (0, 1, 2, ...)
	    let remainder   = kmer_size % 4;    // Determine the remainder within the chunk (0, 1, 2, 3)
		
		for u8_4bp in self.u8_encoded.iter().rev(){
			i += 4;
			if i >= kmer_size {
				//println!("decoding {} bits of this number: {:b}", kmer_size - (i-4), u8_4bp);
				self.u8_to_str( kmer_size - i +4, &u8_4bp, data );
				break;
			}else {
				//println!("decoding 4 bits of this number: {:b}", u8_4bp);
				self.u8_to_str( 4, &u8_4bp,  data );
			}
		}
	}

	pub fn print(&self) {
		println!(">seq\n{}", std::str::from_utf8( &self.storage.clone() ).expect("Invalid UTF-8"));
		println!("I have 'lost' {} u8's from my internal encoded", self.lost );
		print!("[ ");
		for en in &self.u8_encoded{
			print!(", {en:b}");
		}
		println!("]");
	}

}
