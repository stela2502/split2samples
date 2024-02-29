use std::collections::BTreeMap;
use crate::fast_mapper::mapper_entries::second_seq::SecondSeq;
//use crate::errors::SeqError;
//use crate::traits::BinaryMatcher;

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
	pub u8_encoded: Vec::<u8>, // the 2bit encoded array (4 times compressed)
	pub lost:usize, // how many times have I lost 4bp?
	pub shifted:usize, // how many times have we shifted our initial sequence?
	pub kmer_size:usize,
	step_size:usize,
	checker:BTreeMap::<u8, usize>,
	mask:u64, //a mask to fill matching sequences to match the index's kmer_len
	current_position:usize,

}

// Implement the Index trait for MyClass
use std::ops::Index;

impl Index<usize> for IntToStr {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        self.u8_encoded.get(index).unwrap_or(&0_u8)
    }
}


impl PartialEq<Vec<u8>> for IntToStr {
    fn eq(&self, other: &Vec<u8>) -> bool {
        &self.u8_encoded == other
    }
}


impl Iterator for IntToStr {
    type Item = (u16, SecondSeq);

    fn next(&mut self) -> Option<Self::Item> {

        let result = self.seq_at_position(self.current_position);

        if let Some((cell_id, second_seq)) = result {
            // Advance to the next position
            self.current_position += 8;
            Some((cell_id, second_seq))
        } else {
            // If the shift is 7, stop the iteration - we reached the end
            if (self.current_position % 8) + self.step_size >= 8  {
                return None
            }
            // Otherwise, move to the next position
            self.current_position = (self.current_position % 8) + self.step_size;

            // Get the result at the updated position
            let result = self.seq_at_position(self.current_position);

	        if let Some((cell_id, second_seq)) = result {
	            // Advance to the next position
	            self.current_position += 8;
	            Some((cell_id, second_seq))
	        } else {
	        	// reset the iterator for the next use
	        	self.current_position = 0;
	            None
	        }
        }
    }
}

/// Here I have my accessorie functions that more or less any of the classes would have use of.
/// I possibly learn a better way to have them...
impl IntToStr {

	pub fn new(seq:Vec::<u8>, kmer_size:usize) -> Self{
		// 4 of the array u8 fit into one result u8
		//eprintln!("Somtimes I die?! -> processed seq: {:?}", seq);
		let storage:Vec::<u8> = seq.to_vec();
		let long_term_storage = seq.to_vec();

		let u8_encoded = Vec::<u8>::with_capacity( storage.len()/4 +1);
		let lost = 0;
		let shifted = 0;
		let checker = BTreeMap::<u8, usize>::new();
		let size = match kmer_size >= 31 {
		    true => kmer_size,
		    false => 31,
		};
		let mask = !0u64 >> ( (32 - size )  * 2) as u32;
		//let mask: u64 = (1 << (2 * kmer_size)) - 1;

		let mut ret = Self{
            long_term_storage,
			storage,
			u8_encoded,
			lost,
			shifted,
			kmer_size,
			checker,
			mask,
			current_position:0,
			step_size:1,
		};

		ret.regenerate();
		ret
	}

	pub fn step_size( &mut self, size:usize ) {
		if size > 7 {
			eprintln!("step_size can not be larger than 7");
			self.step_size = 7;
		}else if size == 0 {
			eprintln!("step_size can not be zero");
			self.step_size = 1;
		}else {
			self.step_size = size;
		}
	}

    // pub fn iter(&mut self) -> SeqIterator {
    //     SeqIterator {
    //         my_class: self,
    //     }
    // }

    /// returns true if the sequence passes, false if there is a problem like N nucleotides
    /// or too low sequences variability in the first 8bp
    /// or None if the sequence is at it's end.
   /* fn seq_ok(&mut self, short:bool ) -> Result<(u16, SecondSeq), SeqError > {

        let start = self.lost * 4;
		let to = start + 8;

        self.checker.clear();

        if self.storage.len() < to{
        	//println!("not enough data!");
            return Err(SeqError::End)
        }

        // check the initial 8 bp as they will create a key
        for nuc in &self.storage[start..to ] {
            match self.checker.get_mut( nuc ){
                Some( count ) => *count += 1,
                None => {
                    if self.checker.insert( *nuc, 1).is_some(){};
                }
            };
            if *nuc ==b'N'{
            	//println!("N nucleotoides!");
                return Err(SeqError::Ns);
            }
        }

        let short_seq = SecondSeq(self.into_u64_nbp( 8 ), 8_u8);

   		let dnt_s = short_seq.di_nuc_tab();
   		let max_s = dnt_s.iter().max().unwrap_or(&0);
   		// 6 or 7 means two of the 8 nucl were not the same
   		if max_s > &5_i8 {
   			return Err(SeqError::LowComplexity)
   		}
   		if short {
   			return Ok( (short_seq.0 as u16, SecondSeq( 0_u64, 32_u8)) )
   		}

        match self.drop_n(2){// shift 8 bp
        	Some(_) => {},
        	None => {
        		return Err(SeqError::End);
        	}
        };
        let mut sign:u8 = self.kmer_size.try_into().unwrap();
        
        // do we have some sequence left to create a SecondSeq object from?
        // And if how much?
		if self.shifted + self.lost *4 + self.kmer_size > self.storage.len(){
			let missing = (self.shifted + self.lost *4 + self.kmer_size )- self.storage.len() ;
			if missing >= self.kmer_size{
				return Err(SeqError::End)
			}
        	sign = (self.kmer_size - missing).try_into().unwrap();
        }
        let long_seq = SecondSeq( self.into_u64_nbp( self.kmer_size ), sign);
        let dnt_l = long_seq.di_nuc_tab();
   		let max_l = dnt_l.iter().max().unwrap_or(&0);
   		// 6 or 7 means two of the 8 nucl were not the same
   		if *max_l as u8 > ( sign - 5_u8) {
   			return Err(SeqError::LowComplexity)
   		}

        return Ok( (short_seq.0 as u16, long_seq) )

    }


	pub fn next_small(&mut self) -> Option<(u16, SecondSeq)> {

    	//println!("Start with next");

    	match self.seq_ok( true ){ // only small!
    		Ok((short,second) ) => {
    			Some((short,second))
    		},
    		Err( SeqError::End ) => {
    			if self.shifted < 15{
    				match self.shift(){
    					Some(_) => {
    						self.next_small()
    					},
    					None => { None },
    				}
    			}else {
    				None
    			}
    		},
    		Err(_)=>{
    			//SeqError::Ns or SeqError::LowComplexity
    			self.next_small()
    		},
    	}
    }
    */

    pub fn print_second_seq (&self, first:u16, second_seq: SecondSeq ){
    	let mut small_seq :String = Default::default();
        let mut large_seq: String = Default::default();
        self.u16_to_str( 8, &first, &mut small_seq);
        self.u64_to_str( second_seq.1.into(), &second_seq.0, &mut large_seq);
        eprintln!("seq  {}-{}[{}] or {:?}-{:?}", &first, &second_seq.0, &second_seq.1, small_seq, large_seq);
    }


    // This sfucntion returen sither a Some(CellId10x, SecondSeq) - an enhanced u16 and an enhanced u64
    // or none if the u64 would only consist of less than 4 bytes as the remaining sequence is to short.
    pub fn seq_at_position(&self, start:usize ) -> Option<( u16, SecondSeq )> {
    	// Ensure the start position is within bounds
	    if start >= self.storage.len() {
	        return None;
	    }
    	
    	// Determine the range of bytes to consider
    	let start_id = start / 4;
    	let stop:u8;
    	let shift_amount = (start % 4) * 2;

    	// keep the later check from getting negative
    	if self.storage.len() <= 8 + start {
    		return None
    	}

    	let stop_id = if start + 8 + self.kmer_size + shift_amount  < self.storage.len() {
    		stop = self.kmer_size as u8;
		    start_id + 2 + (self.kmer_size + 3 + shift_amount) / 4
		} else {
			//println!("{start} + 8 + {} +3 + {shift_amount}  < {}",self.kmer_size, self.storage.len() );
			stop = (self.storage.len() - 8 - start) as u8;
			//println!( "{stop} = ({} - 8 - {start}) as u8;", self.storage.len());
		    self.u8_encoded.len()
		};

		let min_length = self.kmer_size.min(20);
		
    	if (stop as usize) < min_length  {
    		//println!("seq_at_position - to small u64: {stop} < {min_length}");
    		return None
    	}


    	//println!("I got start_id {start_id} and stop_id {stop_id}");
    	// Create a u128 from the selected byte range
	    let mut u128_value = 0u128;
	    for byte in self.u8_encoded[start_id..stop_id].iter().rev() {
	        u128_value <<= 8;
	        u128_value |= *byte as u128;
	    }

	    //println!( "and I got the the u128 {u128_value:b}");
	    // Shift the u128 value to the right position
	    
	    let shifted_value = u128_value >> shift_amount;
	    
	    //println!("I got the shift amount {shift_amount} [bits] and the u64 len {stop}");

	    // Extract u16 and u64 values
	    let u16_value = (shifted_value & u16::MAX as u128) as u16;
	    let u64_value = (shifted_value >> 16) as u64;

	    //println!("From that u128 I created \nu16: {u16_value:b} and\nu64 {u64_value:b}");

	    // Return the result
    	Some( (u16_value, SecondSeq(u64_value, stop)) )
    }


	pub fn len(&self) -> usize{
		self.u8_encoded.len()
	}

	/// this is needed for the initial mappers
	/// takes the self encoded seqences and converts the last entries into one u16
	pub fn into_u16(&self ) -> u16{

		let mut ret: u16 = 0;

    	//for i in self.len()-1..=self.len()-2 {
    	for i in (0..2).rev(){
	        ret = (ret << 8) | self.u8_encoded.get(i).copied().unwrap_or(0) as u16;
    	}

    	ret
	}	
	/// this is needed for completing the set
	pub fn into_u32(&self ) -> u32{

		let mut ret: u32 = 0;
    	//for i in self.len()-1..=self.len()-4 {
    	for i in (0..4).rev(){
    		ret = (ret << 8) | self.u8_encoded.get(i).copied().unwrap_or(0) as u32;
    	}

    	ret
	}
	/// needed for the secundary mappings
    /// takes the UTF8 encoded sequence and encodes the first 32 into a u64 
	pub fn into_u64(&self ) -> u64{

		let mut ret: u64 = 0;
    	//for i in self.len()-1..=self.len()-8{
    	for i in (0..8).rev(){
	        ret = (ret << 8) | self.u8_encoded.get( i ).copied().unwrap_or(0) as u64;
    	}

    	ret
	}
	/// this is needed for completing the set
	pub fn into_u128(&self ) -> u128{

		let mut ret: u128 = 0;
    	//for i in self.len()-1..=self.len()-16 {
    	for i in (0..16).rev(){
	        ret = (ret << 8) | self.u8_encoded.get(i).copied().unwrap_or(0) as u128;
    	}

    	ret
	}

    /// needed for the secundary mappings
    /// takes the UTF8 encoded sequence and encodes the first kmer_len into a u64 
    pub fn into_u64_nbp(&self, kmer_size:usize ) -> u64{

        let mut ret: u64 = 0;
        let target = kmer_size / 4;
        let residual = kmer_size % 4;

        if residual != 0 {
        	if self.u8_encoded.len() <= target {
        		ret = ret << 8;
        	}
        	else {
        		//println!( "have {} and want {}", self.u8_encoded.len() , target);
        		let mut loc = self.u8_encoded[target];
	        	//println!("\nGet the loc from self.u8_encoded[{}].clone()", target );
	        	//println!("I got this u64 for my u8: {loc:b}, {:b}", self.u8_encoded[target]);
	        	match residual{
	        		0 => (),
	        		3 => loc &= 0b00000011,
	        		2 => loc &= 0b00001111,
	        		1 => loc &= 0b00111111,
	        		_ => panic!("what are you trying to atchieve here? {residual}"),
	        	};
	        	ret = (ret << 8) | loc as u64;
        	}
        	
        }

        for i in (0..target).rev() {
            ret = (ret << 8) | self.u8_encoded.get( i ).copied().unwrap_or(0) as u64;
        }

        

        ret
    }


	pub fn encode_binary( &self, c: u8) -> Base {
    	// might have to play some tricks for lookup in a const
    	// array at some point
    	match c {
        	b'A' | b'a' => A,
        	b'C' | b'c' => C,
	        b'G' | b'g' => G,
        	b'T' | b't' => T,
        	b'N' | b'n' => A, // this is necessary as we can not even load a N containing sequence
	        _ => panic!("cannot decode {c} into 2 bit encoding"),
    	}
	}

	/// drop the last n u8 2 bit encoded sequences (n=1 => 4bp dropped from encoded)
	/// this can be reset using the self.reset() function.
	pub fn drop_n (&mut self, n:usize ) -> Option<()>{
		//let mut  removed:String = "".to_string();
		//self.print();
		if self.u8_encoded.len() < 2{
			return None;
		}
		for _i in 0..n{
			let _a =self.u8_encoded.remove(0);
			//removed.clear();
			//self.u8_to_str( 4, &a, &mut removed);
			//eprintln!("drop_n {i} has dropped {:?}",removed);
		}
		self.lost +=n;
		//self.print();
		Some ( () )
	}

	/// shift the underlying sequence by 1 (frameshift)
	/// can be regenrated using deep_refresh().
	pub fn shift(&mut self )-> Option<()> {
		if self.storage.len() > 1{
			let tmp = self.shifted +1;
			self.storage.remove(0);
			self.regenerate();
			self.shifted = tmp;
		}
		else {
			eprintln!("My sequence is to short to be shifted!");
			return None;
		}
		Some(())
	}


	/// regenerates the complete object with a new Vec::<u8> utf8 encoded
	pub fn from_vec_u8( &mut self, array:Vec::<u8> ) {
		// 4 of the array u8 fit into one result u8
		self.storage = array.to_vec();
		self.long_term_storage = array.to_vec();
		self.current_position=0;
		self.regenerate();
	}

	/// regenerate from the long_term_storage
	pub fn deep_refresh(&mut self ){
		self.storage = self.long_term_storage.to_vec();
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

		let mut ret = Vec::<u8>::with_capacity( target/4 );
		let mut current_byte_id:usize;

		for id in 0..(target/4) {
		//for i in 0..array.len(){
			//let val:u8=0;
			ret.push(0_u8);
			for add in (0..4).rev(){
				ret[id] <<=2;
				current_byte_id = id*4 + add;
				if self.storage.len()> current_byte_id {
					//println!("Trying to push to result {id} / {add} the entry {current_byte_id} ({})",self.storage[current_byte_id] );
					ret[id] |= self.encode_binary(self.storage[current_byte_id]);
				}
				else {
					//println!("Pushing an A instead");
					ret[id] |= self.encode_binary(b'A');
				}
			}
		}

		//self.u8_encoded = ret.iter().rev().copied().collect();
		self.u8_encoded = ret;
		self.lost = 0;
		self.shifted = 0;
	}

	/// <unsigned_int>_to_str functions do exactly that.
	/// they regenerate the initial utf8 encoded String from a IntToStr encoded integer.
	/// u64_to_str - convert a IntToStr encoded u64 to utf8 Vec::<u8> with a length of up to 32bp
	pub fn u64_to_string ( &self, kmer_size:usize, km:&u64) ->String {

	    let array = km.to_le_bytes();
	    let mut data: String = "".to_string();

	    self.u8_array_to_str( kmer_size, array.to_vec(), &mut data);
	    data
	}

	/// <unsigned_int>_to_str functions do exactly that.
	/// they regenerate the initial utf8 encoded String from a IntToStr encoded integer.
	/// u64_to_str - convert a IntToStr encoded u64 to utf8 Vec::<u8> with a length of up to 32bp
	pub fn u64_to_str ( &self, kmer_size:usize, km:&u64,  data:&mut String ){

	    let array = km.to_le_bytes();

	    self.u8_array_to_str( kmer_size, array.to_vec(), data)
	}
	pub fn u16_to_str ( &self, kmer_size:usize, km:&u16,  data:&mut String ){

	    let array = km.to_le_bytes();

	    self.u8_array_to_str( kmer_size, array.to_vec(), data)
	}



	/// Convert any integer value obtained from any of the self.into_u16.to_le_bytes() like functions
	/// and adds them at the end of the mutable data string.
	pub fn to_string (&self, bases:usize, data:&mut String ){
		let mut i = 0;
	    
		for u8_4bp in self.u8_encoded.iter(){
			i += 4;
			if i >= bases {
				//eprintln!("decoding {} bits of this number: {:b}", bases +4 - i, u8_4bp);
				self.u8_to_str( bases +4 - i, u8_4bp, data );
				break;
			}else {
				//println!("decoding 4 bits of this number: {:b}", u8_4bp);
				self.u8_to_str( 4, u8_4bp,  data );
			}
		}
	}

	pub fn u8_to_str ( &self, kmer_size:usize, u8_rep:&u8,  data:&mut String ){

		let mut loc:u8 = *u8_rep;
		//println!("converting u8 {loc:b} to string with {kmer_size} bp.");
		for _i in 0..kmer_size{
			let ch = match loc & 0b11{
	            0 => "A", //0b00
	            1 => "C", // 0b01
	            2 => "G", // 0b10
	            3 => "T", // 0b11
	            _ => "N",
	       	};
	       	*data += ch;
	       	loc >>= 2;
	       	//println!("{ch} and loc {loc:b}");
		}

		//println!("\nMakes sense? {:?}", data);
	}

	pub fn u8_array_to_str ( &self, kmer_size:usize, u8_rep:Vec::<u8>,  data:&mut String ){

		let mut i = 0;

		for u8_4bp in &u8_rep {
			i += 4;
			if i >= kmer_size {
				//println!("decoding {} bits of this number: {:b}", kmer_size - (i-4), u8_4bp);
				self.u8_to_str( kmer_size +4 -i, u8_4bp, data );
				break;
			}else {
				//println!("decoding 4 bits of this number: {:b}", u8_4bp);
				self.u8_to_str( 4, u8_4bp,  data );
			}
		}
	}

	/// This will mask the providied u64 with the internal mask definitivels 'killing' the last bp.
	/// This will convert the masked bp into 'A' or better to say #b00 entries.
	pub fn mask_u64( &self, seq:&u64 ) ->u64{
		//let filled_sed = seq | (!self.mask & (0b11 << (2 * self.kmer_size )));
		//println!("I have the mask {:b}", self.mask);
        seq & self.mask
	}

	/// This will mask the last bp by 'A' or #b00 and only return bp length usabel info.
	pub fn mask_u64_to ( &self, seq:&u64, bp:usize) -> u64  {
		if bp >= 32{
			return *seq;
		}
		let mask = !0u64 >> ( (32 - bp )  * 2) as u32;
		seq & mask
	}


	pub fn print(&self) {
		println!(">seq (n={} [{}*4 + {}])\n{}", self.storage.len(),  self.storage.len()/4,self.storage.len() %4, std::str::from_utf8( &self.storage.clone() ).expect("Invalid UTF-8"));
		println!("I have 'lost' {} u8's from my internal encoded", self.lost );
		print!("[ ");
		for en in &self.u8_encoded{
			print!(", {en:b}");
		}
		println!("]");
		println!("human readable: {:?}", self.u8_encoded );
		println!("binary vecor length {}", self.len());
	}

}
