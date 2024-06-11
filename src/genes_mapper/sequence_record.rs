// sequence_record.rs

/// A really simple copy of needletails SequenceRecord to allow a hard copy of the data
/// This would be necessary for me to used
use std::default::Default;
use core::ops::Add;
//use num_traits::cast::ToPrimitive;
use std::ops::AddAssign;
//use std::collections::HashMap;
//use lazy_static::lazy_static; // Import the lazy_static macro

//use std::cmp::Ordering;
use core::fmt;

// Define the lookup table outside the impl block
/*static LOOKUP_TABLE: Lazy<HashMap<char, u8>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert('A', 0b00);
    map.insert('C', 0b01);
    map.insert('G', 0b10);
    map.insert('T', 0b11);
    map.insert('a', 0b00); // Insert lowercase nucleotides as well
    map.insert('c', 0b01);
    map.insert('g', 0b10);
    map.insert('t', 0b11);
    map
});*/

#[derive(Debug, Clone)]
pub struct SeqRec{
	id: Vec<u8>,
	seq: Vec<u8>,
	qual:Vec<u8>,
}

// Implementing Display trait for SeqRec
impl fmt::Display for SeqRec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "@{}\n{}\n+\n{}",String::from_utf8_lossy( self.id() ), self.as_dna_string(), String::from_utf8_lossy( self.qual() )  )
    }
}

impl PartialEq for SeqRec {
    fn eq(&self, other: &Self) -> bool {
        // Implement your equality comparison logic here
        self.id == other.id && self.seq == other.seq && self.qual == other.qual
    }
}

impl Add for SeqRec {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut id = self.id;
        id.extend(other.id);

        let mut seq = self.seq;
        seq.extend(other.seq);

        let mut qual = self.qual;
        qual.extend(other.qual);

        SeqRec { id, seq, qual }
    }
}

impl AddAssign for SeqRec {
    fn add_assign(&mut self, other: Self) {
        self.id.extend_from_slice(&other.id);
        self.seq.extend_from_slice(&other.seq);
        self.qual.extend_from_slice(&other.qual);
    }
}

impl Default for SeqRec {
    fn default() -> Self {
        SeqRec {
            id: Vec::new(),
            seq: Vec::new(),
            qual: Vec::new(),
        }
    }
}

impl SeqRec{
	pub fn new(id:&[u8], seq:&[u8], qual:&[u8])->Self{
		if  seq.len()!= qual.len() {
			panic!("I need to get three arrays of the same length! And I got arrays of different lenth here!:\n{id:?}, \n{seq:?},\n{qual:?}")
		}
		Self{
			id: id.to_vec(),
			seq: seq.to_vec(),
			qual: qual.to_vec(),
		}

	}

    /*pub fn get_nucleotide_2bit(&self, pos: usize) -> Option<u8> {
        if pos < self.len(){
            Some( self.enc( &self.seq[pos]) )
        }else{
            None
        }
    }*/
	pub fn id(&self) -> &[u8] {
        &self.id
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    pub fn qual(&self) -> &[u8] {
        &self.qual
    }

    pub fn len(&self) -> usize{
        self.seq.len()
    }

    pub fn slice(&self, start:usize, len:usize) -> Option<Self>{
        
        let real_start = if start > self.len(){
            0
        }else {
            start
        };	

    	let used_len = if real_start + len > self.seq.len(){
            #[cfg(debug_assertions)]
            eprintln!( "you try a length of {len} from start {real_start} but the sequence is only {}bp long\nI'll give you the max length.\n{self}",self.seq.len() );
    		self.seq.len() -start
    	}else{
            len
        }; 

        let end  = real_start + used_len;

        let id = String::from_utf8_lossy(self.id.as_slice()).to_string() + &format!("{}+{}",real_start, used_len);
    	Some( 
    		Self::new( &id.into_bytes(),
    		 &self.seq[real_start..end],
    		 &self.qual[real_start..end]
    		)
    	)
    }
    pub fn as_dna_string(&self) -> String{
    	String::from_utf8_lossy( self.seq() ).to_string()
    }

    // Convert the first 32 seq data to a u64
    pub fn to_u64(&self) -> u64 {
    	let mut ret = 0_u64;
    	//for byte in self.seq[0..32.min(self.seq.len())].iter().rev() {
        for byte in self.seq.iter().take( 32.min(self.seq.len()) ).rev() {
    		ret <<= 2;
	        ret |= self.enc::<u64>(byte);
	    }
        ret
    }

    // Convert the first 16 seq data to a u32
    pub fn to_u32(&self) -> u32 {
        let mut ret = 0_u32;
        //eprintln!("I am returning the value for {} nucleotides:" ,16.min(self.seq.len()) );
    	//for byte in self.seq[0..16.min(self.seq.len())].iter().rev() {
        for byte in self.seq.iter().take( 16.min(self.seq.len()) ).rev() {
    		ret <<= 2;
	        ret |= self.enc::<u32>(byte);
	    }
        ret
    }

    // Convert the first 8 seq data to a u16
    pub fn to_u16(&self) -> u16 {
        let mut ret = 0_u16;
    	//for byte in self.seq[0..8.min(self.seq.len())].iter().rev() {
        for byte in self.seq.iter().take( 8.min(self.seq.len()) ).rev() {
    		ret <<= 2;
	        ret |= self.enc::<u16>(byte);
	    }
        ret
    }

    fn enc<T: From<u8>>(&self, base: &u8) -> T {
        /*match LOOKUP_TABLE.get(&base) {
            Some(&encoding) => encoding,
            None => panic!("Nucleotide not found in lookup table"),
        }
        */
        
	    match *base {
	        b'A' | b'a' => T::from(0_u8),
	        b'C' | b'c' => T::from(1_u8),
	        b'G' | b'g' => T::from(2_u8),
	        b'T' | b't' => T::from(3_u8),
	        b'N' | b'n' => T::from(0_u8), // You might want to handle 'N' differently
	        _ => panic!("cannot encode base into 2-bit encoding"),
	    }
	}



}