//fast_mapper/genes_data.rs
use std::collections::HashMap;
use crate::traits::{ BinaryMatcher, Direction, Cell};
use crate::genes_mapper::Cigar;
use core::cmp::max;
use std::hash::{Hash, Hasher};
use core::fmt;
use num_traits::cast::ToPrimitive;

use serde::{Serialize, Deserialize};

use std::fs::File;
use std::io::prelude::*;
use std::io;

// logics copied from https://github.com/COMBINE-lab/kmers/
pub type Base = u8;
pub const A: Base = 0;
pub const C: Base = 1;
pub const G: Base = 2;
pub const T: Base = 3;

const MATCH_SCORE: i32 = 1;
const MISMATCH_SCORE: i32 = -1;
const GAP_PENALTY: i32 = -2;

#[derive(Debug, Clone,Deserialize,Serialize)]
pub struct GeneData {
	/// The data storage - we will never store any seq other than a 2bit encoded one.
	u8_encoded: Vec::<u8>,
	/// store the encoded length of the sequence (needed for subset of matching)
	length:usize,
	/// the name for the sequence (mostly a gene name here)
	name: String,
	/// where came this gene from (chr) not used at the moment
	chr: String,
	/// alternative gene names (like transcript id or family name)
	//other_ids:Vec<String>,
	/// If this is a subseq of the main entry - where on the main entry does this start
	start: usize,
	/// and this is something provate used for the next function
	current_position: usize,
}

impl Hash for GeneData {
    fn hash<H: Hasher>(&self, state: &mut H) {
    	for val in &self.u8_encoded{
    		val.hash(state);
    	}
    }
}

impl PartialEq for GeneData {
    fn eq(&self, other: &Self) -> bool {
        (self.len() == other.len()) && (self.equal_entries( other) == self.len() / 4)
    }
}


impl Iterator for GeneData {
    type Item = (u16, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.key_at_position(self.current_position );
        self.current_position += 1;
        return result
    }
}

// Implementing Display trait for SecondSeq
impl fmt::Display for GeneData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GeneData name {} on chr {} at position {} with {} entries:\n{}",
        	self.name, self.chr, self.start, self.len(), self.as_dna_string() )
    }
}

impl GeneData{
	//pub fn new( seq: &[u8], name:&str, chr:&str, start:usize, other_ids:Vec<String> )-> Self{
	pub fn new( seq: &[u8], name:&str, chr:&str, start:usize ) -> Self {
		let u8_encoded = Self::encode( seq );
		Self{
			u8_encoded,
			length: seq.len(),
			name: name.to_string(),
			chr: chr.to_string(),
			start,
			//other_ids: other_ids,
			current_position:0,
		}
	}

	/// returns the GeneData's data as fasta sequence - >name|chr|start\nseq\n
	pub fn to_fasta(&self) -> String{
		let mut fasta = format!(">{}|{}|{}\n", self.name, self.chr, self.start);
		fasta += &self.as_dna_string();
		fasta += "\n";
		fasta
	}

	// Convert the u8_encoded data to a u64
    pub fn to_u64(&self) -> u64 {
    	let mut ret = 0_u64;
    	for byte in self.u8_encoded[0..8].iter().rev() {
    		ret <<= 8;
	        ret |= *byte as u64;
	    }
        ret
    }

    // Convert the u8_encoded data to a u32
    pub fn to_u32(&self) -> u32 {
        let mut ret = 0_u32;
    	for byte in self.u8_encoded[0..4].iter().rev() {
    		ret <<= 8;
	        ret |= *byte as u32;
	    }
        ret
    }

    // Convert the u8_encoded data to a u16
    pub fn to_u16(&self) -> u16 {
        let mut ret = 0_u16;
    	for byte in self.u8_encoded[0..2].iter().rev() {
    		ret <<= 8;
	        ret |= *byte as u16;
	    }
        ret
    }

    

	// get the number of equal 4bp stretches over the whole length of the element
	pub fn equal_entries(&self, other: &Self ) -> usize{
		let mut ret  = 0;
		for (item1, item2) in self.u8_encoded.iter().zip(other.u8_encoded.iter()) {
			if item1 == item2 {
				ret += 1;
			}
		}
		ret
	}

	pub fn len( &self ) -> usize{
		self.length
	}

	pub fn set_name(&mut self, name:&str) {
		self.name = name.to_string();
	}
	pub fn get_name(&self) -> &str{
		&self.name
	}
	// Getter method for chromosome
    pub fn get_chr(&self) -> &String {
        &self.chr
    }

    // Setter method for chromosome
    pub fn set_chr(&mut self, chr: String) {
        self.chr = chr;
    }

    // Getter method for start position
    pub fn get_start(&self) -> usize {
        self.start
    }

    // Setter method for start position
    pub fn set_start(&mut self, start: usize) {
        self.start = start;
    }

    /// the start here is relative to this DNA fragment, not relative to the source sequence
    /// that the objects start points to.
    pub fn key_at( &self, start:usize) -> u16 {
    	let from = start / 4;
    	let to = (start +3 +8) / 4;
    	let mut u32_value = 0u32;
	    for byte in self.u8_encoded[from..to].iter().rev() {
	        u32_value <<= 8;
	        u32_value |= *byte as u32;
	    }
	    let off = start % 4;
	    (u32_value >> off*2) as u16
    }


	pub fn from_bytes( seq: &[u8] ) -> Self{
		let u8_encoded = Self::encode( seq );
		Self{
			u8_encoded,
			length: seq.len(),
			name: "unnamed".to_string(),
			chr: "unknown".to_string(),
			start:0,
			current_position:0,
		}
	}

	pub fn encode_binary(  base:u8) -> Base {
		match base {
	        b'A' | b'a' => A,
	        b'C' | b'c' => C,
	        b'G' | b'g' => G,
	        b'T' | b't' => T,
	        b'N' | b'n' => A,
	        _ => panic!("cannot encode {base} into 2-bit encoding"),
	    }
	}

	pub fn encode(decoded: &[u8]) -> Vec<u8> {

		let mut encoded = Vec::<u8>::with_capacity((decoded.len() + 3) / 4);

	    for id in 0..((decoded.len()+3) / 4) {
	        encoded.push(0_u8);
	        for add in (0..4).rev() {
	            encoded[id] <<= 2;
	            let current_byte_id = id * 4 + add;
	            if current_byte_id < decoded.len() {
	                encoded[id] |= Self::encode_binary(decoded[current_byte_id]);
	            } else {
	                encoded[id] |= Self::encode_binary(b'A'); // If sequence is shorter, pad with 'A'
	            }
	        }
	    }

	    encoded
	}

	pub fn to_bytes( &self, length:usize ) -> Option<Vec<u8>> {
		Self::decode( &self.u8_encoded, length )
	}

	/// This will return none if there is not enough data in the encoded to fill the requirement!
	pub fn decode(encoded: &[u8], len:usize) -> Option<Vec<u8>> {
		if len > encoded.len() *4{
			return None
		}
	    let mut decoded = Vec::<u8>::with_capacity(len);
	    let mut enc_nuc = 0;

	    for &byte in encoded {
	        for add in 0..4 {
	            let two_bits = (byte >> (add * 2)) & 0b11; // Extract 2 bits at a time, starting from the least significant bits
	            let nucleotide = match two_bits {
	                0b00 => b'A',
	                0b01 => b'C',
	                0b10 => b'G',
	                0b11 => b'T',
	                _ => panic!("Invalid 2-bit encoding detected"),
	            };

	            decoded.push(nucleotide);
	            enc_nuc += 1;

	            if enc_nuc == len {
	                return Some(decoded);
	            }
	        }
	    }

	    Some(decoded)
	}

	/// Slice the sequence of 2-bit representations
    pub fn slice(&self, start: usize, length: usize) -> Option<Self> {
    	if start >= self.len(){
    		return None
    	}
    	let mut end =start + length;

    	if end > self.len(){
    		//eprintln!("You have requested more than I have ({end} > {} ) - giving you only {}! {self}", self.len(), self.len());
    		end = self.len();
    	}
        let mut sliced_data = Vec::with_capacity( (length +3) / 4 );
        let mut current_byte = 0;
        let mut bit_position = 0;

        for i in start..end {
            let byte_index = i / 4;
            let bit_index = i % 4;
            if byte_index >= self.u8_encoded.len(){
            	panic!("WHAAAT? subset {start} to {end} on {self} leeds to an out of range {byte_index}");
            }

            let value = (self.u8_encoded[byte_index] >> ( 2 * bit_index)) & 0b11;
            //println!("I got {value:b} out of the array at self.u8_encoded[{byte_index}] and bit_index {bit_index}");
            current_byte |= value <<  (2 * bit_position); // Corrected bit shifting
            //println!("and my current value looks now like that {current_byte:b} ");

            if bit_position == 3 {
                sliced_data.push(current_byte);
                current_byte = 0;
                bit_position = 0;
            } else {
                bit_position += 1;
            }
        }

        if bit_position != 0 {
            sliced_data.push(current_byte);
        }

        Some(Self {
            u8_encoded: sliced_data,
            length: end - start,
            name: (self.name.to_string() + " slice").to_string(),
            chr : self.chr.to_string(),
            start : self.start + start,
            current_position:0,
        })
    }

    /// This function returns the 8bp following the start (0 indexed).
    pub fn key_at_position(&self, start:usize ) -> Option<(u16, usize)> {

    	// Ensure the start position is within bounds
	    if start + 8 > self.len() {
	        return None;
	    }
    	
    	// Determine the range of bytes to consider
    	let start_id = start / 4;
    	let shift_amount = (start % 4) * 2;
    	let mut u32_value = 0_u32;
    	let to = if shift_amount > 0 {
    		3
    	}else {
    		2
    	};
    	for byte in self.u8_encoded[start_id..(start_id +to)].iter().rev() {
    		//println!("Bytes before <<=8 {u32_value:b}");
	        u32_value <<= 8;
	        //println!("Bytes before |= {u32_value:b}");
	        u32_value |= *byte as u32;
	        //println!("Bytes after both steps {u32_value:b}");
	    }

	    let shifted_value = u32_value >> shift_amount;
	    let u16_value = (shifted_value & u16::MAX as u32) as u16;
	    //println!("next returns a u16 value: {u16_value:b}");
	    Some( (u16_value, start) )
    } 
}

impl BinaryMatcher for GeneData{

	fn max3<T: Ord>(a: T, b: T, c: T) -> T {
        max(a, max(b, c))
    }

	fn get_nucleotide_2bit(&self, pos: usize) -> Option<u8> {
	    if pos >= self.len() {
	        return None; // Position exceeds the length of the encoded sequence
	    }

	    let byte_idx = pos / 4;
	    let bit_offset = pos % 4;
	    if byte_idx >= self.u8_encoded.len(){
	    	return None
	    }
	    //println!("Byte idx {byte_idx}, bit_offset {bit_offset}");
	    Some((self.u8_encoded[byte_idx] >> ( bit_offset * 2)) & 0b11)
	}



	fn as_dna_string(&self) -> String {
		if let Some(decoded) = self.to_bytes( self.len() ){
			String::from_utf8_lossy(&decoded).into_owned()
		}else {
			"".to_string()
		}
	}

	fn di_nuc_abs_diff( &self, other: &Self  ) -> f32 {
		if (self.len() > 500) | ( other.len() > 500) {
			panic!("Please do not do that with large sequences!")
		}

		let a_sum = self.di_nuc_tab();
		let b_sum = other.di_nuc_tab();
		let mut ret = 0;
		for i in 0..a_sum.len(){
            ret += a_sum[i].abs_diff(b_sum[i])
        }
		ret as f32 / 2.0 / self.len().max(other.len()) as f32
	}

	fn tri_nuc_abs_diff( &self, other: &Self  ) -> f32 { 
		if (self.len() > 500) | ( other.len() > 500) {
			panic!("Please do not do that with large sequences!")
		}

		let a_sum = self.tri_nuc_tab();
		let b_sum = other.tri_nuc_tab();
		let mut ret = 0;
		for i in 0..a_sum.len(){
            ret += a_sum[i].abs_diff(b_sum[i])
        }
		ret as f32 / 2.0 / self.len().max(other.len()) as f32
	}

	fn di_nuc_tab (&self ) -> Vec<i8> {
		let mut sum = vec![0_i8;16];
		let order= vec![2,1,0,3];
        for i in 0..self.length-1 {
        	let id = i / 4;
        	let off =  i % 4;
        	//let off = order[i % 4];
        	if off < 3{
        		let acc = (self.u8_encoded[id] >> (off * 2)) & 0b1111;
        		//println!("I try to do a {id} {off} merge and get that 'index' {acc} or {acc:b} ");
        		sum[ acc as usize] +=1;
        	}else {
        		//let acc = ((self.u8_encoded[id] & 0b11) << 2) | ((self.u8_encoded[id + 1] >> 6) & 0b11);
        		//let acc = (((self.u8_encoded[id] >> (off * 2)) & 0b11)<<2) | (self.u8_encoded[id + 1] & 0b11) ;
        		let acc = ((self.u8_encoded[id] >> 6 & 0b11) <<2) | (self.u8_encoded[id + 1] & 0b11) ;
        		//println!("-äI try to do a {id} {off} merge and get that 'index' {acc} or {acc:b} ");
        		sum[acc as usize] +=1;
        	}
        }
        sum
	}

	fn tri_nuc_tab(&self) -> Vec<i8> {
	    let mut sum = vec![0_i8; 64];
	    //let order = vec![ 1, 0, 2, 3];

	    for i in 0..self.length - 2 {
	        let id = i / 4;
	        //let off = order[i % 4];
	        let off = i % 4;
	        if off < 2 {
	            let acc = (self.u8_encoded[id] >> (off * 2)) & 0b111111;
	            //println!("I try to do a {id} {off} merge and get that 'index' {acc} or {acc:b} ");
	            sum[acc as usize] += 1;
	        } else if off == 2 {
	        	let acc = (self.u8_encoded[id] >> 4 & 0b1111) | ((self.u8_encoded[id + 1] & 0b11) << 4);
	            //println!("-I try to do a {id} {off} merge and get that 'index' {acc} or {acc:b} ");
	            sum[acc as usize] += 1;
	        } else {
	        	let acc = (self.u8_encoded[id] >> 6 & 0b11) | ((self.u8_encoded[id + 1] & 0b1111) << 2) ;
	            //println!("--I try to do a {id} {off} merge and get that 'index' {acc} or {acc:b} ");
	            sum[acc as usize] += 1;
	        }
	    }
	    sum
	}


	fn needleman_wunsch(&self, other: &Self, humming_cut: f32, cigar: Option<&mut Cigar> ) -> f32 { 
		let size = self.len().min(other.len());

        if size < 15  || self.tri_nuc_abs_diff(other) > humming_cut{
            return 100.0
        }

        let rows: usize = self.len();
        let cols: usize = other.len();

        let mut matrix = vec![vec![Cell { score: 0, direction: Direction::Diagonal }; cols]; rows];

        // Initialize the first row and column
        for i in 1..rows {
            matrix[i][0].score = matrix[i - 1][0].score + GAP_PENALTY;
            matrix[i][0].direction = Direction::Up;
        }

        for j in 1..cols {
            matrix[0][j].score = matrix[0][j - 1].score + GAP_PENALTY;
            matrix[0][j].direction = Direction::Left;
        }

        // Fill in the matrix
        for i in 1..rows {
	        for j in 1..cols {
	            if let ( Some(a), Some(b) ) = ( self.get_nucleotide_2bit(i-1), other.get_nucleotide_2bit(j-1) ) {

                    let match_score = if a == b { MATCH_SCORE } else { MISMATCH_SCORE };

                    let diagonal_score = matrix[i - 1][j - 1].score + match_score;
                    let up_score = matrix[i - 1][j].score + GAP_PENALTY;
                    let left_score = matrix[i][j - 1].score + GAP_PENALTY;

                    let max_score = Self::max3(diagonal_score, up_score, left_score);

                    matrix[i][j].score = max_score;

                    matrix[i][j].direction = match max_score {
                        _ if max_score == diagonal_score => Direction::Diagonal,
                        _ if max_score == up_score => Direction::Up,
                        _ if max_score == left_score => Direction::Left,
                        _ => unreachable!(),
                    };
                }else {
                    panic!("positions not reachable: {i}/{j} ({:?}/{:?})",  self.get_nucleotide_2bit(i-1), self.get_nucleotide_2bit(j-1));
                }
	        }
	    }

	    if let Some(cig) = cigar {
	    	cig.calculate_cigar( &matrix ,self.get_nucleotide_2bit(rows-1) == other.get_nucleotide_2bit(cols-1));
	    	cig.clean_up_cigar( self, other);
	    	/*if cig.calculate_covered_nucleotides( &cig.cigar ).0 != self.len(){
	    		panic!("This cigar does not describe my mapping correctly! {} has a length of {} and this is me:\n{self}\nand that the other:\n{other}", 
	    			&cig.cigar, cig.calculate_covered_nucleotides( &cig.cigar ).0);
	    	}*/
	    }

        (size as i32 - matrix[rows - 1][cols - 1].score).abs() as f32 / size as f32
	}

}
