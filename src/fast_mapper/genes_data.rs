//fast_mapper/genes_data.rs
use std::collections::HashMap;
use crate::traits::BinaryMatcher ;
use crate::traits::Direction;
use crate::traits::Cell;
use core::cmp::max;

// logics copied from https://github.com/COMBINE-lab/kmers/
pub type Base = u8;
pub const A: Base = 0;
pub const C: Base = 1;
pub const G: Base = 2;
pub const T: Base = 3;

const MATCH_SCORE: i32 = 1;
const MISMATCH_SCORE: i32 = -1;
const GAP_PENALTY: i32 = -2;

#[derive(Debug,PartialEq)]
pub struct GenesData {
	/// The data storage - we will never store any seq other than a 2bit encoded one.
	u8_encoded: Vec::<u8>,
	/// store the encoded length of the sequence (needed for subset of matching)
	length:usize,
	name: String,
	chr: String,
	start: usize,
}

impl GenesData{
	pub fn new( seq: &[u8], name:&str, chr:&str, start:usize )-> Self{
		let u8_encoded = Self::encode( seq );
		Self{
			u8_encoded,
			length: seq.len(),
			name: name.to_string(),
			chr: chr.to_string(),
			start,
		}
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


	pub fn from_bytes( seq: &[u8] ) -> Self{
		let u8_encoded = Self::encode( seq );
		Self{
			u8_encoded,
			length: seq.len(),
			name: "unnamed".to_string(),
			chr: "unknown".to_string(),
			start:0,
		}
	}

	pub fn encode(decoded: &[u8]) -> Vec<u8> {
	    let mut encoded = Vec::<u8>::with_capacity((decoded.len() + 3) / 4);
	    let mut count = 0;
	    let mut byte = 0;
	    for &base in decoded {
	        byte |= match base {
	            b'A' | b'a' => 0b00,
	            b'C' | b'c' => 0b01,
	            b'G' | b'g' => 0b10,
	            b'T' | b't' => 0b11,
	            _ => panic!("cannot encode {base} into 2-bit encoding"),
	        } << (6 - 2 * count); // Shift the nucleotide bits to the correct position in the byte
	        count += 1;
	        if count == 4 {
	            encoded.push(byte);
	            byte = 0;
	            count = 0;
	        }
	    }
	    if count != 0 {
	    //   encoded.push(byte << (6 - 2 * count)); // Pad the remaining bits and add the final byte
	    	encoded.push(byte)
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
	    let mut decoded = Vec::<u8>::with_capacity(encoded.len() * 4);
	    let mut enc_nuc = 0;
	    for &byte in encoded {
	        for i in 0..4 {
	            let two_bits = (byte >> (6- 2 * i)) & 0b11; // Extract 2 bits at a time
	            let nucleotide = match two_bits {
	                0b00 => b'A',
	                0b01 => b'C',
	                0b10 => b'G',
	                0b11 => b'T',
	                _ => panic!("Invalid 2-bit encoding detected"),
	            };

	            decoded.push(nucleotide);
	            // check if we have encoded enough!
	            enc_nuc += 1;
	            if enc_nuc == len {
	            	return Some(decoded) 
	            }
	        }
	    }

	    Some(decoded)
	}

	/// Slice the sequence of 2-bit representations
    pub fn slice(&self, start: usize, end: usize) -> Self {
        let mut sliced_data = Vec::with_capacity((end +3 - start) / 4 );
        let mut current_byte = 0;
        let mut bit_position = 0;

        for i in start..end {
            let byte_index = i / 4;
            let bit_index = i % 4;

            let value = (self.u8_encoded[byte_index] >> (6 - 2 * bit_index)) & 0b11;
            current_byte |= value << (6 - 2 * bit_position);

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

        Self {
            u8_encoded: sliced_data,
            length: end - start,
            name: (self.name.to_string() + " slice").to_string(),
            chr : self.chr.to_string(),
            start : self.start + start,
        }
    }

    pub fn get_nucleotide_2bit(&self, pos: usize) -> Option<u8> {
	    if pos > self.len() {
	        return None; // Position exceeds the length of the encoded sequence
	    }

	    let byte_idx = pos / 4;
	    let bit_offset = pos % 4;

	    Some((self.u8_encoded[byte_idx] >> (6 - bit_offset * 2)) & 0b11)
	}

	fn max3<T: Ord>(a: T, b: T, c: T) -> T {
        max(a, max(b, c))
    }

}

impl BinaryMatcher for GenesData{
	fn to_string(&self) -> String {
		if let Some(decoded) = self.to_bytes( self.len() ){
			String::from_utf8_lossy(&decoded).into_owned()
		}else {
			"".to_string()
		}
	}

	fn di_nuc_abs_diff( &self, other: &Self  ) -> f32 {
		if (self.len() > 100) | ( other.len() > 100) {
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
		if (self.len() > 100) | ( other.len() > 100) {
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
        	let off = order[i % 4];
        	if off < 3{
        		let acc = (self.u8_encoded[id] >> (off * 2)) & 0b1111;
        		sum[ acc as usize] +=1;
        	}else {
        		let acc = ((self.u8_encoded[id] & 0b11) << 2) | ((self.u8_encoded[id + 1] >> 6) & 0b11);
        		//let acc = (self.u8_encoded[id] & 0b11) | ((self.u8_encoded[id + 1] >> 6) & 0b11)  ;
        		//println!("-äI try to do a {id} {off} merge and get that 'index' {acc} or {acc:b} ");
        		sum[acc as usize] +=1;
        	}
        }
        sum
	}

	fn tri_nuc_tab(&self) -> Vec<i8> {
	    let mut sum = vec![0_i8; 64];
	    let order = vec![ 1, 0, 2, 3];

	    for i in 0..self.length - 2 {
	        let id = i / 4;
	        let off = order[i % 4];
	        if off < 2 {
	            let acc = (self.u8_encoded[id] >> (off * 2)) & 0b111111;
	            //println!("I try to do a {id} {off} merge and get that 'index' {acc} or {acc:b} ");
	            sum[acc as usize] += 1;
	        } else if off == 2 {
	            let acc = ((self.u8_encoded[id] & 0b1111) << 2) | ((self.u8_encoded[id + 1] >> 6) & 0b11);
	            //println!("-I try to do a {id} {off} merge and get that 'index' {acc} or {acc:b} ");
	            sum[acc as usize] += 1;
	        } else {
	            let acc = ((self.u8_encoded[id] & 0b11) << 4) | ((self.u8_encoded[id + 1] >> 4) & 0b1111);
	            //println!("--I try to do a {id} {off} merge and get that 'index' {acc} or {acc:b} ");
	            sum[acc as usize] += 1;
	        }
	    }
	    sum
	}

	fn needleman_wunsch(&self, other: &Self, humming_cut: f32 ) -> f32 { 
		let size = self.len().min(other.len());

        if size < 15  || self.tri_nuc_abs_diff(other) > humming_cut{
            return 100.0
        }

        let rows: usize = size;
        let cols: usize = size;

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
	            let a = self.get_nucleotide_2bit(i-1);
	            let b = other.get_nucleotide_2bit(j - 1);
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
	        }
	    }

        // Uncomment the following lines to print the alignment matrix
        /*for i in 0..rows {
            for j in 0..cols {
                print!("{:4} ", matrix[i][j].score);
            }
            println!();
        }*/

        /*
        println!("Can that be cut short(di_diff {}, tri_diff {}, NW {}) : \n{self} vs \n{other}", 
            self.di_nuc_abs_diff(other),  
            self.tri_nuc_abs_diff(other),  
            (size as i32 - matrix[rows - 1][cols - 1].score).abs() as f32 / size as f32 );
        */
        (size as i32 - matrix[rows - 1][cols - 1].score).abs() as f32 / size as f32
	}

}
