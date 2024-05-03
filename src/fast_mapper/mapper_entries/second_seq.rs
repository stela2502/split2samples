
use std::hash::{Hash, Hasher};
use core::fmt;

use crate::traits::{ BinaryMatcher, Cell, Direction };
use crate::genes_mapper::Cigar;
use std::cmp::{max, Ord};

const MATCH_SCORE: i32 = 1;
const MISMATCH_SCORE: i32 = -1;
const GAP_PENALTY: i32 = -2;

/*static MASK2: u64 = 0xAAAAAAAAAAAAAAAA;
static MASK3: u64 = 0x5555555555555555;
const HAMMING_LOOKUP: [u32; 4] = [0, 1, 1, 1];*/


#[derive(Debug, Copy, Clone)]
pub struct SecondSeq( pub u64, pub u8);

impl PartialEq for SecondSeq {
    fn eq(&self, other: &Self) -> bool {
        
        let length:u8 = if other.1 > self.1{
            self.1
        }else {
            other.1
        };

        if length < 20 {
            return false
        }

        let mask: u64 = if length >= 32 {
            u64::MAX
            // Rest of your code using the mask
        } else {
            (1 << ((length as u64) *2) ) - 1
        };

        //eprintln!("I'll compare \n{:b} to \n{:b}", (self.0 & mask) , (other.0 & mask));
        (self.0 & mask) == (other.0 & mask)
    }
}

impl Eq for SecondSeq {}

impl Hash for SecondSeq {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}



impl BinaryMatcher for SecondSeq {

    fn len(&self) -> usize{
        self.1 as usize
    }

    fn max3<T: Ord>(a: T, b: T, c: T) -> T {
        max(a, max(b, c))
    }

    fn as_dna_string(&self) -> String {
        let mut data = String::new();
        //println!("converting u64 {loc:b} to string with {kmer_size} bp.");
        for i in 0..self.1.min(32) {
            // Use a mask (0b11) to extract the least significant 2 bits
            let ch = match (self.0 >> (i * 2)) & 0b11 {
                0b00 => "A",
                0b01 => "C",
                0b10 => "G",
                0b11 => "T",
                _ => "N",
            };
            data += ch;
            
            //println!("{ch} and loc {loc:b}");
        }

        //println!("\nMakes sense? {:?}", data);
        data
    }

    fn tri_nuc_tab (&self ) -> Vec<i8>{
        let mut sum = vec![0_i8;64];

        for i in 0..self.1.min(29){
            let a = (self.0 >> (i * 2)) & 0b111111;
            sum[a as usize] +=1;
        }
        sum
    }
    /// dnt will be a Vec usize, a table where the position of the vec represents the di-nucleotides from #b0000 to #b1111,
    /// and the value the number of occurances in that u64.
    fn di_nuc_tab (&self ) -> Vec<i8>{
        let mut sum = vec![0_i8;16];

        for i in 0..self.1.min(30){
            let a = (self.0 >> (i * 2)) & 0b1111;
            sum[a as usize] +=1;
        }
        sum
    }

    fn tri_nuc_abs_diff( &self, other: &SecondSeq  ) -> f32 {
        let size = self.min_length( other );

        let mut a_sum = [0_i8; 64];
        let mut b_sum = [0_i8; 64];

        for i in 0..size.min(29){
            let a = (self.0 >> (i * 2)) & 0b111111;
            let b = (other.0 >> (i * 2)) & 0b111111;
            a_sum[a as usize] +=1;
            b_sum[b as usize] +=1;
        }
        let mut ret = 0;
        for i in 0..64{
            ret += a_sum[i].abs_diff(b_sum[i]) as usize
        }
        ret as f32 / 2.0 / size as f32
    }

    fn di_nuc_abs_diff( &self, other: &SecondSeq  ) -> f32 {
        let size = self.min_length( other );

        let mut a_sum = [0_i8; 16];
        let mut b_sum = [0_i8; 16];

        for i in 0..size.min(30){
            let a = (self.0 >> (i * 2)) & 0b1111;
            let b = (other.0 >> (i * 2)) & 0b1111;
            a_sum[a as usize] +=1;
            b_sum[b as usize] +=1;
        }
        let mut ret = 0;
        for i in 0..16{
            ret += a_sum[i].abs_diff(b_sum[i]) as usize
        }
        ret as f32 / 2.0 / size as f32
    }

    fn get_nucleotide_2bit(&self, pos: usize) -> Option<u8> {
        if pos > self.1 as usize {
            return None; // Position exceeds the length of the encoded sequence
        }
        //println!("Byte idx {byte_idx}, bit_offset {bit_offset}");
        Some(((self.0 >> (pos * 2)) & 0b11) as u8)
    }

    /// Almost a needleman_wunsch implementation. It just returns the difference from the expected result
    /// comparing the sequences in there minimal defined length. Similar to the hamming_distance function.
    /// for sequences shorter than 15 bp this fails and returns 100.0
    fn needleman_wunsch(&self, other: &SecondSeq, humming_cut:f32, cigar: Option<&mut Cigar> ) -> f32 {

        let size = self.min_length(other).min(33);

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
                if let (Some(a), Some(b)) = ( self.get_nucleotide_2bit(j-1), other.get_nucleotide_2bit(i-1)){
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
                    panic!("Sequence not defined {i}/{j} self: {self}\nother: {other}");
                }
                
            }
        }

        // Trace back the alignment path
        if let Some( cig ) = cigar{
            cig.calculate_cigar( &matrix, self.get_nucleotide_2bit(rows-1) == other.get_nucleotide_2bit(cols-1) );
            cig.clean_up_cigar( self, other);
        };

        (size as i32 - matrix[rows - 1][cols - 1].score).abs() as f32 / size as f32
    }
}


// Implementing Display trait for SecondSeq
impl fmt::Display for SecondSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SecondSeq (u64: {} or {:b}, u8: {})", BinaryMatcher::as_dna_string(self), self.0, self.1)
    }
}

impl SecondSeq {

    const BYTE_SIZE: usize = std::mem::size_of::<u64>() + std::mem::size_of::<u8>();

    fn max3<T: Ord>(a: T, b: T, c: T) -> T {
        max(a, max(b, c))
    }

    pub fn to_le_bytes(&self) -> [u8; Self::BYTE_SIZE ] {
        let mut bytes = [0; Self::BYTE_SIZE ];
        bytes[..std::mem::size_of::<u64>()].copy_from_slice(&self.0.to_le_bytes());
        bytes[std::mem::size_of::<u64>()..].copy_from_slice(&self.1.to_le_bytes());
        bytes
    }

    pub fn from_le_bytes(bytes: [u8; Self::BYTE_SIZE]) -> Option<Self> {
        if bytes.len() == Self::BYTE_SIZE {
            let mut u64_bytes = [0; std::mem::size_of::<u64>()];
            let mut u8_bytes = [0; std::mem::size_of::<u8>()];

            u64_bytes.copy_from_slice(&bytes[..std::mem::size_of::<u64>()]);
            u8_bytes.copy_from_slice(&bytes[std::mem::size_of::<u64>()..]);

            let u64_val = u64::from_le_bytes(u64_bytes);
            let u8_val = u8::from_le_bytes(u8_bytes);

            Some(SecondSeq(u64_val, u8_val))
        } else {
            None
        }
    }
    

    /// returns the minimal size the two SecondSeq obejcts are defined for
    /// Something between 32 and 0 bp.
    pub fn min_length( &self, other: &SecondSeq) -> usize {

        if other.1 > self.1{
            self.1 as usize
            //mask = (1 << (self.1 as u64) *2 ) - 1;
        }else {
            other.1 as usize
            //mask = (1 << (other.1 as u64) *2 ) - 1;
        }
    }

    

    /// the == takes a mininmal matching region into a account
    /// This same function soes not do that.
    pub fn same(&self, other: &Self ) -> bool {
        self.0 == other.0
    }

    pub fn print_second_seq(&self) {
        println!("Contents of SecondSeq:");
        println!("u64: {:b} or {:?}", self.0, BinaryMatcher::as_dna_string(self) );
        println!("{} sig bp", self.1);
    }

    pub fn fuzzy_match(&self, other:&SecondSeq, max_dist:f32 ) -> bool {

        //return self.hamming_distance( other ) <= max_dist.try_into().unwrap()
        self.needleman_wunsch( other, max_dist, None ) <= max_dist
    }

    
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_second_seq_equality() {
        let seq1 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq2 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq3 = SecondSeq(0b1111000010101010101010101010101010101010, 20);
        let seq4 = SecondSeq(0b101010101010101010101010101010101010101010101010, 24);

        assert_eq!(seq1, seq2, "same sequences"); // Test equal sequences
        assert_ne!(seq1, seq3, "different sequences"); // Test unequal sequences
        assert_eq!(seq1, seq4,  "same, but one is longer");
    }

    // #[test]
    // fn test_second_table() {
    //     let seq1 = SecondSeq(0b1010101011010101, 32);
    //     let mut exp = std::collections::HashMap::new();
    //     exp.insert('A', 24);
    //     exp.insert('C', 3);
    //     exp.insert('G', 4);
    //     exp.insert('T', 1);

    //     assert_eq!(seq1.table() , exp,  "table did return the right counts");
    // }

    #[test]
    fn test_second_seq_hashing() {
        let mut map = std::collections::HashMap::new();
        let seq1 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq2 = SecondSeq(0b1111000010101010101010101010101010101010, 20);

        map.insert(seq1, "Value for seq1");

        assert_eq!(map.get(&seq1), Some(&"Value for seq1")); // Test retrieval by equal key
        assert_eq!(map.get(&seq2), None); // Test retrieval by different key
    }

    // #[test]
    // fn test_humming2() {
    //     let seq1 = SecondSeq(0b101010, 20);
    //     let seq2 = SecondSeq(0b101010, 20);
    //     assert_eq!( seq1.hamming_distance( &seq2 ), 0 );
    //     let seq3 = SecondSeq(0b011010, 20);
    //     assert_eq!( seq1.hamming_distance( &seq3 ), 1 );
    //     let seq4 = SecondSeq(0b001010, 20);
    //     assert_eq!( seq1.hamming_distance( &seq4 ), 1 );
    //     let seq5 = SecondSeq(0b011001, 20);
    //     assert_eq!( seq1.hamming_distance( &seq5 ), 2 );
    //     let seq6 = SecondSeq(0b0, 20);
    //     assert_eq!( seq1.hamming_distance( &seq6 ), 3 );
    // }

    #[test]
    fn test_di_nuc_abs_diff() {
        let seq1 = SecondSeq(0b101010, 20);
        let seq2 = SecondSeq(0b101010, 20);
        assert_eq!( seq1.di_nuc_abs_diff( &seq2 ), 0.0 );
        let seq3 = SecondSeq(0b011001, 20);
        assert_eq!( seq1.di_nuc_abs_diff( &seq3 ), 3.0/20.0 );
        assert_eq!( seq3.di_nuc_abs_diff( &seq1 ), 3.0/20.0 );
        let seq3 = SecondSeq(0b011001, 15);
        assert_eq!( seq3.di_nuc_abs_diff( &seq1 ), 3.0/15.0 );
    }
}
