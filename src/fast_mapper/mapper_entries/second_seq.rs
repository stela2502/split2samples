
use std::hash::{Hash, Hasher};
use core::fmt;

use crate::traits::{ BinaryMatcher, Cell, Direction};

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
        let mask: u64;
        let length:u8;
        if other.1 > self.1{
            length =  self.1
        }else {
            length = other.1;
        }
        if length < 20 {
            return false
        }

        if length >= 32 {
            mask = u64::MAX;
            // Rest of your code using the mask
        } else {
            mask = (1 << (length as u64) *2 ) - 1;
        }
        //eprintln!("I'll compare {:b} to {:b}", (self.0 & mask) , (other.0 & mask));
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

    fn to_string(&self) -> String {
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

        let mut a_sum = vec![0_i8;64];
        let mut b_sum = vec![0_i8;64];

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

        let mut a_sum = vec![0_i8;16];
        let mut b_sum = vec![0_i8;16];

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

    /// Almost a needleman_wunsch implementation. It just returns the difference from the expected result
    /// comparing the sequences in there minimal defined length. Similar to the hamming_distance function.
    /// for sequences shorter than 15 bp this fails and returns 100.0
    fn needleman_wunsch(&self, other: &SecondSeq ) -> f32 {

        let size = self.min_length(other).min(33);

        if size < 15  || self.tri_nuc_abs_diff(other) > 0.6{
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
        let mut a: u64;
        let mut b: u64;

        for i in 1..rows {
            for j in 1..cols {
                a = (self.0 >> ((i-1) * 2)) & 0b11;
                b = (other.0 >> ((j-1) * 2)) & 0b11;
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

        /*println!("Can that be cut short(di_diff {}, tri_diff {}, NW {}) : \n{self} vs \n{other}", 
            self.di_nuc_abs_diff(other),  
            self.tri_nuc_abs_diff(other),  
            (size as i32 - matrix[rows - 1][cols - 1].score).abs() as f32 / size as f32 );
        */
        (size as i32 - matrix[rows - 1][cols - 1].score).abs() as f32 / size as f32
    }

    /// calculate the base flips between two u64 sequences
    /// stops after having detected 4 different bases.
    fn hamming_distance(self, other: &SecondSeq) -> u32 {
        
        //let mask:u64;
        let size = self.min_length(other);

        let mut a: u64;
        let mut b: u64;
        let mut ret: u32 = 0;
        for i in 0..size{
            a = (self.0 >> (i * 2)) & 0b11;
            b = (other.0 >> (i * 2)) & 0b11;
            if a != b {
                ret +=1;
            }
            //ret += HAMMING_LOOKUP[(a ^ b) as usize];
            if ret == 4{
                 break;
            }
        }
        //println!("hamming dist was {ret}");
        return ret
        
        /*
        // quite much slower!
        let size = usize::min(self.1 as usize, other.1 as usize);

        let mut a_shifted = self.0;
        let mut b_shifted = other.0;
        
        let mut ret: u32 = 0;

        for _ in 0..size {
            let a_value = a_shifted & 0b11;
            let b_value = b_shifted & 0b11;
            if a_value != b_value {
                ret +=1;
            }
            a_shifted >>= 2;
            b_shifted >>= 2;
        }

        ret
        */
    }

    fn table(&self) -> std::collections::HashMap<char, u32> {
        let mut a_cnt = 0;
        let mut c_cnt = 0;
        let mut g_cnt = 0;
        let mut t_cnt = 0;
        let sequence = self.0;

        for i in 0..self.1 {
            let pair = (sequence >> (i * 2)) & 0b11;
            match pair {
                0b00 => a_cnt += 1,
                0b01 => c_cnt += 1,
                0b10 => g_cnt += 1,
                0b11 => t_cnt += 1,
                _ => {} // Handle invalid pairs if needed
            }
        }

        let mut counts = std::collections::HashMap::new();
        counts.insert('A', a_cnt);
        counts.insert('C', c_cnt);
        counts.insert('G', g_cnt);
        counts.insert('T', t_cnt);

        counts
    }


}


// Implementing Display trait for SecondSeq
impl fmt::Display for SecondSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SecondSeq (u64: {} or {:b}, u8: {})", BinaryMatcher::to_string(self), self.0, self.1)
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
        println!("u64: {:b} or {:?}", self.0, BinaryMatcher::to_string(self) );
        println!("{} sig bp", self.1);
    }

    pub fn fuzzy_match(&self, other:&SecondSeq, max_dist:f32 ) -> bool {

        //return self.hamming_distance( other ) <= max_dist.try_into().unwrap()
        self.needleman_wunsch( other ) <= max_dist.try_into().unwrap()
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

    #[test]
    fn test_second_table() {
        let seq1 = SecondSeq(0b1010101011010101, 32);
        let mut exp = std::collections::HashMap::new();
        exp.insert('A', 24);
        exp.insert('C', 3);
        exp.insert('G', 4);
        exp.insert('T', 1);

        assert_eq!(seq1.table() , exp,  "table did return the right counts");
    }

    #[test]
    fn test_second_seq_hashing() {
        let mut map = std::collections::HashMap::new();
        let seq1 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq2 = SecondSeq(0b1111000010101010101010101010101010101010, 20);

        map.insert(seq1, "Value for seq1");

        assert_eq!(map.get(&seq1), Some(&"Value for seq1")); // Test retrieval by equal key
        assert_eq!(map.get(&seq2), None); // Test retrieval by different key
    }

    #[test]
    fn test_humming2() {
        let seq1 = SecondSeq(0b101010, 20);
        let seq2 = SecondSeq(0b101010, 20);
        assert_eq!( seq1.hamming_distance( &seq2 ), 0 );
        let seq3 = SecondSeq(0b011010, 20);
        assert_eq!( seq1.hamming_distance( &seq3 ), 1 );
        let seq4 = SecondSeq(0b001010, 20);
        assert_eq!( seq1.hamming_distance( &seq4 ), 1 );
        let seq5 = SecondSeq(0b011001, 20);
        assert_eq!( seq1.hamming_distance( &seq5 ), 2 );
        let seq6 = SecondSeq(0b0, 20);
        assert_eq!( seq1.hamming_distance( &seq6 ), 3 );
    }

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
