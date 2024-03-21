use crate::traits::{BinaryMatcher, Cell, Direction};
use crate::genes_mapper::Cigar;
use core::cmp::max;

use std::hash::{Hash, Hasher};

use core::fmt;

/// This is in fact only a u32, but I'll attach all ma matching function to that

#[derive(Debug, Clone, Copy)]
pub struct CellId10x(pub u32);


const MATCH_SCORE: i32 = 1;
const MISMATCH_SCORE: i32 = -1;
const GAP_PENALTY: i32 = -2;


impl PartialEq for CellId10x {
    fn eq(&self, other: &Self) -> bool {
    	self.0 == other.0
    }
}

impl Eq for CellId10x {}

impl Hash for CellId10x {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

// Implementing Display trait for CellId10x
impl fmt::Display for CellId10x {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}",self.as_dna_string() )
    }
}


impl CellId10x {

	const BYTE_SIZE: usize = std::mem::size_of::<u32>();

    pub fn to_le_bytes(&self) -> [u8; Self::BYTE_SIZE ] {
    	self.0.to_le_bytes()
    }

    pub fn from_le_bytes( bytes: [u8; Self::BYTE_SIZE ] ) -> Option<Self> {
         Some(CellId10x(u32::from_le_bytes(bytes)))
    }

    
}

impl BinaryMatcher for CellId10x {

    fn max3<T: Ord>(a: T, b: T, c: T) -> T {
        max(a, max(b, c))
    }

    fn get_nucleotide_2bit(&self, pos: usize) -> Option<u8> {
        if pos > 16 {
            return None; // Position exceeds the length of the encoded sequence
        }
        //println!("Byte idx {byte_idx}, bit_offset {bit_offset}");
        Some(((self.0 >> ( pos * 2)) & 0b11) as u8)
    }

	fn as_dna_string(&self) -> String {
        let mut data = String::new();
        //println!("converting u64 {loc:b} to string with {kmer_size} bp.");
        for i in 0..16 {
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

        for i in 0..14{
            let a = (self.0 >> (i * 2)) & 0b111111;
            sum[a as usize] +=1;
        }
        sum
    }

    fn di_nuc_tab (&self ) -> Vec<i8>{
        let mut sum = vec![0_i8;16];

        for i in 0..15{
            let a = (self.0 >> (i * 2)) & 0b1111;
            sum[a as usize] +=1;
        }
        sum
    }

    fn tri_nuc_abs_diff( &self, other: &Self  ) -> f32 {
 
        let mut a_sum = [0_i8; 64];
        let mut b_sum = [0_i8; 64];

        for i in 0..14{
            let a = (self.0 >> (i * 2)) & 0b111111;
            let b = (other.0 >> (i * 2)) & 0b111111;
            a_sum[a as usize] +=1;
            b_sum[b as usize] +=1;
        }
        let mut ret = 0;
        for i in 0..64{
            ret += a_sum[i].abs_diff(b_sum[i]) as usize
        }
        ret as f32 / 2.0 / 16.0
    }


    fn di_nuc_abs_diff( &self, other: &Self  ) -> f32 {

        let mut a_sum = [0_i8; 16];
        let mut b_sum = [0_i8; 16];

        for i in 0..15{
            let a = (self.0 >> (i * 2)) & 0b1111;
            let b = (other.0 >> (i * 2)) & 0b1111;
            a_sum[a as usize] +=1;
            b_sum[b as usize] +=1;
        }
        let mut ret = 0;
        for i in 0..16{
            ret += a_sum[i].abs_diff(b_sum[i]) as usize
        }
        ret as f32 / 2.0 / 16.0
    }

    
    /// Almost a needleman_wunsch implementation. It just returns the difference from the expected result
    /// comparing the sequences in there minimal defined length. Similar to the hamming_distance function.
    /// for sequences shorter than 15 bp this fails and returns 100.0
    fn needleman_wunsch(&self, other: &Self, _humming_cut:f32, cigar: Option<&mut Cigar> ) -> f32 {

        let rows: usize = 16;
        let cols: usize = 16;

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
                    panic!("positions not reachable: {i}/{j}");
                }
            }
        }

        // Trace back the alignment path
        if let Some(cig) = cigar {
            cig.calculate_cigar( &matrix );
            cig.clean_up_cigar( self, other);
        }

        (16.0 - (matrix[rows - 1][cols - 1].score).abs() as f32) / 16.0
    }
    
}