/// here I plann to define all my traits - or the one I have planned for now :-D
use crate::errors::CellIdError;
use crate::genes_mapper::Cigar;
use crate::genes_mapper::SeqRec;

use core::fmt;


pub trait Index : Sync{

	fn get(&self, seq: &[u8]) ->  Option< usize >;
	fn add(&mut self, seq: &[u8],unique_name: &str, name: &str, class_ids: Vec<String> ) -> usize;
	fn get_id( &self, name: String ) -> usize;
	fn print( &self );
	fn change_start_id ( &mut self, new_start :usize );
	fn names_len(&self) -> usize;
	fn names(&self) -> Vec<String>;
	fn to_header_n( &self, names: &[String] ) -> std::string::String;
	fn max_id( &self ) -> usize; // return the max:id foir the sparse export of the data
	fn write_index( &mut self, path: String ) -> Result< (), &str>;
	fn load_index( &mut self, path: String ) -> Result< (), &str>;
}

pub trait CellIndex: Sync{
	//                                    former u32   , u64
	fn to_cellid (&self, r1: &SeqRec  )-> Result<( u32, u64, SeqRec, SeqRec ), CellIdError>;
}


/// likely necessary/helpful for the BinaryMatcher::needleman_wunsch
#[derive(Clone, Copy)]
pub struct Cell {
    pub score: i32,
    pub direction: Direction,
}

#[derive(Clone, Copy, PartialEq)]
pub enum Direction {
    Diagonal,
    Up,
    Left,
}

impl fmt::Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let direction_str = match self {
            Direction::Diagonal => "Diagonal",
            Direction::Up => "Up",
            Direction::Left => "Left",
        };
        write!(f, "{}", direction_str)
    }
}



pub trait BinaryMatcher : Sync + std::fmt::Display {
	fn max3<T: Ord>(a: T, b: T, c: T) -> T;
	fn get_nucleotide_2bit(&self, pos: usize) -> Option<u8>;
	fn as_dna_string(&self) -> String ;
	fn di_nuc_abs_diff( &self, other: &Self  ) -> f32;
	fn tri_nuc_abs_diff( &self, other: &Self  ) -> f32;
	fn di_nuc_tab (&self ) -> Vec<i8>;
	fn tri_nuc_tab (&self ) -> Vec<i8>;
	fn needleman_wunsch(&self, other: &Self, humming_cut: f32, cigar: Option<&mut Cigar> ) -> f32;
	fn len(&self) -> usize;
}
