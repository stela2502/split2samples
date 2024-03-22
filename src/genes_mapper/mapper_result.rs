// mapper_result.rs

use crate::genes_mapper::Cigar;

use core::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct MapperResult{
	gene_id: usize,
	start: usize,
	save: bool,
	cigar: Option<String>,
	mapq: u8,
	score: usize,
	edit_dist: usize,
	gene_name:String,
}

// Implementing Display trait for MapperResult
impl fmt::Display for MapperResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MapperResult ( gene_id {}, start {}, save {}, cigar {:?} )", self.gene_id, self.start, self.save, self.cigar )
    }
}

impl MapperResult{
	pub fn new( gene_id:usize,  start: usize, save:bool, cigar:Option<String>, mapq:u8, score:usize, edit_dist:usize, name:&str ) -> Self{
		Self{
			gene_id,
			start,
			save,
			cigar,
			mapq,
			score,
			edit_dist,
			gene_name: name.to_string(),
		}
	}

	pub fn get_name(&self) -> &str{
		&self.gene_name
	}

	pub fn gene_id(&self) -> usize {
        self.gene_id
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn save(&self) -> bool {
        self.save
    }

    pub fn cigar(&self) -> &Option<String> {
        &self.cigar
    }

    pub fn mapq(&self) -> u8 {
        self.mapq
    }
    pub fn score(&self) -> usize {
        self.score
    }
    pub fn edit_dist(&self) -> usize {
        self.edit_dist
    }
}