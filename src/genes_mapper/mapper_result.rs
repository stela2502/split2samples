// mapper_result.rs

use crate::genes_mapper::Cigar;
//use crate::genes_mapper::CigarEndFix;

use core::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct MapperResult{
	gene_id: usize,
	start: usize,
	save: bool,
	cigar: Option<Cigar>,
	mapq: u8,
	score: usize,
    nw: f32,
	edit_dist: f32,
	gene_name:String,
    /// the length of the database entry - is necessary to estimate if a short match is good enough (e.g. sampleid or AB tag)
    db_length: usize,

}

// Implementing Display trait for MapperResult
impl fmt::Display for MapperResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MapperResult ( gene_id {}, gene_name {}, start {}, save {}, cigar {:?}, nw {}, mapq {}, edit_dist {} )\n", 
            self.gene_id, self.gene_name, self.start, self.save, self.cigar, self.nw, self.mapq, self.edit_dist )
    }
}

impl Default for MapperResult {
    fn default() -> Self {
        MapperResult {
            gene_id: 0,
            start: 0,
            save: false,
            cigar: None,
            mapq: 0,
            score: 0,
            nw:0.0,
            edit_dist: 0.0,
            gene_name: String::default(),
            db_length: 0,
        }
    }
}


impl MapperResult{
	pub fn new( gene_id:usize,  start: usize, save:bool, cigar:Option<Cigar>, mapq:u8, nw:f32,score:usize, edit_dist:f32, name:&str, db_length:usize ) -> Self{

        Self{
			gene_id,
			start,
			save,
			cigar,
			mapq,
			score,
            nw,
			edit_dist,
			gene_name: name.to_string(),
            db_length,
		}
	}

    pub fn db_length(&self) -> usize{
        self.db_length
    }

    pub fn get_nw(&self) -> f32{
        self.nw
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

    pub fn cigar(&self) -> &Option<Cigar> {
        &self.cigar
    }

    pub fn get_cigar(&self) -> Cigar{
        match &self.cigar{
            Some(cig) => cig.clone(),
            None => Cigar::default(),
        }
    }

    pub fn mapq(&self) -> u8 {
        self.mapq
    }
    pub fn score(&self) -> usize {
        self.score
    }
    pub fn edit_dist(&self) -> f32 {
        self.edit_dist
    }
}