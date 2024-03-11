// genes_mapper/gene_link.rs
use std::collections::HashMap;


#[derive(PartialEq, Debug)]
pub struct GeneLink {
	// store all chr position combos for this u16 8bp area.
	data: Vec< (usize, usize)>,
}

impl GeneLink{
	pub fn new()->Self {
		Self{
			data: Vec::with_capacity(10),
		}
	}

	pub fn is_empty(&self) -> bool{
		self.data.is_empty()
	}

	pub fn add(&mut self,  gene_id:usize, start:usize ) {
		self.data.push( ( gene_id, start) );
	}

	pub fn get( &self, res:&mut HashMap<usize, usize>, pos:usize ) {
		for entry in &self.data {
			let value = res.entry(entry.0).or_insert(entry.1 - pos );
			if *value != entry.1 - pos {
				eprintln!("Re-match a gene, but with a fifferent start value! {} + {} and new {}  ",entry.0, value, entry.1 - pos);
			}
		}
	}

}



