// genes_mapper/gene_link.rs
use std::collections::HashMap;
use serde::{Serialize, Deserialize};

#[derive(PartialEq, Debug, Clone,Deserialize,Serialize)]
pub struct GeneLink {
	// store all chr position combos for this u16 8bp area.
	data: Vec<(usize, usize)>,
}

impl GeneLink{
	pub fn new()->Self {
		Self{
			data: Vec::with_capacity(10),
		}
	}

	pub fn clone( &self ) -> Self {
		Self{
			data: self.data.clone()
		}
	}

	pub fn is_empty(&self) -> bool{
		self.data.is_empty()
	}

	pub fn add(&mut self,  gene_id:usize, start:usize ) {
		self.data.push( ( gene_id, start) );
	}

	pub fn get( &self, res:&mut HashMap<(usize, i32), usize>, pos:i32 ) {
		for (gene_id, start_on_gene) in &self.data {
			res.entry(( *gene_id, *start_on_gene as i32 - pos )).and_modify(|e| {
		        // the count for this combo
		            *e += 1; // Increment the total_amount
		    }).or_insert( 1 );
		}
	}


	pub fn data(&self) -> std::slice::Iter<(usize, usize)>{
		self.data.iter()
	}

}



