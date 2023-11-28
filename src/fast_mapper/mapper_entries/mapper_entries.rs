use std::collections::HashSet;

use crate::fast_mapper::mapper_entries::NameEntry;

#[derive(Debug,PartialEq)]
pub struct MapperEntry{
	pub map:Vec::<(u64, NameEntry)>, // the data storage
	only:usize,
	hamming_cut: u32, // the bit difference up to which a match between two 32bp regions would still be acceptable.
}

impl MapperEntry{
	pub fn new( allocate: usize) -> Self{
		let mut all = allocate;
		if allocate < 4 {
			all = 4
		}
		let map = Vec::<(u64, NameEntry)>::with_capacity(all);
		let only =0;
		let hamming_cut = 4;
		Self{
			map,
			only,
			hamming_cut,
		}
	}

	pub fn collapse_classes (&mut self ){
		// afterwards we should only have the mapper entries left that actually contained unique information.
		self.map.retain_mut( |(_second_seq, name_entry)|  match name_entry.best_name_4_entries(){
			Some(gene_id) => {
				name_entry.data = vec![gene_id];
				name_entry.classes = Vec::<Vec<usize>>::with_capacity(0);
				true
			},
			None => false,
		} );
	}

	pub fn possible_ids(&self) -> Vec<(usize, usize)>{
		let mut ids = HashSet::<(usize, usize)>::with_capacity(5);
		for entry in &self.map {
			for gid in &entry.1.data{
				ids.insert( *gid);
			}
		}
		ids.into_iter().collect()
	}

	/// add a match pair 8bp + <sign> bp.
	pub fn add( &mut self, seq:u64, id:(usize,usize), classes:Vec<usize>) -> bool{

		for i in 0..self.map.len() {
			if self.map[i].0 == seq {
				self.only = 0;
				return  self.map[i].1.add( id, classes.clone());
			}
		}
		// now we have no match to the seq and therefore need to add one
		let mut name_entry = NameEntry::new();
		if name_entry.add( id , classes.clone() ){
			self.map.push( ( seq, name_entry ) );
			if self.map.len() == 1{
				self.only = id.0;
			}else {
				self.only = 0;
			}
		}
		
		true
	}

	/// calculate the bit flips between two u64 sequences
	pub fn hamming_distance(x: u64, y: u64) -> u32 {
    	(x ^ y).count_ones()
	}

	/// get is the exact match whereas find is a somewhat fuzzy match.
	/// So if get does not find anything at all - try find instead.
	pub fn get( &self,seq:&u64 ) -> Option<&NameEntry> {

		for i in 0..self.map.len() {
			if &self.map[i].0 == seq {
				// if self.map[i].1.data.len() > 1{
				// 	eprintln!("Ooops - we have a get in more than one gene: {:?}", self.map[i].1.data);
				// }
				return Some( &self.map[i].1 )
			}
		}
		// now we have an initial match, but no secondary.
		None
	}

	pub fn get_mut( &mut self,seq:&u64 ) -> Option<&mut NameEntry> {

		for i in 0..self.map.len() {
			if &self.map[i].0 == seq {
				// if self.map[i].1.data.len() > 1{
				// 	eprintln!("Ooops - we have a get in more than one gene: {:?}", self.map[i].1.data);
				// }
				return Some( &mut self.map[i].1 )
			}
		}
		// now we have an initial match, but no secondary.
		None
	}
	/// finds the most likely matching entry in our set of sequences.
	/// This now matches - if the u64 does not match in total the u8 4bp kmers instead.
	/// But not continuousely as that would need me to convert them back to string.
	/// If this becomes necessary it can be added later.
	pub fn find (&self, seq:&u64 ) -> Option<&NameEntry> {
		// for i in 0..self.map.len() {
		// 	if &self.map[i].0 == seq {
		// 		if self.map[i].1.data.len() > 1{
		// 			eprintln!("Ooops - we have a find in more than one gene: {:?}", self.map[i].1.data);
		// 		}
		// 		return Some( &self.map[i].1 )
		// 	}
		// }
		// now we have an initial match, but no secondary.
		for i in 0..self.map.len() {
			//eprintln!("Hamming distance below {} - returning {:?}", self.hamming_cut, self.map[i].1.data );
			if MapperEntry::hamming_distance( self.map[i].0, *seq ) < self.hamming_cut{
				// let mut a:String = String::from("");
				// let mut b:String = String::from("");
				// tool.u64_to_str(32, seq, &mut a );
				// tool.u64_to_str(32, &self.map[i].0, &mut b);
				// eprintln!("I have a match between sequence \nA{:?}\nB{:?}\nhamming_dist: {}\nsignificant_bp: {significant_bp}\ngenes: {:?}",a, b,MapperEntry::hamming_distance( self.map[i].0, *seq ), self.map[i].1.data );
			
				return Some( &self.map[i].1 )
			}
			// else {
			// 	let mut a:String = String::from("");
			// 	let mut b:String = String::from("");
			// 	tool.u64_to_str(32, seq, &mut a );
			// 	tool.u64_to_str(32, &self.map[i].0, &mut b);
			// 	eprintln!("I have no match between sequence \nA{:?}\nB{:?}\nhamming_dist: {}",a, b,MapperEntry::hamming_distance( self.map[i].0, *seq ) );
			// }
		}
		// if this still has not worked - is the RNA fragments possibly shorter than our genomic sequence here?
		for i in 0..self.map.len() {
			let our = self.map[i].0.to_be_bytes();
			let other = seq.to_be_bytes();
			let mut sum = 0;
			for id in 0..8{
				if our[id] == other[id] {
					sum += 1;
				}
			}
			if sum > 2 { // at least 24 bp exact match + position
				return Some( &self.map[i].1 )
			}
		}
		// so now we could have a frameshift due to a long stretch of a single nucleotide
		// if necessary implement that later on.
		// For now I think we should just stop here.
		return None
	}

	pub fn print(&self) {
		if  self.has_data() {
			println!("I have {} u64 matching sequences:", self.map.len() );
			for entry in &self.map{
				println!( "\tThe sequence {} links to the {}", entry.0, entry.1.to_string() );
			}
		}
	}

	

	pub fn has_data(&self) -> bool{
		! self.map.is_empty()
	}
	/// first - how many entries
	/// second how many entries with one gene
	/// third how many entries with more than one gene
	pub fn info (&self) -> [usize;3] {
		let mut ret: [usize;3] = [0,0,0];
		for i in 0..self.map.len(){
			ret[0] +=1;
			ret[1+ (self.map[i].1.data.len() > 1) as usize] += 1;
			// match entry.map.len() > 1{
			// 	true => {
			// 		ret[2] += 1;
			// 	},
			// 	false => {
			// 		ret[1] += 1;
			// 	}
			// }
		}
		ret
	}

}

