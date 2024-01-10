use std::collections::HashSet;
use std::collections::HashMap;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// A mapper entry is a simplistic storage for the fast mapper class.
/// It mainly consists of a list of u64 sequences that can be used for mapping (32bp binary)


#[derive(Debug,PartialEq)]
pub struct NameEntry{
	// the gene ids storage - they identify the genes in the fast_mapper::names_store vector
	// in version 4 there is a change in data structure - the gene_id is now a tuple of
	// .0 = gene_id and .1 being the 'level' of the id - 0 -gene/transcript, 1 -> family, ->2 class etc.
	// It depends on how many different name variables the resp. create_index program took into account.
	pub data:Vec<(usize, usize)>,
	// a short vector with class info. This is needed for the TE scanner.
	// the ids point to the fast_mapper::
	pub classes:Vec<Vec<usize>>,
	// I want to quickly find if a dataset has alreadxy been added.
	hashes: HashSet<u64>, 
	pos:usize,// the position to return with next
}

impl NameEntry{
	pub fn new() -> Self{
		let data = Vec::with_capacity(4);
		let classes = Vec::with_capacity(4);
		let hashes = HashSet::new();
		let pos = 0;
		Self{ 
			data,
			classes,
			hashes,
			pos
		}
	}

	pub fn next(&mut self, significant_bp:&usize) -> Option<u64> {

		// we have returned all values
		if self.pos == self.data.len(){
			return None
		}

		if significant_bp >= &32{
			self.pos += 1;
			return Some(u64::MAX)
		}
		//else if significant_bp < &32{
		//	println!("I am returning a map for {} {} -> {}",self.significant_bp[self.pos], significant_bp ,self.significant_bp[self.pos].min(*significant_bp));
		//}
		//println!( "NameEntry::next() - I have {} significant bits for the iteration {}",self.significant_bp[self.pos]*2, self.pos);
		let mask: u64 = (1 << 32.min(*significant_bp) *2) -1;
		// println!("the returned map {mask:b}");
        self.pos += 1;
        Some(mask)
    }

    pub fn reset(&mut self) {
    	self.pos=0;
    }

    fn hash_classes(tup: &(usize, usize), vals: &Vec<usize> ) -> u64 {
	    let mut hasher = DefaultHasher::new();

	    tup.0.hash(&mut hasher);
	    tup.1.hash(&mut hasher);
	    // Hash the usize values
	    for val in vals {
	    	val.hash(&mut hasher);
	    }
	    // Finalize and return the hash value
	    hasher.finish()
	}

	pub fn add( &mut self, tup:(usize, usize), classes:Vec<usize>) -> bool {

		let hash = Self::hash_classes(&tup, &classes );

		if self.hashes.contains( &hash){
			// the combo of gene id, mapping class, [gene tags] has already been added here.
			// eprintln!("I alread had this data in my data storage: {:?}, {:?} -> {}: {:?}", tup, classes, hash, self.hashes.clone() );
			return false
		}
		//eprintln!("Mapper_entries: I add {:?} and {:?} with hash {}", tup, classes, hash);
		self.data.push( tup );
		self.classes.push(classes);
		self.hashes.insert( hash );
		true
	}

	/// this requires the classes vectors to be filled with optional names for the genes.
	/// like family or class names.
	/// If the classes vector is not populated and we have more than one gene matching here it will break.
	/// The classes vector should be filled with incresingly common entries line [ transcript_id, gene_id, family_id, class_id ]
	pub fn best_name_4_entries( &mut self ) -> Option<(usize, usize)> {

		// if this 40bp mapper links to only one gene the most likely return value is this gene_id
		if self.data.len() == 1{
			//println!("best_name_4_entries - I only have one! {}", self.data[0]);
			return Some(self.data[0]);
		}

		let mut counter: HashMap<usize, usize> = HashMap::new();
		if self.classes[0].is_empty(){
			panic!("This function can only be used while creating an index - I am missing the gene classes info here!")
		}
		// for every position in the classes vectors (should be sorted by expected rarity many ... view)
		for yid in 0..self.classes[0].len(){
			counter.clear();
			// check every classes vector and collect the id of the same rarity class
			for xid in 0..self.data.len(){
				if ! xid < self.classes.len() {
					self.classes.push(vec![0;self.classes[0].len()]);
				}
				*counter.entry(self.classes[xid][yid]).or_insert(0) += 1;
			}
			if counter.len() == 1{
				// we only got one ID in one rarity class - that is what is the best match!
				let best = *counter.keys().next().unwrap();
				//println!("best_name_4_entries - aaand the best one is {}",best);
				return Some( ( best, yid) )
			}
		}
		// so here we have not obtained any unique hit - best approach is possibly to return a None then - or?!
		None
	}

	pub fn get( &self ) -> Vec<(usize, usize)>{
		self.data.clone()
	}

	pub fn to_string( &self ) -> String {
		format!("NameEntry linking to {} gene_ids: {:?}", self.data.len(), self.data )
	}

	pub fn contains(&self, gene_id:usize) -> bool{
		for i in 0..self.data.len(){
			if self.data[i].0 == gene_id{
				return true
			}
		}
		return false
	}


}
