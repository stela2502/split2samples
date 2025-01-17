use std::collections::HashSet;
use std::collections::HashMap;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use core::fmt;
//use std::mem;

use crate::fast_mapper::mapper_entries::SecondSeq;

/// A mapper entry is a simplistic storage for the fast mapper class.
/// It mainly consists of a list of u64 sequences that can be used for mapping (32bp binary)

#[derive(Debug )]
pub struct NameEntry{
	// the gene ids storage - they identify the genes in the fast_mapper::names_store vector
	// in version 4 there is a change in data structure - the gene_id is now a tuple of
	// .0 = gene_id and .1 being the 'level' of the id - 0 -gene/transcript, 1 -> family, ->2 class etc.
	// It depends on how many different name variables the resp. create_index program took into account.
	pub key: SecondSeq,
	pub data:Vec<(usize, usize)>,
	// a short vector with class info. This is needed for the TE scanner.
	// the ids point to the fast_mapper::
	pub classes:Vec<Vec<usize>>,
	// I want to quickly find if a dataset has alreadxy been added.
	hashes: HashSet<u64>, 
	pos:usize,// the position to return with next
	min_level: usize, // where to start with checking for unique names? 
	pub keep:bool, // used in the filtering process
}


// Implementing Clone manually for NameEntry because it contains a HashSet
impl Clone for NameEntry {
    fn clone(&self) -> Self {
        Self {
        	key: self.key,
            data: self.data.clone(),
            classes: self.classes.clone(),
            hashes: self.hashes.clone(),
            pos: self.pos,
            min_level : self.min_level,
            keep: self.keep, //when clearing out unneccessary ones.
        }
    }
}



impl PartialEq for NameEntry {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}


impl Eq for NameEntry {}

impl Hash for NameEntry {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.key.hash(state);
    }
}

// Implementing Display trait for SecondSeq
impl fmt::Display for NameEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "NameEntry linking to {} gene_ids: {:?}", self.data.len(), self.data )
    }
}


impl NameEntry{
	pub fn new( key:SecondSeq) -> Self{
		let data = Vec::with_capacity(4);
		let classes = Vec::with_capacity(4);
		let hashes = HashSet::new();
		let pos = 0;
		let keep = true;
		Self{ 
			key,
			data,
			classes,
			hashes,
			pos,
			min_level:0,
			keep,
		}
	}

	// Method to calculate memory size
    pub fn memory_size(&self) -> usize {
    	0
        /*let size = mem::size_of::<SecondSeq>() // Size of SecondSeq field 
                 + mem::size_of::<usize>() *3 // Size of usize fields (length, min_level and pos)
                 + self.data.capacity() * mem::size_of::<(usize, usize)>() // Size of data vector's elements
                 + mem::size_of::<Vec<Vec<usize>>>() // Size of Vec<Vec<usize>>
                 + self.classes[0].capacity() * mem::size_of::<usize>()
                 + mem::size_of::<HashSet<u64>>() + self.hashes.capacity() * mem::size_of::<u64>() // Size of HashSet<u64>
                 + mem::size_of::<bool>(); // Size of bool field (keep)
        size*/
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
		let mask: u64 = (1 << (32.min(*significant_bp) *2)) -1;
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
			//println!("I alread had this data in my data storage: gene id and class{:?}, more classes {:?} -> this combos hash: {}: my hashes: {:?}", tup, classes, hash, self.hashes.clone() );
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
	/// This basically returns a new "self.data" entry.
	pub fn best_name_4_entries( &self ) -> Option<((usize, usize), Vec<usize>)> {

		if ! self.keep {
			return None
		}
		if self.classes[0].is_empty(){
			if self.data.len() == 1{
				return Some( (self.data[0], vec![self.data[0].0] ) )
			}else {
				return None
			}
		}

		// if this 40bp mapper links to only one gene the most likely return value is this gene_id
		if self.data.len() == 1{
			//println!("best_name_4_entries - I only have one! {}", self.data[0]);
			if self.min_level == 0{
				return Some((self.data[0], self.classes[0].clone()) )
			}
		}

		let mut counter: HashMap<usize, usize> = HashMap::new();

		let min_y = self.data.iter().map(|&(_, second)| second).max().unwrap_or(0);
		
		// for every position in the classes vectors (should be sorted by expected rarity many ... view)

		for yid in min_y..self.classes[0].len(){
			counter.clear();
			// check every classes vector and collect the id of the same rarity class
			for xid in 0..self.data.len(){
				*counter.entry(self.classes[xid][yid]).or_insert(0) += 1;
			}
			if counter.len() == 1{
				// we only got one ID in one rarity class - that is what is the best match!
				let best = *counter.keys().next().unwrap();
				//println!("best_name_4_entries - aaand the best one is {}",best);
				return Some( ( ( best, yid), self.classes[0].clone() ) )
			}
		}
		// so here we have not obtained any unique hit - best approach is possibly to return a None then - or?!
		None
	}

	pub fn get( &self ) -> Vec<(usize, usize)>{
		self.data.clone()
	}


	pub fn contains(&self, gene_id:usize) -> bool{
		for i in 0..self.data.len(){
			if self.data[i].0 == gene_id{
				return true
			}
		}
		false
	}

	pub fn same(&self, seq:&SecondSeq) -> bool{
		self.key.same( seq )
	}


}


#[cfg(test)]
mod tests {
	use crate::fast_mapper::mapper_entries::NameEntry;
	use crate::fast_mapper::mapper_entries::SecondSeq;
    #[test]
    fn check_add() {
    	let seq= SecondSeq(12345_u64, 20_u8);
    	let mut this = NameEntry::new( seq );
    	// say we Have a gene x with family Y and class Z
    	this.add( (0_usize,0_usize), vec!(0_usize, 1_usize, 2_usize));
    	assert_eq!( this.best_name_4_entries(), Some(((0_usize,0_usize), vec!(0_usize, 1_usize, 2_usize)) ) ); 	
    }
    #[test]
    fn check_best_name_4_entries_1() {
    	let seq= SecondSeq(12345_u64, 20_u8);
    	let mut this = NameEntry::new(seq);
    	// say we Have a gene x with family Y and class Z
    	this.add( (0_usize,0_usize), vec!(0_usize, 1_usize, 2_usize));
    	this.add( (3_usize,0_usize), vec!(3_usize, 1_usize, 2_usize));
    	assert_eq!( this.best_name_4_entries(), Some(((1_usize,1_usize), vec!(0_usize, 1_usize, 2_usize)) ) ); 	
    }
    #[test]
    fn check_best_name_4_entries_2() {
    	let seq= SecondSeq(12345_u64, 20_u8);
       	let mut this = NameEntry::new( seq );
    	// say we Have a gene x with family Y and class Z
    	this.add( (0_usize,0_usize), vec!( 0_usize, 1_usize, 2_usize));
    	this.add( (3_usize,0_usize), vec!( 3_usize, 4_usize, 2_usize));
    	assert_eq!( this.best_name_4_entries(), Some(((2_usize,2_usize), vec!(0_usize, 1_usize, 2_usize))) ); 	
    }
    #[test]
    fn check_same_not_entered_twice() {
    	let seq= SecondSeq(12345_u64, 20_u8);
       	let mut this = NameEntry::new( seq );
    	// say we Have a gene x with family Y and class Z
    	assert_eq!( this.add( (0_usize,0_usize), vec!( 0_usize, 1_usize, 2_usize)), true );
    	assert_eq!( this.add( (0_usize,0_usize), vec!( 0_usize, 1_usize, 2_usize)), false);

    	assert_eq!( this.data.len(), 1);
    	assert_eq!( this.best_name_4_entries(), Some(((0_usize,0_usize), vec!(0_usize, 1_usize, 2_usize))) ); 	
    }
}