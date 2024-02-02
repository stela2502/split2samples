use std::collections::HashSet;

use crate::fast_mapper::mapper_entries::second_seq::SecondSeq;
use crate::fast_mapper::mapper_entries::NameEntry;

#[derive(Debug,PartialEq)]
pub struct MapperEntry{
	//pub map:Vec::<(SecondSeq, NameEntry)>, // the old data storage 
	pub map: HashSet<NameEntry>, // to speed up the fast_mapper merge function
	only:usize,
	hamming_cut: u32, // the bit difference up to which a match between two 32bp regions would still be acceptable.
	needleman_wunsch_cut: f32,
}

impl MapperEntry{
	pub fn new( ) -> Self{
		/*let mut all = allocate;
		if allocate < 4 {
			all = 4
		}*/
		let map = HashSet::new();
		Self{
			map,
			only :0,
			hamming_cut :2,
			needleman_wunsch_cut: 0.4 // you want 0.3 there to not get a lot of crap - but I need more values - I need to try this.
		}
	}

	pub fn collapse_classes (&mut self ){
		// afterwards we should only have the mapper entries left that actually contained unique information.
		let mut new_map = HashSet::new();

	    for name_entry in &self.map {
	        match name_entry.best_name_4_entries() {
	        	// gene_id is a (usize, usize)
	            Some(gene_id) => {
	            	let mut new_entry = NameEntry::new( name_entry.key );
	                new_entry.add( gene_id.0, gene_id.1);
	                new_map.insert(new_entry);
	            }
	            None => {
	            }
	        }
	    }

	    // Remove entries marked for removal
	    self.map = new_map;
	    
	}

	pub fn possible_ids(&self) -> Vec<(usize, usize)>{
		let mut ids = HashSet::<(usize, usize)>::with_capacity(5);
		for entry in &self.map {
			for gid in &entry.data{
				ids.insert( *gid);
			}
		}
		ids.into_iter().collect()
	}

	/// add a u64 kmer_size fragment as second key.
	/// here we also need the sign amount of bp that contain informative sequence
	pub fn add( &mut self, seq:SecondSeq, id:(usize,usize), classes:Vec<usize>) -> bool{

		let mut new_entry = NameEntry::new( seq );

	    for entry in &self.map {
	        if entry.key.same(&new_entry.key) {
	        	for i in 0..entry.data.len(){
	        		_= new_entry.add( entry.data[i], entry.classes[i].clone() );
	        	}
	            break;
	        }
	    }

	    let ret = new_entry.add(id, classes.clone());
	    self.map.remove( &new_entry );
	    self.map.insert( new_entry );
	    
	    ret
	}


	/// get is the exact match whereas find is a somewhat fuzzy match.
	/// So if get does not find anything at all - try find instead.
	pub fn get( &self,seq:&SecondSeq ) -> Option<Vec<&NameEntry>> {

		let mut ret : Vec::<&NameEntry> = vec![];
		
		let new_entry = NameEntry::new( *seq );

		match &self.map.get( &new_entry ){
			Some(entry) => {
				ret.push(entry);
				Some(ret)
			},
			None => {
				None
			}
		}

	}

	/// finds the most likely matching entry in our set of sequences.
	/// This now matches - if the u64 does not match in total the u8 4bp kmers instead.
	/// But not continuousely as that would need me to convert them back to string.
	/// If this becomes necessary it can be added later.
	pub fn find (&self, seq:&SecondSeq ) -> Option<Vec<&NameEntry>> {
		let mut ret : Vec::<&NameEntry> = vec![];

		let mut dists : Vec::<f32> = vec![];
		let mut min_dist: f32 = f32::MAX;
		for name_entry in self.map.iter() {
			//eprintln!("Hamming distance below {} - returning {:?}", self.hamming_cut, self.map[i].1.data );
			//let dist = self.map[i].0.hamming_distance( seq );
			let dist = name_entry.key.needleman_wunsch( seq );
			//eprintln!("Distance is {dist}");
			//if dist <= self.hamming_cut {
			if dist <= self.needleman_wunsch_cut {
				// look at the matches that are almost rejected.
				/*if dist > self.needleman_wunsch_cut * 0.9 {
					println!( "{seq} fastq did match to \n{} database with {} - should that be right?\n", self.map[i].0, dist);
				}*/
				
				ret.push( name_entry );
				dists.push( dist );
				if dist < min_dist{
					min_dist = dist
				}
			}
		}
		return match ret.len() {
			0 =>{
				//eprintln!("I find nothing here!");
				None
			}
			1 =>{
				//eprintln!("I find exactly one!");
				Some(ret)
			},
			len =>{
				let mut ret2: Vec::<&NameEntry> = vec![];
				for i in 0..len{
					if dists[i] == min_dist{
						ret2.push( ret[i] )
					}
					
				}
				//eprintln!("I found multiple and {:?} with the lowest difference {}",ret2, min_dist );
				Some(ret2)
			}
		}

	}

	pub fn print(&self) {
		if  self.has_data() {
			println!("I have {} u64 matching sequences:", self.map.len() );
			for entry in &self.map{
				println!( "\tThe sequence {} links to the {}", entry.key, entry.to_string() );
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
		for entry in &self.map{
			ret[0] +=1;
			ret[1+ (entry.data.len() > 1) as usize] += 1;
		}
		ret
	}

}


#[cfg(test)]
mod tests {
	use crate::fast_mapper::mapper_entries::NameEntry;
	use crate::fast_mapper::mapper_entries::MapperEntry;
	use crate::fast_mapper::mapper_entries::SecondSeq;

    #[test]
    fn check_add() {
    	let mut this = MapperEntry::new();
    	// say we Have a gene x with family Y and class Z
    	this.add(SecondSeq(32_u64, 20) , (0_usize,0_usize), vec!(0_usize, 1_usize, 2_usize));

    	let name_entry = NameEntry::new( SecondSeq(32_u64, 20) );
    	let result = match this.map.get(&name_entry){
    		Some(entry) => entry,
    		None => &name_entry,
    	};

    	assert_eq!(result.classes.len(), 1);
    	this.collapse_classes();

    	let result = match this.map.get(&name_entry){
    		Some(entry) => entry,
    		None => &name_entry,
    	};
    	assert_eq!( this.map.len(), 1);
    	assert_eq!(result.data.len(), 1);
    	assert_eq!(result.classes.len(), 1);

    }
    
    #[test]
    fn check_collapse_classes_1() {
    	let mut this = MapperEntry::new();
    	// say we Have a gene x with family Y and class Z
    	this.add( SecondSeq(32_u64, 20), (0_usize,0_usize), vec!(0_usize, 1_usize, 2_usize));
    	this.add( SecondSeq(32_u64, 20), (3_usize,0_usize), vec!(3_usize, 1_usize, 2_usize));
    	let name_entry = NameEntry::new( SecondSeq(32_u64, 20) );
    	let result = match this.map.get(&name_entry){
    		Some(entry) => entry,
    		None => &name_entry,
    	};
    	assert_eq!(result.classes.len(), 2);
		this.collapse_classes();
		let result = match this.map.get(&name_entry){
    		Some(entry) => entry,
    		None => &name_entry,
    	};
    	assert_eq!(result.classes.len(), 1);
    	assert_eq!(result.data.len(), 1);
    	assert_eq!(result.data[0], (1_usize, 1_usize));
    }
    #[test]
    fn check_collapse_classes_2() {
    	let mut this = MapperEntry::new();
    	// say we Have a gene x with family Y and class Z
    	this.add( SecondSeq(32_u64, 20), (0_usize,0_usize), vec!( 0_usize, 1_usize, 2_usize));
    	this.add( SecondSeq(32_u64, 20), (3_usize,0_usize), vec!( 3_usize, 4_usize, 2_usize));
    	let name_entry = NameEntry::new( SecondSeq(32_u64, 25) );
    	let result = match this.map.get(&name_entry){
    		Some(entry) => entry,
    		None => &name_entry,
    	};

    	assert_eq!(result.classes.len(), 2);
    	this.collapse_classes();

    	let result = match this.map.get(&name_entry){
    		Some(entry) => entry,
    		None => &name_entry,
    	};
    	assert_eq!(result.classes.len(), 1);
    	assert_eq!(result.data.len(), 1);
    	assert_eq!(result.data[0], (2_usize, 2_usize));
    }
    
}
