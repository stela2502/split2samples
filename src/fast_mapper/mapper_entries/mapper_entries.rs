use std::collections::HashMap;
use std::collections::HashSet;
use crate::traits::BinaryMatcher;
use std::mem;

use crate::fast_mapper::mapper_entries::second_seq::SecondSeq;
use crate::fast_mapper::mapper_entries::NameEntry;

#[derive(Debug,PartialEq)]
pub struct MapperEntry{
	//pub map:Vec::<(SecondSeq, NameEntry)>, // the old data storage 
	pub map: HashMap<SecondSeq, NameEntry>, // to speed up the fast_mapper merge function
	//hamming_cut: u32, // the bit difference up to which a match between two 32bp regions would still be acceptable.
	//needleman_wunsch_cut: f32,
}

impl Default for MapperEntry {
    fn default() -> Self {
        Self::new()
    }
 }

impl MapperEntry{

	pub fn len(&self) -> usize{
		self.map.len()
	}

	pub fn new( ) -> Self{
		/*let mut all = allocate;
		if allocate < 4 {
			all = 4
		}*/
		let map = HashMap::new();
		Self{
			map,
		}
	}

	// Method to calculate memory size
    pub fn memory_size(&self) -> usize {
        let mut size = mem::size_of::<HashMap<SecondSeq, NameEntry>>(); // Size of HashMap<SecondSeq, NameEntry>
        
        // Iterate over each entry in the map and sum up the memory size of NameEntry instances
        for name_entry in self.map.values() {
            size += name_entry.memory_size();
        }
        
        // Add the size of other fields in MapperEntry
        size += mem::size_of::<u32>() // Size of u32 field (hamming_cut)
              + mem::size_of::<f32>(); // Size of f32 field (needleman_wunsch_cut)
        
        size
    }

	pub fn collapse_classes (&mut self ){
		// afterwards we should only have the mapper entries left that actually contained unique information.

		for key in self.map.keys().cloned().collect::<Vec<_>>(){
	    //for (key, name_entry) in &self.map {
	    	if let Some(name_entry) = self.map.get_mut(&key) {
		    	match name_entry.best_name_4_entries() {
		        	// gene_id is a (usize, usize)
		            Some(gene_id) => {
		            	*name_entry = NameEntry::new( name_entry.key );
		            	name_entry.add( gene_id.0, gene_id.1);
		            }
		            None => {
		            	let mut new_entry = NameEntry::new( name_entry.key );
		            	//new_entry.add( (gene_id.0.0, 20), gene_id.1 ); // I assume I will never get 20 classes of a gene - This will always FAIL!
		            	new_entry.keep = false;
		            	self.map.insert( key, new_entry);
		            }
		        }
		    }
		}

	}

	pub fn possible_ids(&self) -> Vec<(usize, usize)>{
		let mut ids = HashSet::<(usize, usize)>::with_capacity(5);
		for entry in self.map.values() {
			for gid in &entry.data{
				ids.insert( *gid);
			}
		}
		ids.into_iter().collect()
	}

	/// add a u64 kmer_size fragment as second key.
	/// here we also need the sign amount of bp that contain informative sequence
	pub fn add( &mut self, seq:SecondSeq, id:(usize,usize), classes:Vec<usize>) -> bool{

	    for (_, entry) in self.map.iter_mut() {
	        if entry.key.same(&seq) {
	        	return entry.add( id, classes );
	        }
	    }

	    // we did not have the same entry here!
	    let mut new_entry = NameEntry::new( seq );
	    _ = new_entry.add(id, classes.clone());
	    self.map.insert( seq, new_entry );
	    true
	}


	/// get is the exact match whereas find is a somewhat fuzzy match.
	/// So if get does not find anything at all - try find instead.
	pub fn get( &self,seq:&SecondSeq) -> Option<(Vec<&NameEntry>, f32)> {
		
		
		let ret:Vec::<&NameEntry> = self.map.iter().filter_map( |(key, name_entry)| 
				if key.same( seq ){
					//println!("Match in mapper_entry!");
					Some(name_entry)
				}
				else {
					//println!("No match in mapper_entry!");
					None
				} 
			).collect();
		
		/*
		let mut ret : Vec::<&NameEntry> = vec![];
		for (key, name_entry) in self.map.iter() {
			if key.same( &seq ) {
				//println!("the id {seq} has a 100% matching sequence: {key} ");
				ret.push(name_entry);
			}
		}
		*/
		if ret.len() == 1 {
			Some((ret, 0.0))
		}
		else {
			None
		}
	}

	/// finds the most likely matching entry in our set of sequences.
	/// It returns the vector of matching gene ids and the needleman_wunsch
	/// matching value it got for them.
	pub fn find (&self, seq:&SecondSeq, nw_cut:f32, humming_cut:f32  ) -> Option<(Vec<&NameEntry>, f32)> {
		let mut ret : Vec::<&NameEntry> = vec![];

		let mut dists : Vec::<f32> = vec![];
		let mut min_dist: f32 = f32::MAX;
		for name_entry in self.map.values().filter_map(| name_entry | {
		    if name_entry.key.tri_nuc_abs_diff(seq) < humming_cut {
		        Some(name_entry)
		    } else {
		        None
		    }
		}) {
		    let dist = name_entry.key.needleman_wunsch(seq, nw_cut, None);
		    if dist <= nw_cut {
		        ret.push(name_entry);
		        dists.push(dist);
		        if dist < min_dist {
		            min_dist = dist;
		        }
		    }
		}
		match ret.len() {
			0 =>{
				//eprintln!("I find nothing here!");
				None
			}
			1 =>{
				//eprintln!("I find exactly one!");
				Some((ret, min_dist))
			},
			len =>{
				let mut ret2: Vec::<&NameEntry> = vec![];
				for i in 0..len{
					if dists[i] == min_dist{
						ret2.push( ret[i] )
					}
					
				}
				//eprintln!("I found multiple and {:?} with the lowest difference {}",ret2, min_dist );
				Some((ret2, min_dist))
			}
		}

	}

	pub fn print(&self) {
		if  self.has_data() {
			println!("I have {} u64 matching sequences:", self.map.len() );
			for entry in self.map.values(){
				println!( "\tThe sequence {} links to the {}", entry.key, entry );
			}
		}
	}

	

	pub fn has_data(&self) -> bool{
		let mut ok = false;
		for name_entry in self.map.values() {
			if name_entry.keep {
				ok = true;
				break
			}
		}
		ok
	}

	pub fn with_data( &self) -> usize{
		let mut with_data = 0;
		for name_entry in self.map.values(){
			if name_entry.keep {
				with_data +=1 
			}
		}
		with_data
	}
	/// first - how many entries
	/// second how many entries with one gene
	/// third how many entries with more than one gene
	pub fn info (&self) -> [usize;3] {
		let mut ret: [usize;3] = [0,0,0];
		for entry in self.map.values(){
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
    	let result = match this.map.get(&name_entry.key){
    		Some(entry) => entry,
    		None => &name_entry,
    	};

    	assert_eq!(result.classes.len(), 1);
    	this.collapse_classes();

    	let result = match this.map.get(&name_entry.key){
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
    	let result = match this.map.get(&name_entry.key){
    		Some(entry) => entry,
    		None => &name_entry,
    	};
    	assert_eq!(result.classes.len(), 2);
		this.collapse_classes();
		let result = match this.map.get(&name_entry.key){
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
    	let result = match this.map.get(&name_entry.key){
    		Some(entry) => entry,
    		None => &name_entry,
    	};

    	assert_eq!(result.classes.len(), 2);
    	this.collapse_classes();

    	let result = match this.map.get(&name_entry.key){
    		Some(entry) => entry,
    		None => &name_entry,
    	};
    	assert_eq!(result.classes.len(), 1);
    	assert_eq!(result.data.len(), 1);
    	assert_eq!(result.data[0], (2_usize, 2_usize));
    }
    
}
