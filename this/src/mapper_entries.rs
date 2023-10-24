use std::collections::HashSet;
use std::collections::HashMap;

/// A mapper entry is a simplistic storage for the fast mapper class.
/// It mainly consists of a list of u64 sequences that can be used for mapping (32bp binary)


#[derive(Debug,PartialEq)]
pub struct NameEntry{
	// the gene ids storage - they identify the genes in the fast_mapper::names_store vector
	pub data:Vec<usize>,
	// a short vector with class info. This is needed for the TE scanner.
	// the ids point to the fast_mapper::
	pub classes:Vec<Vec<usize>>, 
	pub significant_bp:Vec<usize>, // how much data should be matched against here(e.g. sample IDs are quite short)
	pos:usize,// the position to return with next
}

impl NameEntry{
	pub fn new() -> Self{
		let data = Vec::with_capacity(4);
		let classes = Vec::with_capacity(4);
		let significant_bp = Vec::with_capacity(4);
		let pos = 0;
		Self{ 
			data,
			classes,
			significant_bp,
			pos
		}
	}

	pub fn next(&mut self, significant_bp:&usize) -> Option<u64> {

		// we have returned all values
		if self.pos == self.data.len(){
			return None
		}

		if self.significant_bp[self.pos] == 32 &&  significant_bp >= &32{
			self.pos += 1;
			return Some(u64::MAX)
		}
		//else if significant_bp < &32{
		//	println!("I am returning a map for {} {} -> {}",self.significant_bp[self.pos], significant_bp ,self.significant_bp[self.pos].min(*significant_bp));
		//}
		//println!( "NameEntry::next() - I have {} significant bits for the iteration {}",self.significant_bp[self.pos]*2, self.pos);
		let mask: u64 = (1 << self.significant_bp[self.pos].min(*significant_bp) *2) -1;
		// println!("the returned map {mask:b}");
        self.pos += 1;
        Some(mask)
    }

    pub fn reset(&mut self) {
    	self.pos=0;
    } 

	pub fn add( &mut self, gene_id:usize, sign:usize, classes:Vec<usize>) {
		if sign == 0 {
			panic!("This does not make sense - I can not use a u64 with {sign} significant bits!");
		}
		if ! self.data.contains( &gene_id ) {
			self.data.push(gene_id);
			self.significant_bp.push(sign);
			self.classes.push(classes);
		}else {
			eprintln!("This mapping region already links to the gene {} - ignored", gene_id );
		}
	}

	/// this requires the classes vectors to be filled with optional names for the genes.
	/// like family or class names.
	/// If the classes vector is not populated and we have more than one gene matching here it will break.
	/// The classes vector should be filled with incresingly common entries line [ transcript_id, gene_id, family_id, class_id ]
	pub fn best_name_4_entries( &mut self ) -> Option<usize> {

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
				return Some( best )
			}
		}
		// so here we have not obtained any unique hit - best approach is possibly to return a None then - or?!
		None
	}

	pub fn get( &self ) -> Vec<usize>{
		self.data.clone()
	}

	pub fn to_string( &self ) -> String {
		format!("NameEntry linking to {} gene_ids: {:?}", self.data.len(), self.data )
	}

	pub fn contains(&self, gene_id:usize) -> bool{
		for i in 0..self.data.len(){
			if self.data[i] >= gene_id{
				return true
			}
		}
		return false
	}


}


#[derive(Debug,PartialEq)]
pub struct MapperEntry{
	pub map:Vec::<(u64, NameEntry)>, // the data storage
	only:usize,
	hamming_cut: u32, // the bit difference up to which a match between two 32bp regions would still be acceptable.
}

impl MapperEntry{
	pub fn new() -> Self{
		let map = Vec::<(u64, NameEntry)>::with_capacity(4);
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

	pub fn possible_ids(&self) -> Vec<usize>{
		let mut ids = HashSet::<usize>::with_capacity(5);
		for entry in &self.map {
			for gid in &entry.1.data{
				ids.insert( *gid);
			}
		}
		ids.into_iter().collect()
	}

	/// add a match pair 8bp + <sign> bp.
	pub fn add( &mut self, seq:u64, id:usize, sign:usize, classes:Vec<usize>) {

		let mut added = false;
		for i in 0..self.map.len() {
			if self.map[i].0 == seq {
				self.map[i].1.add( id, sign, classes.clone());
				added = true;
				self.only = 0;
			}
		}
		if ! added{
			let mut name_entry = NameEntry::new();
			name_entry.add( id , sign, classes.clone() );
			self.map.push( ( seq, name_entry ) );
			if self.map.len() == 1{
				self.only = id;
			}else {
				self.only = 0;
			}
		}
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
		for i in 0..self.map.len() {
			if &self.map[i].0 == seq {
				if self.map[i].1.data.len() > 1{
					eprintln!("Ooops - we have a find in more than one gene: {:?}", self.map[i].1.data);
				}
				return Some( &self.map[i].1 )
			}
		}
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

#[cfg(test)]
mod tests {

    use crate::mapper_entries::MapperEntry;
    use crate::int_to_str::IntToStr;

    #[test]
    fn check_geneids() {
        let mut mapper = MapperEntry::new();

        //let tool = IntToStr::new(b"AGCTGTGAGACTCTTCACACTATCATCATTATTCGGAGG".to_vec(), 16);
        mapper.add(12, 4, 16 );
        mapper.add(45, 3, 16 );

        assert_eq!( mapper.get(&12), Some(vec![4]) );
        assert_eq!( mapper.get(&45), Some(vec![3]) );

		mapper.add(12, 14, 16 );
		assert_eq!( mapper.get(&12), Some(vec![4, 14]) );


        assert_eq!( mapper.get(&14), None );

        assert_eq!( mapper.info(), [2,1,1] );
    }

    #[test]
    fn hamming_distance() {
        let mut mapper = MapperEntry::new();

        let a:u64 = 0b11010101;
        let b:u64 = 0b10110001;

        assert_eq!( MapperEntry::hamming_distance(a, b), 3 as u32 );
    }


}