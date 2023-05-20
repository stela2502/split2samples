use std::collections::HashSet;


/// A mapper entry is a simplistic storage for the fast mapper class.
/// It mainly consists of a list of u64 sequences that can be used for mapping (32bp binary)


#[derive(Debug,PartialEq)]
pub struct NameEntry{
	pub data:Vec<usize>, // the data storage
	pub significant_bp:Vec<usize>, // how much data should be matched against here(e.g. sample IDs are quite short)
	pos:usize,// the position to return with next
}

impl NameEntry{
	pub fn new() -> Self{
		let data = Vec::new();
		let significant_bp = Vec::new();
		let pos = 0;
		Self{ 
			data,
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
		// else if significant_bp < &32{
		// 	println!("I am returning a map for {} {} -> {}",self.significant_bp[self.pos], significant_bp ,self.significant_bp[self.pos].min(*significant_bp));
		// }
		//println!( "NameEntry::next() - I have {} significant bits for the iteration {}",self.significant_bp[self.pos]*2, self.pos);
		let mask: u64 = (1 << self.significant_bp[self.pos].min(*significant_bp) *2) -1;
		// println!("the returned map {mask:b}");
        self.pos += 1;
        Some(mask)
    }

    pub fn reset(&mut self) {
    	self.pos=0;
    } 

	pub fn add( &mut self, id:usize, sign:usize) {
		if sign == 0 {
			panic!("This does not make sense - I can not use a u64 with {sign} significant bits!");
		}
		if ! self.data.contains( &id ) {
			self.data.push(id);
			self.significant_bp.push(sign);
		}
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
	only:usize
}

impl MapperEntry{
	pub fn new() -> Self{
		let map = Vec::<(u64, NameEntry)>::new();
		let only =0;
		Self{
			map,
			only
		}
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

	pub fn add( &mut self, seq:u64, id:usize, sign:usize) {
		match self.find(&seq, 32) {
			Some( i ) =>  {
				self.map[i].1.add( id, sign );
			},
			None => {
				let mut name_entry = NameEntry::new();
				name_entry.add( id , sign);
				self.map.push( (seq, name_entry) );
			}
		};
		self.only = id;
	}

	/// finds the most likely matching entry in our set of sequences.
	/// This now matches - if the u64 does not match in total the u8 4bp kmers instead.
	/// But not continuousely as that would need me to convert them back to string.
	/// If this becomes necessary it can be added later.
	pub fn find (&mut self, seq:&u64, significant_bp:usize) -> Option<usize> {
		for i in 0..self.map.len() {
			self.map[i].1.reset(); // a NameEntry
			let mut bitmask = self.map[i].1.next(&significant_bp);
			// if self.map[i].1.contains(467){
			// 	while let Some(mask) = bitmask{
			// 		if mask != u64::MAX {
			// 			println!("I am processing using this bitmask: {:b}", mask);
			// 			let mut loc_seq = seq.clone();
			// 			loc_seq &= mask;
			// 			let mut loc_entry = self.map[i].0.clone();
			// 			loc_entry &= mask;
			// 			println!( "I had the seq {seq} and am now using {loc_seq} to match to my entry orig{} changed {}", *&self.map[i].0 , loc_entry);
			// 			if loc_seq == loc_entry {
			// 				println!("reducinge significant bits lead to a match to  {:?}", &self.map[i].1.data );
			// 				return Some(i)
			// 			}
			// 		}
			// 		else {
			// 			println!("Checking {} vs {}", seq,  &self.map[i].0 );
			// 			if seq == &self.map[i].0 {
			// 				println!("And I have a match to {i}");
			// 				return Some(i)
			// 			}
			// 		}
			// 		println!("But this did not turn out to be a match");
			// 		bitmask = self.map[i].1.next(&significant_bp);
			// 	}
			// }
			// else {
				while let Some(mask) = bitmask{
					if mask != u64::MAX {
						//println!("I am processing using this bitmask: {:b}", mask);
						let mut loc = seq.clone();
						loc &= mask;
						//println!( "I have the seq \n{seq:b} and am using this to map:\n{loc:b}");
						if loc == *&self.map[i].0 {
							//println!("reducinge significant bits lead to a match to  {:?}", &self.map[i].1.data );
							return Some(i)
						}
					}
					else {
						//println!("Checking {} vs {}", seq,  &self.map[i].0 );
						if seq == &self.map[i].0 {
							//println!("And I have a match to {i}");
							return Some(i)
						}
					}
					//println!("Some problem?");
					bitmask = self.map[i].1.next(&significant_bp);
				}
			//}
			
			// if seq == &self.map[i].0 {
			// 	//println!("Match to the internal seq {}", self.map[i].0 );
			// 	return Some(i);
			// }
		}
		// now check if we could use the to_le_bytes() on both u64's and find the one with the best sub-match
		// these sequences might have a polyA tail - cut that!
		let seq_u8 = seq.clone().to_le_bytes();

		let mut seq_other:[u8;8];
		let mut count = Vec::<usize>::with_capacity(self.map.len() );
		let mut max = 0;
		let mut id = 0;
		let mut entry_id = 0;
		let mut c:usize;
		for entry in &self.map {
			seq_other = entry.0.clone().to_le_bytes();
			c=0;
			//println!("I try to match the other {} ({:?}) to mine: {} or {:?}",seq, seq_u8, entry.0, seq_other);
			for i in 0..8{
				if seq_u8[i] == seq_other[i] {
					if seq_u8[i] == 0{
						continue; // not match AAAA as they could easily be 'not existent' in one of the oligos
					}
					//println!("\t\t\tmatch {}", seq_u8[i]);
					c +=1;
				}	
			}
			if max < c {
				max = c;
				entry_id = id;
			}
			count.push(c);
			id +=1;
		}
		id = usize::MAX;
		if max >1 {
			for i in 0..count.len() {
				if count[i] == max{
					if id != usize::MAX{
						// more than one entry has top matches
						// a classical multimapper to the end
						// return None!
						//println!("count found multiple matches - A multi mapper! -> returning None");
						return None
					}
					id = i;
				}
			}
			//println!("New count has identified a possible gene: {:?}", self.map[id].1.data );
			return Some(entry_id);
		}

		// Still no match - come on - a frame shift - I bet it's true
		let mut other:u64;
		id = 0;
		for entry in &self.map {
			other = entry.0.clone();
			other >>=2; //shift in the one direction
			seq_other = other.to_le_bytes();
			c=0;
			for i in 0..8{
				if seq_u8[i] == seq_other[i] {
					if seq_u8[i] == 0{
						continue; // not match AAAA as they could easily be 'not existent' in one of the oligos
					}
					//println!("\t\t\tmatch {}", seq_u8[i]);
					c +=1;
				}
			}
			if max < c {
				max = c;
				entry_id = id;
			}
			id +=1;
			count.push(c);
		}
		if max >1 {
			//println!("We found one match with a >>=!");
			return Some(entry_id);
		}
		// Still no match - come on - a frame shift - I bet it's true
		let mut other:u64;
		id = 0;
		for entry in &self.map {
			other = entry.0.clone();
			other <<=2; //shift in the one direction
			seq_other = other.to_le_bytes();
			c=0;
			for i in 0..8{
				if seq_u8[i] == seq_other[i] {
					if seq_u8[i] == 0{
						continue; // not match AAAA as they could easily be 'not existent' in one of the oligos
					}
					// println!("\t\t\tmatch {} to gene(s) {:?}", seq_u8[i], entry.1.data );
					c +=1;
				}
			}
			if max < c {
				max = c;
				entry_id = id;
			}
			id +=1;
			count.push(c);
		}
		if max >1 {
			//println!("We found one match with a <<=!");
			return Some(entry_id);
		}

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

	pub fn get( &mut self,seq:&u64, significant_bp:usize ) -> Option<Vec<usize>> {
		match self.find(seq, significant_bp){
			Some(id) => {
				return Some(self.map[id].1.data.clone())
			},
			None => return None,
		};
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

     #[test]
    fn check_geneids() {
        let mut mapper = MapperEntry::new();

        mapper.add(12, 4 );
        mapper.add(45, 3);

        assert_eq!( mapper.get(&12), Some(vec![4]) );
        assert_eq!( mapper.get(&45), Some(vec![3]) );

		mapper.add(12, 14 );
		assert_eq!( mapper.get(&12), Some(vec![4, 14]) );


        assert_eq!( mapper.get(&14), None );

        assert_eq!( mapper.info(), [2,1,1] );
    }
}