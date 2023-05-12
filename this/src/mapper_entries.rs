//use std::collections::BTreeMap;


/// A mapper entry is a simplistic storage for the fast mapper class.
/// It mainly consists of a list of u64 sequences that can be used for mapping (32bp binary)


#[derive(Debug,PartialEq)]
pub struct NameEntry{
	pub data:Vec<usize>, // the data storage
}

impl NameEntry{
	pub fn new() -> Self{
		let data = Vec::new();
		Self{ 
			data
		}
	}

	pub fn add( &mut self, id:usize) {
		if ! self.data.contains( &id ) {
			self.data.push(id);
		}
	}
	pub fn get( &self ) -> Vec<usize>{
		self.data.clone()
	}

	pub fn to_string( &self ) -> String {
		format!("NameEntry linking to {} gene_ids: {:?}", self.data.len(), self.data )
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

	pub fn add( &mut self, seq:u64, id:usize) {
		match self.find(&seq) {
			Some( id ) =>  {
				self.map[id].1.add( id );
			},
			None => {
				let mut name_entry = NameEntry::new();
				name_entry.add( id );
				self.map.push( (seq, name_entry) );
			}
		};
		self.only = id;
	}

	/// finds the most likely matching entry in our set of sequences.
	/// This now matches - if the u64 does not match in total the u8 4bp kmers instead.
	/// But not continuousely as that would need me to convert them back to string.
	/// If this becomes necessary it can be added later.
	pub fn find (&self, seq:&u64 ) -> Option<usize> {
		for i in 0..self.map.len() {
			if seq == &self.map[i].0 {
				println!("Match to the internal seq {}", self.map[i].0 );
				return Some(i);
			}
		}
		// now check if we could use the to_le_bytes() on both u64's and find the one with the best sub-match
		let seq_u8 = seq.clone().to_le_bytes();
		let mut seq_other:[u8;8];
		let mut count = Vec::<usize>::with_capacity(self.map.len() );
		let mut max = 0;
		let mut id = 0;
		let mut entry_id = 0;
		for entry in &self.map {
			seq_other = entry.0.clone().to_le_bytes();
			let mut c = 0;
			println!("I try to match the other {} to mine: {} or {:?}",seq, entry.0, seq_u8);
			for i in 0..8{
				if seq_u8[i] == seq_other[i] {
					println!("\t\t\tmatch {}", seq_u8[i]);
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
			return Some(entry_id);
		}
		return None
	}

	pub fn print(&self) {
		if  self.has_data() {
			//println!("I have {} u64 matching sequences:", self.map.len() );
			for entry in &self.map{
				println!( "\tThe sequence {} links to the {}", entry.0, entry.1.to_string() );
			}
		}
	}

	pub fn get( &self,seq:&u64 ) -> Option<Vec<usize>> {
		match self.find(seq){
			Some(id) => return Some(self.map[id].1.data.clone()),
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
        mapper.add(44, 3);

        assert_eq!( mapper.get(&12), Some(vec![4]) );
        assert_eq!( mapper.get(&44), Some(vec![3]) );

		mapper.add(12, 14 );
		assert_eq!( mapper.get(&12), Some(vec![4, 14]) );


        assert_eq!( mapper.get(&14), None );

        assert_eq!( mapper.info(), [2,1,1] );
    }
}