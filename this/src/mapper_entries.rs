use std::collections::BTreeMap;


/// A mapper entry is a simplistic storage for the fast mapper class.
/// It mainly consists of a list of u64 sequences that can be used for mapping (32bp binary)


#[derive(Debug,PartialEq)]
pub struct NameEntry{
	pub map:Vec<usize>, // the data storage
}

impl NameEntry{
	pub fn new() -> Self{
		let map = Vec::new();
		Self{ 
			map 
		}
	}

	pub fn add( &mut self, id:usize) {
		self.map.push(id);
	}
	pub fn get( &self ) -> Vec<usize>{
		self.map.clone()
	}
	pub fn get_mut( &mut self ) -> &mut Vec<usize>{
		&mut self.map
	}

}


#[derive(Debug,PartialEq)]
pub struct MapperEntry{
	pub map:BTreeMap<u64, NameEntry>, // the data storage
	only:usize
}

impl MapperEntry{
	pub fn new() -> Self{
		let  map = BTreeMap::new();
		let only =0;
		Self{
			map,
			only
		}
	}

	pub fn add( &mut self, seq:u64, id:usize) {
		match self.map.get_mut( &seq ){
			Some( name_entry ) =>  {
				name_entry.add( id );
			},
			None => {
				let mut name_entry = NameEntry::new();
				name_entry.add( id );
				self.map.insert( seq, name_entry );
			}
		};
		self.only = id;
	}

	pub fn get( &self,seq:&u64 ) -> Option<Vec<usize>> {
		Some(self.map.get( seq )?.get())
	}

	pub fn has_data(&self) -> bool{
		! self.map.is_empty()
	}
	/// first - how many entries
	/// second how many entries with one gene
	/// third how many entries with more than one gene
	pub fn info (&self) -> [usize;3] {
		let mut ret: [usize;3] = [0,0,0];
		for entry in self.map.values(){ //.cloned().collect() {
			ret[0] +=1;
			ret[1+ (entry.map.len() > 1) as usize] += 1;
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

	pub fn print( &self ){

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