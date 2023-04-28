use std::collections::BTreeMap;


/// A mapper entry is a simplistic storage for the fast mapper class.
/// It mainly consists of a list of u64 sequences that can be used for mapping (32bp binary)

#[derive(Debug,PartialEq)]
pub struct MapperEntry{
	pub map:BTreeMap<u64, usize>, // the data storage
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
		self.map.insert( seq, id );
		self.only = id;
	}

	pub fn get( &self,seq:&u64 ) -> Option<usize> {
		if self.map.len() == 1{
			return Some(self.only)
		}
		self.map.get( seq ).copied()
	}

	pub fn has_data(&self) -> bool{
		! self.map.is_empty()
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

        assert_eq!( mapper.get(&12), Some(4) );
        assert_eq!( mapper.get(&44), Some(3) );

        assert_eq!( mapper.get(&14), None );
    }
}