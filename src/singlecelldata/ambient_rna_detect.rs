use std::collections::BTreeMap;

use crate::singlecelldata::cell_data::GeneUmiHash;

/// Ambient RNA is RNA that more or less floats through the library creation space,
/// being captured by a multitiude of cells. This seams to lead to multiple cells
/// obtaining the same Gene+UMI combination.
/// This class can stack up the GeneUmiHashes and count their occurance over the whole dataset.

#[derive(Debug)]
pub struct AmbientRnaDetect{    
    /// The data storage 
    ambient: BTreeMap< GeneUmiHash, usize>,
    /// Has the data been finalized
    checked: bool
}


impl Default for AmbientRnaDetect {
    fn default() -> Self {
        Self::new()
    }
}

// here the functions
impl AmbientRnaDetect{
	pub fn new() -> Self {
		let ambient = BTreeMap::new();
		Self{
			ambient,
			checked: false,
		}
	}

	pub fn add( &mut self, data:GeneUmiHash) -> usize  {

		let counts = self.ambient.entry(data).or_insert(0);
		*counts +=1;
		*counts
	}

	pub fn finalize( &mut self, counts:usize ) {

		if self.checked{
			println!("This had already been checked");
		}else {
			self.ambient.retain( |&_key, size| *size >= counts );
			self.checked = true;
		}

	}

	pub fn is_ambient( &self, data:&GeneUmiHash) -> bool {

		if ! self.checked {
			panic!("This AmbientRnaDetect object has to be finalized first!")
		}
		self.ambient.contains_key( data )
	}

	pub fn clone( &self ) -> Self {
		if ! self.checked{
			panic!("Pleasedo not clone a not finalized object!");
		}
		
		Self {
            ambient: self.ambient.clone(),
            checked: true,
        }
	}

	pub fn size( &self ) -> usize{
		self.ambient.len()
	}

	pub fn entries( &self) -> Vec<&GeneUmiHash> {
		self.ambient.keys().collect()
	}
}


