use std::collections::BTreeMap;



pub struct IndexedGenes{
	/// to store the name to external id of the genes
	names: BTreeMap<String, usize>,
	ids_to_name: Vec<String>,
	offset: usize,
}


impl IndexedGenes{

	pub fn new ( data: &BTreeMap<String, usize>, offset: usize ) -> Self {
		let mut ids_to_name = vec![String::new()];//;data.len() + offset];
		let mut names = BTreeMap::new();
		for (name, id) in data {
			ids_to_name[*id] = name.to_string(); 
			names.insert( name.to_string(), *id );
		}
		Self{
			names,
			ids_to_name,
			offset,
		}
	}

	/// returns the external id - the ID used in the single_cell_data class.
	pub fn ids_for_gene_names(&self, names:&Vec<String> ) -> Vec<usize>{
		let mut ret = Vec::<usize>::with_capacity( names.len() );
		for name in names{
			if let Some(id) = self.names.get( name ){
				ret.push(*id + self.offset)
			}else{
				panic!("The gene {name} is not found in this IndexedGenes object with names line {:?}", 
					&self.ids_to_name[0..3]);
			}
		}
		ret
	}

	/// return all genes we describe here
	pub fn get_all_gene_names(&self) -> Vec<String> {
		self.ids_to_name.clone()
	}

	/// get the header for dens matrices
	pub fn to_header_n( &self, names:&Vec<String> ) -> String{
		let mut ret= Vec::<std::string::String>::with_capacity( names.len() +5 );
        //println!( "I get try to push into a {} sized vector", self.names.len());
        for name in names {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push( name.to_string() ) ;
        }
        ret.push("AssignedSampleName".to_string());
        ret.push("FractionTotal".to_string());
        ret.push("n".to_string());
        ret.push("dist to nr.2 [%max]".to_string());
        "CellID\t".to_owned()+&ret.join("\t")
	}
}
