use std::collections::BTreeMap;
use core::fmt;

use regex::Regex;


pub struct IndexedGenes{
	/// to store the name to external id of the genes
	names: BTreeMap<String, usize>,
	ids_to_name: Vec<String>,
	offset: usize,
}

// Implementing Display trait for IndexedGenes
impl fmt::Display for IndexedGenes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    	let first = 5.min( self.names.len() );
    	let names: Vec<String> = (0..first).map( |i| self.ids_to_name[ i ].to_string() ).collect();
        write!(f, "IndexedGenes for {} entries and gene names like {:?}",  self.names.len(), names )
    }
}

impl IndexedGenes{

	pub fn empty ( off:Option<usize> ) -> Self {

		Self{
			names: BTreeMap::new(),
			ids_to_name: Vec::with_capacity( 80_000 ),
			offset : match off {
				Some(o) => o,
				None => 0,
			},
		}

	}

	/// End pattern match subsetting of the gene names structure
	pub fn subset(&self, regex: &Regex, offset: usize ) -> IndexedGenes {
        // Use enumerate to build both ids_to_name and names
        let mut new_names = BTreeMap::new();
        let new_ids_to_name: Vec<String> = self
            .ids_to_name
            .iter()
            .filter(|name| regex.is_match(name))
            .map(|name| {
                let new_index = new_names.len(); // Current size of new_names gives the new index
                new_names.insert(name.clone(), new_index); // Populate new_names
                name.clone() // Return the name for the Vec
            })
            .collect();

        IndexedGenes {
            names: new_names,
            ids_to_name: new_ids_to_name,
            offset: 0, // Reset offset for the new subset
        }
    }


	/// merges the other list into this list and returns a Vec of new ids for the other data 
	pub fn merge (&mut self, other: &IndexedGenes) -> Vec<usize>{
		let mut rec = Vec::<usize>::with_capacity( other.len() );
		for gene in &other.ids_to_name {
			rec.push ( self.get_gene_id( gene ) )
		}
		rec
	}


	/// Returns the gene ID for the given gene name.
    /// If the gene does not exist, assigns a new ID.
    pub fn get_gene_id(&mut self, gene: &str) -> usize {

        // Check if the gene already exists in the map
        if let Some(&gene_id) = self.names.get(gene) {
            return gene_id;
        }

        // Gene doesn't exist; assign a new ID
        let new_id = self.offset + self.ids_to_name.len();
        self.names.insert(gene.to_string(), new_id);
        self.ids_to_name.push(gene.to_string());
        new_id
    }

    pub fn len( &self ) -> usize{
    	self.names.len()
    }

	pub fn new ( data: &BTreeMap<String, usize>, offset: usize ) -> Self {
		let mut ids_to_name = vec![String::new(); data.len() + offset ];
		let mut names = BTreeMap::new();
		let mut max_entry = 0_usize;
		let _ = data.iter().for_each(| (_,val) | if *val > max_entry {max_entry = *val});
		println!("I have gotten a max extry of {max_entry} and have a vec with {} available spaces",ids_to_name.len() );
		if max_entry > (data.len() + offset) {
			panic!("This is a library error - why is my largest id {max_entry} when it should be {}?!?",  data.len() + offset);
		}

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
	pub fn ids_for_gene_names(&self, names:&Vec::<String> ) -> Vec<usize>{
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
