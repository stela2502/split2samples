use crate::traits::Index;

use crate::errors::MappingError;

use crate::genes_mapper::gene_data::GeneData;
use crate::genes_mapper::gene_link::GeneLink;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

use crate::traits::BinaryMatcher;

use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::BTreeMap;

use core::fmt;

#[derive(Debug,PartialEq)]
pub struct GenesMapper{
	genes: Vec<GeneData>,
	mapper: Vec<GeneLink>,
	/// make sure no gene is added with the same sequence twice
	gene_hashes: HashSet<u64>,
	/// how many of my mapper entries have data?
	with_data: usize,
	/// if this is not the only mapper used - where should I start to count my genes?
	offset:usize,
	names: BTreeMap<String, usize>,
	/// to which hw value should a match be accepted?
	highest_nw_val: f32,
	/// how many u8's need to be the same for a read to match a gene?
	min_matches: usize,
	report4: Option<HashSet<usize>>, // report for a gene?
}

// Implementing Display trait for SecondSeq
impl fmt::Display for GenesMapper {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GenesMapper with {} mapper entries for {} genes", self.with_data, self.genes.len() )
    }
}

impl GenesMapper{

	pub fn new( offset: usize) -> Self{
		Self{
			genes: Vec::<GeneData>::with_capacity( 40_000 ),
			mapper: Vec::<GeneLink>::with_capacity(u16::MAX as usize),
			gene_hashes: HashSet::new(),
			with_data: 0,
			offset,
			names: BTreeMap::new(),
			highest_nw_val: 0.3,
			min_matches: 10, //40 bp exact match
			report4: None, // report for a gene?
		}
	}

	pub fn add(&mut self, seq: &[u8], name: String, chr: String, start:usize ) -> usize{
		let mut gene_data = GeneData::new( seq, &name, &chr, start );
		// Hash the gene_data object to obtain its hash value
	    let mut hasher = DefaultHasher::new();
	    gene_data.hash(&mut hasher);
	    let hash_value = hasher.finish();
		if self.gene_hashes.contains( &hash_value ){
			eprintln!("The sequence for gene {gene_data} has already been added before");
			return 0
		}
		
		self.gene_hashes.insert( hash_value );
		
		// add the sequence info and names
		let gene_id = self.genes.len();
		self.names.insert( name.to_string(), gene_id );
		self.genes.push(gene_data.clone());

		// add the mapper keys
		let mut keys = 0;
		while let Some( (key, offset )) = gene_data.next(){
			// check the key and ignore too simple ones:
			if &self.di_nuc_tab_length( &key ) < &4{ // could be as simple as AAAACAAA or ACACACAC
				continue;
			}
			if self.mapper[key as usize].is_empty() {
                // will add in the next step so
                self.with_data +=1;
            }
            self.mapper[key as usize].add( gene_id, offset  );
            keys += 1;
		}

		keys
	}

	/// For the multicore processing I need this object to be mergable
    /// All genes from the other object need to be incorporated and all primaryid - secondaryid compbos need to be collected.
    /// All other classes have to be copied, too. So we need to use our own incorporate_match_combo function.
    pub fn merge( &mut self, other: &Self ) {

        'main: for other_data in other.genes {
        	let mut gene_data = other_data.clone();
        	
        	let mut hasher = DefaultHasher::new();
		    gene_data.hash(&mut hasher);
		    let hash_value = hasher.finish();
			if self.gene_hashes.contains( &hash_value ){
				eprintln!("The sequence for gene {gene_data} has already been added before");
				continue 'main;
			}else {
				self.gene_hashes.insert( hash_value );
			}
			let gene_id = self.genes.len();
			self.genes.push(gene_data.clone());
			while let Some((key, offset)) = gene_data.next(){
				// check the key and ignore too simple ones:
				if &self.di_nuc_tab_length( &key ) < &4{ // could be as simple as AAAACAAA or ACACACAC
					continue;
				}
				if self.mapper[key as usize].is_empty() {
	                // will add in the next step so
	                self.with_data +=1;
	            }
	            self.mapper[key as usize].add( gene_id, offset );
			}
		}
	}

	pub fn get(&self, seq: &[u8] ) ->  Result< Vec<usize>, MappingError >{ 
		let mut gene_data = GeneData::new( seq, "read", "read2", 0 );
		let mut res = HashMap::<usize, usize>::new();

		while let Some((key, start)) = gene_data.next(){
			if ! self.mapper[key as usize].is_empty() {
				self.mapper[key as usize].get( &mut res, start );
			}
			for (gene_id, start) in &res {
				if let Some(cmp) = self.genes[*gene_id].slice( *start, gene_data.len() ){
					if cmp.needleman_wunsch( &gene_data, 0.3 ) < self.highest_nw_val {
						return Ok( vec![*gene_id + self.offset] )
					}
				}
			}
		}

		Err(MappingError::NoMatch)
	}

	pub fn get_strict(&self, seq: &[u8] ) ->  Result< Vec<usize>, MappingError >{ 
		let mut gene_data = GeneData::new( seq, "read", "read2", 0 );
		let mut res = HashMap::<usize, usize>::new();

		while let Some((key, start)) = gene_data.next(){
			if ! self.mapper[key as usize].is_empty() {
				self.mapper[key as usize].get( &mut res, start );
			}
			for (gene_id, start) in &res {
				if let Some(cmp) = self.genes[*gene_id].slice( *start, gene_data.len() ){
					if cmp.equal_entries( &gene_data ) < self.min_matches {
						return Ok( vec![*gene_id + self.offset] )
					}
				}
			}
		}

		Err(MappingError::NoMatch)
	}

	pub fn report4( &mut self, genes: &[&str]) {
        let mut hash = HashSet::<usize>::new();
        let mut added = false;
        for gname in genes{
            if let Some(gene_id) = self.extern_id_for_gname( gname ){
                hash.insert(gene_id);
                added = true;
            }
        }
        if added {
            self.report4 = Some( hash );
        }else {
            self.report4 = None;
        }
    }

    pub fn report4gene( &self, geneid:&[usize] ) -> bool {
        if let Some(hash) = &self.report4{
            for gid in geneid{
                if hash.contains(gid){
                    return true
                }
            }
            false
        }else {
            false
        }
    }

    pub fn extern_id_for_gname(&self, gname: &str ) -> Option<usize>{
    	self.names.get( gname ).map(|id| id + self.offset)
    }

	pub fn get_gene_count(&self) -> usize{
		self.genes.len()
	}

	pub fn get_all_gene_names( &self ) -> Vec<String> {
		self.names.keys().map(|name| name.to_string()).collect()
	}

	fn di_nuc_tab_length (&self, val: &u16 ) -> usize {
		let mut sum = vec![0_i8;16];

		for off in 0..7{
			let dinuc = (val >> (off * 2)) & 0b1111;
			sum[ dinuc as usize] +=1;
		}
        sum.iter().filter(|&x| *x != 0).count()
	}

	pub fn set_min_matches( &mut self, val:usize) {
		self.min_matches = val;
	}

	pub fn set_highest_nw_val(&mut self, val:f32){
		self.highest_nw_val = val;
	}
	
	fn get_id( &self, name: String ) -> usize{
		0
	}
	pub fn print( &self ){
		println!("{self}");
	}
	fn change_start_id ( &mut self, new_start :usize ){
		if self.offset == 0 {
            self.offset = new_start;
        }else {
            panic!("You try to change the start id twice - not supported!");
        }
	}
	fn names_len(&self) -> usize{
		0
	}
	fn names(&self) -> Vec<String>{
		vec!["Implement me!".to_string()]
	}
	fn to_header_n( &self, names: &[String] ) -> std::string::String{
		let mut ret= Vec::<std::string::String>::with_capacity( self.genes.len() +4 );
        //println!( "I get try to push into a {} sized vector", self.names.len());
        for obj in self.genes {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push(  obj.get_name().to_string() ) ;
        }
        ret.push("Most likely name".to_string());
        ret.push("Faction total".to_string());
        ret.push("dist to nr.2 [%max]".to_string());
        "CellID\t".to_owned()+&ret.join("\t")
	}
	fn max_id( &self ) -> usize { // return the max:id for the sparse export of the data
		self.genes.len()
	}
	pub fn write_index( &mut self, path: String ) -> Result< (), &str>{
		// this index needs to store the genes vectors and names. Genes binary and names as comma separated list?

		// the mapper would probably also make sense to store. Let's check the time we need to re-create that before storing it.
		Err("Not implemented")
	}
	pub fn load_index( &mut self, path: String ) -> Result< (), &str>{
		Err("Not implemented")
	}

	pub fn write_index_txt( &mut self, path: String ) -> Result< (), &str>{
		// this index needs to store the genes vectors and names. Genes binary and names as comma separated list?

		// the mapper would probably also make sense to store. Let's check the time we need to re-create that before storing it.
		Err("Not implemented")
	}

}