use crate::traits::Index;

use crate::errors::MappingError;
use crate::errors::GeneSelectionError;

use crate::genes_mapper::gene_data::GeneData;
use crate::genes_mapper::gene_link::GeneLink;
use crate::genes_mapper::Cigar;
use crate::genes_mapper::MapperResult;

use crate::singlecelldata::IndexedGenes;
use crate::traits::BinaryMatcher;


use std::hash::{Hash, Hasher};

use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::BTreeMap;

use regex::Regex;

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
	/// to which hw value should a match be accepted default to 0.2?
	highest_nw_val: f32,
	/// a pre-filter to the needleman wunsch test default to 0.6
	highest_humming_val: f32,
	/// how many u8's need to be the same for a read to match a gene?
	min_matches: usize,
	/// print the matching sequences for a gene in this list
	report4: Option<HashSet<usize>>, // report for a gene?
	/// additional mathcing and adding infos printed
	debug: bool,
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
			mapper: vec![GeneLink::new(); u16::MAX as usize],
			gene_hashes: HashSet::new(),
			with_data: 0,
			offset,
			names: BTreeMap::new(),
			highest_nw_val: 0.2,
			highest_humming_val: 0.6,
			min_matches: 10, //40 bp exact match
			report4: None, // report for a gene?
			debug:false,
		}
	}


	pub fn set_min_matches( &mut self, val:usize) {
		self.min_matches = val;
	}

	/// Function to set highest_nw_val
    pub fn set_highest_nw_val(&mut self, value: f32) {
        self.highest_nw_val = value;
    }

    /// Function to set highest_humming_val
    pub fn set_highest_humming_val(&mut self, value: f32) {
        self.highest_humming_val = value;
    }

	/// the IndexedGenes replace the full mapper class in the data export
    /// This allowes for multiple Indices to all feed the same data structure.
	pub fn as_indexed_genes(&self) -> IndexedGenes{
		IndexedGenes::new( &self.names, self.offset )
	}

	pub fn debug( &mut self, debug: Option<bool> ) -> bool{
		if let Some(dbg) = debug{
			self.debug = dbg;
		}
		self.debug
	}

	fn key_to_string( key:&u16 ) -> String{
        let mut data = String::new();
        for i in 0..8 {
            let ch = match (key >> (i * 2)) & 0b11 {
                0b00 => "A",
                0b01 => "C",
                0b10 => "G",
                0b11 => "T",
                _ => "N",
            };
            data += ch;
        }
        data
	}
	fn reject_key(&self,  key:&u16 ) -> bool{
		if &self.di_nuc_tab_length( &key ) < &4{ // could be as simple as AAAACAAA or ACACACAC
			if self.debug{
				println!("this key was rejected: {}", Self::key_to_string( &key ) );
			}
			return true;
		}
		if key > &65534 {
			if self.debug{
				println!("this key was rejected: {}", Self::key_to_string( &key ) );
			}
			return true
		}
		false
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
			if self.reject_key( &key ){
				continue
			}
			if self.debug{
				println!("I add the key {} for the gene {gene_data}", Self::key_to_string(&key) );
			}
			// or the add the info:
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

        'main: for other_data in &other.genes {
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

	fn slice_objects ( &self, start:i32, change_start: &GeneData, change_end: &GeneData ) -> Option<( GeneData, GeneData )> {
		
		let abs_start = start.abs() as usize;
		match start < 0{
			true => {
				let obj_a = match change_end.slice( abs_start, (change_end.len() - abs_start ).min( change_start.len() )){
					Some(val) => val,
					None => {
						eprintln!("I could not slice change_start start {start} length {} and seq {change_start}", (change_end.len() - abs_start ).min( change_start.len() ));
						return None
					},
				};
				let obj_b = match change_start.slice( 0, obj_a.len() ){
					Some(val) => val,
					None => {
						eprintln!("I could not slice change_end start 0 length {} and seq {change_end}", obj_a.len());
						return None
					},
				};
				Some((obj_b, obj_a))
			},
			false => {
				let obj_a = match change_start.slice( abs_start, (change_start.len() - abs_start ).min( change_end.len() )){
					Some(val) => val,
					None => {
						eprintln!("I could not slice change_start start {start} length {} and seq {change_start}", (change_start.len() - abs_start ).min( change_end.len() ));
						return None
					},
				};
				let obj_b = match change_end.slice( 0, obj_a.len() ){
					Some(val) => val,
					None => {
						eprintln!("I could not slice change_end start 0 length {} and seq {change_end}", obj_a.len());
						return None
					},
				};
				Some((obj_a, obj_b))
			}
		}	
	}

	/// this function adjusts the GeneData objects to get rid of terminal I's.
	fn fix_database_4_cigar( &self, gene_id:&usize, database: &GeneData, read: &GeneData, cigar: &mut Cigar, match_length_problems:&Regex ){

		if let Some(capture) = match_length_problems.captures( &cigar.cigar ) {
			let num_bp: usize = capture[1].parse().unwrap();
			eprintln!("updateing the cigar string {cigar}");
			match &capture[2] {
				"I" => {
					// the database was num_bp too large
					let adjusted_db = database.slice(0, database.len() - num_bp ).unwrap();
					let uw = read.needleman_wunsch( &adjusted_db, self.highest_humming_val, Some( cigar) );
					eprintln!("Type I: I compare the read {read} to original db {database} and the new databse {adjusted_db} \nand obtained the new cigar '{cigar}'");
				},
				"D" => {
					// the database was num_bp too short
					let adjusted_db = &self.genes[*gene_id].slice( database.get_start() - self.genes[*gene_id].get_start() , database.len() + num_bp ).unwrap();
					let uw = read.needleman_wunsch( &adjusted_db, self.highest_humming_val, Some( cigar) );
					eprintln!("Type D: I compare the read {read} to original db {database} and the new databse {adjusted_db} \nand obtained the new cigar '{cigar}'");
				},
				&_ => (),
			};
			eprintln!("updataed cigar string   = '{cigar}'");
		}
	}

	pub fn get(&self, read_data: &GeneData, res_vec: &Vec<(usize, Vec<i32>)>, cellid:u32 ) ->  Result< Vec<MapperResult>, MappingError >{

	    let mut ret = Vec::<MapperResult>::new();
	    let mut cigar = Cigar::new("");
		for (gene_id, start) in res_vec {

	    	match Self::all_values_same( start ){
	    		Ok(()) => {
	    			if let Some((read, database)) = self.slice_objects( start[0] , &self.genes[*gene_id], &read_data){
						if self.debug{
							println!("using the match gene id {gene_id} and start at {} I'll compare these two sequences", start[0]);
							println!("read \n{read} to database\n{database}\n");
						}
						if read.needleman_wunsch( &database, self.highest_humming_val, None ) < self.highest_nw_val {
							if self.report4this_gene( gene_id ) || self.debug{
								let nw = read.needleman_wunsch( &database, self.highest_humming_val, Some(&mut cigar) );
								println!("get        got gene_id {gene_id} or {} for the cell {cellid}, the cigar {cigar}, the start {}, the nw_value {nw} and the read {}", 
									self.genes[*gene_id].get_name(), start[0], read_data.to_dna_string() );
								ret.push( MapperResult::new( *gene_id + self.offset, start[0] as usize, true, Some(format!("{}",cigar)),
									 cigar.qual(), nw*read.len() as usize, cigar.edit_distance() ) );
							}else {
								// 
								ret.push( MapperResult::new( *gene_id + self.offset, start[0] as usize, false, None, 0, 0, 0 ) );
							}
							
						}
					}
					
				},
				Err(GeneSelectionError::NotSame) =>{
		    		let counts = Self::table( start );
		    		//eprintln!("I have multiple matching 8bp's but they match at different positions relative to the gene start!\n{start:?}");
		    		for (start, count) in counts{
		    			if count < 3 {
	    					break;
	    				}
	    				// if the objects are slicable and depending whether start > 0 or not we need to slice them differently
	    				if let Some((read, database)) = self.slice_objects( start , &self.genes[*gene_id], &read_data){

							if self.debug{
								println!("using the match start {start} count {count} I'll compare these two sequences");
								println!("read \n{read} to database\n{database}\n");
							}
							if read.needleman_wunsch( &database, self.highest_humming_val, None ) < self.highest_nw_val {
								if self.report4this_gene( gene_id ) || self.debug{
									// need to get the cigar string
									let nw= read.needleman_wunsch( &database, self.highest_humming_val, Some(&mut cigar) );
									println!("get        got gene_id {gene_id} or {} for the cell {cellid}, the cigar {cigar}, the start {}, the nw_value {nw} and the read {}", 
										self.genes[*gene_id].get_name(), start, read_data.to_dna_string() );
									ret.push( MapperResult::new( *gene_id + self.offset, start as usize, true, Some(format!("{}",cigar)), 
										cigar.qual(), nw*read.len() as usize, cigar.edit_distance() ) );
								}else {
									// 
									ret.push( MapperResult::new( *gene_id + self.offset, start as usize, false, None ,0, 0, 0 ) );
								}
							}
						}	
	    			}
	    		},
		    	_=> { // not enough matches in the first place!
		    	},
		    };
		}

		if self.debug{
			println!("mapping the read {read_data}\new obtained a initial mapping result {res_vec:?}\n and in the end END we found {ret:?}");
		}
		if ret.is_empty(){
			Err(MappingError::NoMatch)
		}else {
			Ok(ret)
		}
	}


	fn all_values_same(vec: &[i32]) -> Result<(), GeneSelectionError> {
		if vec.len() < 3 {
			return Err(GeneSelectionError::TooView);
		}
	    if let Some(first) = vec.first() {
	        if vec.iter().all(|x| x == first){
	        	Ok(())
	        }else {
	        	Err(GeneSelectionError::NotSame)
	        }
	    } else {
	        // If the vector is empty, technically all values are the same (trivially true)
	        // but useless here as filtered out before
	        Err(GeneSelectionError::TooView) // this should not be checked - too bad an initial match!!
	    }
	}

	/// this will return a sorted vector of (<start on source>:i32, count:usize)
	/// Highest counts first
	fn table(vec: &[i32]) -> Vec<(i32,usize)>{
		let count_map = vec.iter().fold(HashMap::new(), |mut map, &val| { 
	        *map.entry(val).or_insert(0) += 1; 
	        map 
    	});
    	let mut sorted_counts: Vec<_> = count_map.into_iter().collect();
    	sorted_counts.sort_by(|&(_, count1), &(_, count2)| count2.cmp(&count1));
    	sorted_counts
	}


	pub fn get_strict(&self, seq: &[u8], cellid:u32 ) ->  Result< Vec<MapperResult>, MappingError >{ 
		let mut read_data = GeneData::new( seq, "read", "read2", 0 );
		let mut res = HashMap::<usize, Vec<i32>>::new();

		let mut i = 0;
		while let Some((key, start)) = read_data.next(){
			if self.reject_key(&key){
				continue;
			}

			i+=1;
			if ! self.mapper[key as usize].is_empty() {
				self.mapper[key as usize].get( &mut res, start as i32 );
			}
			//println!("after {i} iterations ({}) the res looks like that {res:?}", Self::key_to_string( &key) );
		}
		
		
		if res.is_empty(){
			return Err(MappingError::NoMatch)
		}
		// Convert the HashMap into a vector of key-value pairs
	    let mut res_vec: Vec<_> = res.into_iter().collect();

	    // Sort the vector by the length of the internal arrays
	    // longest first
	    res_vec.sort_by(|(_, a), (_, b)| b.len().cmp(&a.len()));

	    let mut ret = Vec::<MapperResult>::new();
	    let mut cigar= Cigar::new("");
	    let match_length_problems = Regex::new(r"(\d+)(D|I)$").unwrap();

	    if let Some (entry) = &res_vec.first() {
	    	let gene_id = &entry.0;
	    	let start = &entry.1;
	    	match Self::all_values_same( start ){
	    		Ok(()) => {
	    			if let Some((read, database)) = self.slice_objects( start[0], &self.genes[*gene_id], &read_data ){
						if self.debug{
							println!("using the match {entry:?} I'll compare these two sequences");
							println!("read \n{read} to database\n{database}\n");
						}
						if database.equal_entries( &read ) >= self.min_matches {
							if self.report4this_gene( gene_id ) || self.debug{
								let nw = read.needleman_wunsch( &database, self.highest_humming_val, Some(&mut cigar) );
								// gene_id:usize, databsase: &GeneData, read: &GeneData, cigar: &mut Cigar, match_length_problems:Regex
								self.fix_database_4_cigar( gene_id, &database, &read, &mut cigar, &match_length_problems );
							
								println!("get_strict got gene_id {gene_id} or {} for the cell {cellid}, the cigar {cigar}, the start {}, the nw_value {nw} and the read {}", 
									self.genes[*gene_id].get_name(), start[0], read_data.to_dna_string() );
								ret.push( MapperResult::new( *gene_id + self.offset, start[0] as usize, true, Some(format!("{}",cigar)), 
									cigar.qual(), nw*read.len() as usize, cigar.edit_distance() ) );
							}else {
								// 
								ret.push( MapperResult::new( *gene_id + self.offset, start[0] as usize, false, None, 0, 0, 0 ) );
							}
						}
					};			
				},
				Err(GeneSelectionError::NotSame) =>{
		    		let counts = Self::table( start );
		    		if self.debug{
		    			println!("I have multiple matching 8bp's but they match at different positions relative to the gene start!\n{counts:?}");
		    		}
		    	
	    			//eprintln!("But there seams to be some higly likely entries: {counts:?}");
	    			for (start, count) in counts{
	    				if count < 3 {
	    					break;
	    				}
	    				if let Some((read, database)) = self.slice_objects( start, &self.genes[*gene_id], &read_data ){
							if self.debug{
								println!("using the match start {start} count {count} I'll compare these two sequences");
								println!("read \n{read} to database\n{database}\n");
							}
							if database.equal_entries( &read ) >= self.min_matches {
								if self.report4this_gene( gene_id ) || self.debug{
									let nw = read.needleman_wunsch( &database, self.highest_humming_val, Some(&mut cigar) );
									self.fix_database_4_cigar( gene_id, &database, &read, &mut cigar, &match_length_problems );
									
									println!("get_strict got gene_id {gene_id} or {} for the cell {cellid}, the cigar {cigar}, the start {}, the nw_value {nw} and the read {}", 
										self.genes[*gene_id].get_name(), start, read_data.to_dna_string() );
									ret.push( MapperResult::new( *gene_id + self.offset, start as usize, true, Some(format!("{}",cigar)), 
										cigar.qual(), nw*read.len() as usize, cigar.edit_distance() ) );
								}else {
									// 
									ret.push( MapperResult::new( *gene_id + self.offset, start as usize, false, None, 0 ,0, 0 ) );
								}
							}	
						}else {
							// I could not slice?!
							if self.debug{
								println!("I could not slice using {start} and read\n{read_data} and database\n{}",self.genes[*gene_id]);
							}
						}
		    		}
		    		
		    	},
		    	_=> { // not enough matches in the first place!
		    	},
		    };
		}

		if self.debug{
			println!("mapping the read {read_data}\nwe obtained a initial mapping result {res_vec:?}\n and in the end found {ret:?}");
		}
		if ret.is_empty(){
			self.get( &read_data, &res_vec, cellid )
			//Err(MappingError::NoMatch)
		}else {
			Ok(ret)
		}
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

    pub fn report4gene( &self, geneid:&Vec<usize> ) -> bool {
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

    pub fn report4this_gene( &self, geneid:&usize ) -> bool {
        if let Some(hash) = &self.report4{
            if hash.contains(geneid){
                true
            }else {
            	false
            }

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
        for obj in &self.genes {
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
		eprintln!("GenesMapper::write_index - NOT IMPLEMENTED!!!");
		Ok(())
	}
	pub fn load_index( &mut self, path: String ) -> Result< (), &str>{
		eprintln!("GenesMapper::load_index - NOT IMPLEMENTED!!!");
		Ok(())
	}

	pub fn write_index_txt( &mut self, path: String ) -> Result< (), &str>{
		// this index needs to store the genes vectors and names. Genes binary and names as comma separated list?

		// the mapper would probably also make sense to store. Let's check the time we need to re-create that before storing it.
		
		Ok(())
	}

}