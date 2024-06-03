//use crate::traits::Index;

use crate::errors::MappingError;
use crate::errors::GeneSelectionError;

use crate::genes_mapper::gene_data::GeneData;
use crate::genes_mapper::gene_link::GeneLink;
use crate::genes_mapper::Cigar;
use crate::genes_mapper::MapperResult;
use crate::genes_mapper::MultiMatch;
use crate::genes_mapper::CigarEndFix;

use crate::singlecelldata::IndexedGenes;
use crate::traits::BinaryMatcher;

use crate::genes_mapper::NeedlemanWunschAffine;

use std::hash::{Hash, Hasher};

use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::BTreeMap;

use std::fs::File;
use std::io::Read;

use serde::{Serialize, Deserialize};

use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::BufWriter;
use std::io::Write;

//use std::process::exit;

use core::fmt;

#[derive(Debug,PartialEq,Deserialize,Serialize)]
pub struct GenesMapper{
	genes: Vec<GeneData>,
	mapper: Vec<GeneLink>,
	/// make sure no gene is added with the same sequence twice and report the orig name of the entry
	gene_hashes: BTreeMap<u64, String>,
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
	/// if there are extremely small entries in this database we (sample tags or antibody tags) we need to know that
	small_entries:bool,
	/// an internal version that should be iterated one anything important changes.
	version:usize,

}

// Implementing Display trait for SecondSeq
impl fmt::Display for GenesMapper {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    	let first = 5.min(self.genes.len());
    	let names: Vec<String> = (0..first).map( |i| self.genes[i].get_name().to_string() ).collect();
        write!(f, "GenesMapper with {} mapper entries (offset {}) for {} genes like {:?}", self.with_data, self.offset, self.genes.len(), names )
    }
}

impl GenesMapper{

	pub fn new( offset: usize) -> Self{
		Self{
			genes: Vec::<GeneData>::with_capacity( 40_000 ),
			mapper: vec![GeneLink::new(); u16::MAX as usize],
			gene_hashes: BTreeMap::new(),
			version: 8,
			with_data: 0,
			offset: offset,
			names: BTreeMap::new(),
			highest_nw_val: 0.2,
			highest_humming_val: 0.6,
			min_matches: 10, //40 bp exact match
			report4: None, // report for a gene?
			debug:false,
			small_entries: false,
		}
	}

	pub fn len(&self) -> usize {
		self.genes.len()
	}

	pub fn depth(&self) -> usize {
		let mut sum = 0;
		for mapper in &self.mapper{
			if ! mapper.is_empty() {
				sum+=1;
			}
		}
		sum
	}

	pub fn set_small_entries( &mut self ) {
		self.small_entries = true;
	}

	/// this function could possibly be a good idea to get rid of repeat infos
	pub fn make_index_te_ready(&mut self ) {
		eprintln!("make_index_te_ready - Please implement me!")
	}

	// sam_header will return both the header string as well as the fasta database to that
	pub fn sam_header(&self) -> Option<(String, String)>{
		
		if let Some(hash) = &self.report4{
			let mut vec = Vec::<String>::with_capacity( hash.len() );
			let mut fasta = "".to_string();
			for &id in hash {
				let gene_name = &self.genes[id].get_unique_name();
				let gene_len = self.genes[id].len();
				let formatted_line = format!("@SQ\tSN:{}\tLN:{}", gene_name, gene_len);
				vec.push(formatted_line);
				fasta+= &self.genes[id].to_fasta();

			}
			let ret = vec.join("\n") + "\n";
			Some((ret, fasta))
		}else {
			let mut vec = Vec::<String>::with_capacity( self.genes.len() );
			let mut fasta = "".to_string();
			for  gene_obj in &self.genes {
				let gene_name = gene_obj.get_unique_name();
				let gene_len = gene_obj.len();
				let formatted_line = format!("@SQ\tSN:{}\tLN:{}", gene_name, gene_len);
				vec.push(formatted_line);
				fasta+= &gene_obj.to_fasta();
			}
			let ret = vec.join("\n") + "\n";
			Some((ret, fasta))
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
    	println!("I'll try to index my genes using the offset {}", self.offset);
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


	/// the main add function to add a gene into the index. 
	pub fn add(&mut self, seq: &[u8], unique_name:&str, name: &str, chr: &str, start:usize ) -> usize{

		let mut gene_data = GeneData::new( seq, &unique_name, name, chr, start );
		// Hash the gene_data object to obtain its hash value
		let mut hasher = DefaultHasher::new();
		gene_data.hash(&mut hasher);
		let hash_value = hasher.finish();
		if let Some(other_name) = self.gene_hashes.get( &hash_value ){
			#[cfg(debug_assertions)]
			eprintln!("The sequence for gene {gene_data} has already been added before: {other_name}");
			return 0
		}
		
		self.gene_hashes.insert( hash_value, name.to_string() );
		
		// add the sequence info and names
		// the above test only checks if we have the same gene + sequence in the database.
		// therefore we here need to check the name once more!
		let gene_id;
		if self.names.contains_key( name ){
			let mut i =1;
			let mut name2 = format!("{}.{}", name.to_string(),i);
			while self.names.contains_key( &name2 ){
				i +=1;
				name2 = format!("{}.{}", name,i);
			}
			gene_id = self.genes.len();
			self.names.insert( name2.to_string(), gene_id );
			gene_data = GeneData::new( seq, unique_name, &name2, chr, start );
			self.genes.push(gene_data.clone());
		}else {
			gene_id = self.genes.len();
			self.names.insert( name.to_string(), gene_id );
			self.genes.push(gene_data.clone());
		}

		

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
    		if let Some(other_name) = self.gene_hashes.get( &hash_value ){
    			#[cfg(debug_assertions)]
    			eprintln!("The sequence for gene {gene_data} has already been added before {other_name}");
    			continue 'main;
    		}else {
    			self.gene_hashes.insert( hash_value, gene_data.get_name().to_string() );
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

				if self.debug{
					println!("slicing with start <0 {start} and \n{change_start}\nand\n{change_end}");
					println!("I got you this return values:\n{obj_b}\n{obj_a}\ngood?\n")
				}
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
				#[cfg(debug_assertions)]
				if self.debug{
					println!("slicing with start >0 {start} and \n{change_start}\nand\n{change_end}");
					println!("I got you this return values:\n{obj_a}\n{obj_b}\ngood?\n")
				}
				Some((obj_a, obj_b))
			}
		}	
	}

	// /// this function adjusts the GeneData objects to get rid of terminal I's.
	// fn fix_database_4_cigar( &self, gene_id:&usize, database: &GeneData, read: &GeneData, cigar: &mut Cigar, match_length_problems:&Regex ){

	// 	if let Some(capture) = match_length_problems.captures( &cigar.cigar ) {
	// 		let num_bp: usize = capture[1].parse().unwrap();
	// 		//eprintln!("updateing the cigar string {cigar}");
	// 		match &capture[2] {
	// 			"I" => {
	// 				// the database was num_bp too large
	// 				let adjusted_db = database.slice(0, database.len() - num_bp ).unwrap();
	// 				let _uw = read.needleman_wunsch( &adjusted_db, self.highest_humming_val, Some( cigar) );
	// 				//eprintln!("Type I: I compare the read {read} to original db {database} and the new databse {adjusted_db} \nand obtained the new cigar '{cigar}'");
	// 			},
	// 			"D" => {
	// 				// the database was num_bp too short
	// 				let adjusted_db = &self.genes[*gene_id].slice( database.get_start() - self.genes[*gene_id].get_start() , database.len() + num_bp ).unwrap();
	// 				let _uw = read.needleman_wunsch( &adjusted_db, self.highest_humming_val, Some( cigar) );
	// 				//eprintln!("Type D: I compare the read {read} to original db {database} and the new databse {adjusted_db} \nand obtained the new cigar '{cigar}'");
	// 			},
	// 			&_ => (),
	// 		};
	// 		//eprintln!("updataed cigar string   = '{cigar}'");
	// 	}
	// }

	

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

	fn as_dna_string(val:&u16) -> String {
        let mut data = String::new();
        //println!("converting u64 {loc:b} to string with {kmer_size} bp.");
        for i in 0..8 {
            // Use a mask (0b11) to extract the least significant 2 bits
            let ch = match (val >> (i * 2)) & 0b11 {
                0b00 => "A",
                0b01 => "C",
                0b10 => "G",
                0b11 => "T",
                _ => "N",
            };
            data += ch;
            
            //println!("{ch} and loc {loc:b}");
        }

        //println!("\nMakes sense? {:?}", data);
        data
    }

    /*fn query_database( &self, seq: &[u8] ) -> Vec<((usize, i32), usize)> {
    	//if self.small_entries {
    		return self.query_database_small( seq );
    	//}
		let mut res = HashMap::<(usize, i32), usize>::new();// Vec<i32>>::new();
		let mut read_data = GeneData::new( seq,"read_t", "read", "read2", 0 );
		let mut i = 0;
		let mut failed = 0;

		while let Some((key, start)) = read_data.next(){

			//println!("after {i} (matching) iterations ({}) the res looks like that {res:?}", Self::key_to_string( &key) );
			if self.reject_key(&key){
				continue;
			}
			
			if ! self.mapper[key as usize].is_empty() {
				i+=1;
				failed = 0;
				#[cfg(debug_assertions)]
				if self.debug {
					let values: Vec<_> = self.mapper[key as usize].data().collect();
					let m = 9.min( values.len() );
					println!("I am testing this key: {} and got {:?} ...", Self::as_dna_string(&key) , &values[0..m] );
				}
				// this will directly modify the res HashMap as it can also find more than one!
				self.mapper[key as usize].get( &mut res, start as i32 );
			}
			else {
				failed +=1;
				if failed == 3 {
					match read_data.iterator_skip_to_next_frame() {
						Ok(_) => failed = 0,
						Err(_) => break,
					}
				}
				#[cfg(debug_assertions)]
				if self.debug{
					println!("I am testing this key: {} and got Nothing", Self::as_dna_string(&key) );
				}
			}
			if i == 30 {
				break;
			}
		}
		let mut res_vec: Vec<_> = res.into_iter().collect();
		res_vec.retain(|(_, counts)| *counts > 4);
	    // longest first
	    res_vec.sort_by(|( (a_key, _a_rel_start), a_counts), ((b_key, _b_rel_start), b_counts)| {
	    	match b_counts.cmp(&a_counts) {
		        std::cmp::Ordering::Equal => a_key.cmp(b_key),
		        other => other,
		    }
		});

		res_vec
    }

	fn query_database_small( &self, seq: &[u8] ) -> Vec<((usize, i32), usize)> {
    */
    
    fn query_database( &self, seq: &[u8] ) -> Vec<((usize, i32), usize)> {
		let mut res = HashMap::<(usize, i32), usize>::new();// Vec<i32>>::new();
		let mut read_data = GeneData::new( seq, "read_t", "read", "read2", 0 );
		let mut i = 0;

		#[cfg(debug_assertions)]
		println!("I am runing query_database_small");

		while let Some((key, start)) = read_data.next(){

			//println!("after {i} (matching) iterations ({}) the res looks like that {res:?}", Self::key_to_string( &key) );
			if self.reject_key(&key){
				continue;
			}
			
			if ! self.mapper[key as usize].is_empty() {
				i+=1;
				#[cfg(debug_assertions)]
				if self.debug {
					let values: Vec<_> = self.mapper[key as usize].data().collect();
					let m = 9.min( values.len() );
					println!("I am testing this key: {} and got {:?} ...", Self::as_dna_string(&key) , &values[0..m] );
				}
				// this will directly modify the res HashMap as it can also find more than one!
				self.mapper[key as usize].get( &mut res, start as i32 );
			}
			else {
				#[cfg(debug_assertions)]
				if self.debug{
					println!("I am testing this key: {} and got Nothing", Self::as_dna_string(&key) );
				}
			}
			if i == 30 {
				break;
			}
		}
		let mut res_vec: Vec<_> = res.into_iter().collect();

		#[cfg(debug_assertions)]
		{
			res_vec.retain(|(_, counts)| *counts > 2);
			res_vec.sort_by(|( (a_key, _a_rel_start), a_counts), ((b_key, _b_rel_start), b_counts)| {
	    	match b_counts.cmp(&a_counts) {
		        std::cmp::Ordering::Equal => a_key.cmp(b_key),
		        other => other,
		    }
		});
			println!("And I got the initial results {res_vec:?}");
		}

		res_vec.retain(|(_, counts)| *counts > 4);
	    // longest first
	    res_vec.sort_by(|( (a_key, _a_rel_start), a_counts), ((b_key, _b_rel_start), b_counts)| {
	    	match b_counts.cmp(&a_counts) {
		        std::cmp::Ordering::Equal => a_key.cmp(b_key),
		        other => other,
		    }
		});

		res_vec
    }

    pub fn get_strict(&self, seq: &[u8], _cellid:u32, nwa: &mut NeedlemanWunschAffine ) ->  Result< Vec<MapperResult>, MappingError >{ 
		let read_data = GeneData::new( seq,"read_t", "read", "read2", 0 );
		// store the gene id, the relative start on that and the count for this combo
		let res_vec = self.query_database( seq );

		if res_vec.is_empty() {
			return Err(MappingError::NoMatch)
		}

		let mut helper = MultiMatch::new();
		
	    //res_vec.sort_by(|(_, a), (_, b)| b.len().cmp(&a.len()));
	    #[cfg(debug_assertions)]
		if self.debug {
			let to=9.min(res_vec.len());
			println!("I have collected these initial matches (first 10): {:?}", &res_vec[0..to] );
		}
	    let mut cigar= Cigar::new("");
	    cigar.set_debug( nwa.debug() ); // propagate the debug setting from the nwa object
	    //panic!("remind me what I get here: {res_vec:?}");

		

		// collect all possible matches
		#[allow(unused_variables)]
	    for ((gene_id, start), count) in &res_vec {

	    	if let Some((read, database)) = self.slice_objects( *start, &self.genes[*gene_id], &read_data ){
	    		cigar.clear();
	    		if (read.len() as f32) < (read_data.len() as f32 * 0.8) && (read.len() as f32) < (self.genes[*gene_id].len() as f32 * 0.9) {
	    			// this database match is a little short!
	    			#[cfg(debug_assertions)]
					if self.debug{
						println!("Mappable sequence is too short: {}", read_data.len() );
					}
	    			continue;
	    		}
				let nw = &nwa.needleman_wunsch_affine( &read, &database, self.highest_humming_val  );

				#[cfg(debug_assertions)]
				if self.debug{
					println!("#################################################################################");
					println!("using the match to gene #{gene_id} with rel start {start} and {count} key matches supporting that, I'll compare these two sequences:\n");
					println!("read \n{read}\nto database\n{database}\n");
					println!("the alignement:\n{}",nwa.to_string( &read, &database, self.highest_humming_val ));
				}
				
				cigar.convert_to_cigar( &nwa.cigar_vec() );
				cigar.clean_up_cigar(&read, &database);

				if nw.abs() < self.highest_nw_val  {
					#[cfg(debug_assertions)]
					if self.debug{
						println!("################## And I deem this match intereting");
					}
					if cigar.mapping_quality() > 20 && cigar.state_changes() < 10  {
						helper.push( 
							MapperResult::new( 
									*gene_id + self.offset, *start.max(&0) as usize, 
								true, Some(cigar.clone()), 
								cigar.mapping_quality(), *nw, (nw*read.len() as f32) as usize,
								cigar.edit_distance(), self.genes[*gene_id].get_name(), self.genes[*gene_id].len()
							)
						);
					}

					cigar.clear();
				}
				#[cfg(debug_assertions)]
				if self.debug{
					println!("I got this nw: {nw} and the matches: {}",helper);
				}
			};			

		} // end populating helper
		#[cfg(debug_assertions)]
		if self.debug{
			//let res_lines = res_vec.iter().map(|obj| format!("{:?}", obj)).collect::<Vec<String>>().join("\n");
			println!("mapping the read {read_data}\nwe obtained an initial mapping result \n{helper}\nand in the end found \n{:?}",
				helper.get_best( seq.len()) );
			println!("#################################################################################");
		}
		match helper.get_best( seq.len()){
			Ok(val) => {
				//println!("I found a best result! {}", &val);
				if let Some(cigar) = val.cigar(){
					if cigar.mapping_quality() > 20 && cigar.fixed != Some(CigarEndFix::Both) {
						#[cfg(debug_assertions)]
						println!("get_strict got an accepted match (get_strict) {}", val );
						Ok(vec![val])
					}else {
						#[cfg(debug_assertions)]
						println!("secondary checks failed: {} <= 20 || {:?} != Some(CigarEndFix::Both) ", cigar.mapping_quality(), cigar.fixed);
						Err(MappingError::NoMatch)
						//self.get( &read_data, &res_vec, cellid, nwa)
					}
				}else {
					#[cfg(debug_assertions)]
					println!("The match had no acciociated cigar!?!" );
					Err(MappingError::NoMatch)
					//self.get( &read_data, &res_vec, cellid, nwa)
				}			
			},
			Err(_) => {
				#[cfg(debug_assertions)]
				println!("No best mapper identified!");
				Err(MappingError::NoMatch)
				//self.get( &read_data, &res_vec, cellid, nwa)
			}, //that is kind of OK
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

	pub fn get_gene(&self, gene_id:usize) -> Option<&GeneData> {
		if gene_id < self.genes.len(){
			Some(&self.genes[gene_id])
		}else {
			None
		}
	}

	pub fn get_gene_count(&self) -> usize{
		self.genes.len()
	}

	pub fn get_max_gene_id(&self) -> usize{
		self.get_gene_count() + self.offset
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
	
	pub fn print( &self ){
		println!("{self}");
	}

	/*
	fn change_start_id ( &mut self, new_start :usize ){
		if self.offset == 0 {
            self.offset = new_start;
        }else {
            panic!("You try to change the start id twice - not supported!");
        }
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
	}*/
	pub fn write_index( &mut self, path: &str ) -> Result< (), String>{
		// this index needs to store the genes vectors and names. Genes binary and names as comma separated list?
		// Serialize the vector to binary
		let serialized = bincode::serialize(&self).unwrap();

	    // Write the binary data to a file
	    let mut file = File::create(path.to_string() + "/index.bin").unwrap();
	    file.write_all(&serialized).unwrap();

		// the mapper would probably also make sense to store. Let's check the time we need to re-create that before storing it.
		eprintln!("GenesMapper binary Index written to {}/index.bin", path);

		match self.write_index_txt( path ){
			Ok(_) => {},
			Err(e) => {panic!("I did not manage to write the fasta data {e}")}
		}; //will be necessary for the downstream analyses
		Ok(())
	}
	pub fn load_index(path: &str) -> Result<Self, String> {
	    let index_file_path = format!("{}/index.bin", path);
	    
	    // Attempt to open the file
	    let mut file = match File::open(&index_file_path){
	    	Ok(f) => f,
	    	Err(e) => {
	    		eprintln!("File not found {}", path);
	    		return Err("File not found".to_string() );
	    	}
	    };

	    // Read the file into a buffer
	    let mut buffer = Vec::new();
	    match file.read_to_end(&mut buffer){
	    	Ok(_) => (),
	    	Err(e) => {
	    		eprintln!("Failed to read the index file: {:?}\nYou might need to re-index the genome due to a software update?", e);
	    		return  Err("Failed to read the index file".to_string() );
	    	}
	    };

	    // Deserialize the buffer
	    let deserialized: Self = match bincode::deserialize(&buffer){
	    	Ok(this) => this,
	    	Err(e) => {
	    		eprintln!("Error deserializing the index: {:?}", e);
	    		return Err( "Error deserializing the index".to_string() )
	    	}
	    };

	    // Check if the index is empty
	    if deserialized.genes.is_empty() {
	        return Err(format!("The read index is empty {}", path));
	    }

	    eprintln!("GenesMapper binary Index loaded from {}/index.bin\n{}", path, deserialized);
	    Ok(deserialized)
	}

	/// this will simply export the mapper entries as gzipped fasta database
	pub fn write_index_txt( &mut self, path: &str ) -> Result< (), String>{

		let f1 = match File::create(path.to_string() + "/indexed_sequences.fa.gz"){
			Ok(file) => file,
			Err(err) => panic!("The file {}/indexed_sequences.fa.gz cound not be created: {err}",path )
		};
		let file1 = GzEncoder::new(f1, Compression::default());
		let mut buff1 = BufWriter::new( file1 );
		for gene in &self.genes{
			match write!(buff1, "{}", gene.to_fasta() ){
				Ok(_) => (),
				Err(err) => panic!("write_index_txt encountered an error writing the genes: {err:?}"),
			}
		}
		println!("Human readable sequences stored in {}/indexed_sequences.fa.gz", path);
		Ok(())
	}

}