use crate::fast_mapper::FastMapper;
use crate::mapper_entries::MapperEntry;
use std::collections::HashSet;

use crate::gene::Gene;

use regex::Regex;

/// This gene family represents exactly that. Multiple gene entries that we want as much info on as possible.
/// We are looking for gene specific mappers, family specific mappers and possibly also something in between 
pub struct GeneFamily{
	pub chrom:Vec<String>, // the cromosomes the family entries lie on
	pub family: Vec<Gene>, // the family entries
	pub name:String, // the gene symbol
}

impl GeneFamily{
	pub fn new( gene:Gene, name: String ) -> Self {
		let family = Vec::<Gene>::with_capacity( 100 ); 
		let mut chrom = Vec::<String>::with_capacity( 5 ); 
		chrom.push( gene.chrom.clone() );
		Self{
			chrom,
			family,
			name,
		}
	}
	/// index this gene family - This function will work differently if you give it a filled in index or an empty one.
	/// The empty index will be filled with family specific (highest overlap) or parially family mapping ( less than max_per_mapper)
	/// different family mambers have this sequence or member specififc with one family member specific sequence combinatrions.
	/// Depending on the class of the mapping the gene name will either be family.name::member::name for the singlets, family.name::n
	/// for all multi mappers. This can of cause add some issues downstream, but that is a problem for another day.
	pub fn index( &self, index:&mut FastMapper, max_area:usize, seq_records:&HashSet<String> , max_per_mapper:usize ) {
		let mut loc_idx = FastMapper::new( index.kmer_len.clone() ); // use the other mappers kmer length
		let mut max_area_loc = max_area;

		let chr = Regex::new(r"^chr").unwrap();

		for gene in &self.family{
			if max_area == 0{
				max_area_loc = gene.end - gene.start;
			}
			match seq_records.get( &gene.chrom.to_string() ){
	            Some(seq) => {
	                gene.add_to_index( seq.as_bytes(), &mut loc_idx, max_area_loc );
	                //println!("The genes detected: {:?}", index.names_store );
	            },
	            None => {
	                if chr.is_match ( &gene.chrom.to_string() ){
	                    match seq_records.get( &gene.chrom.to_string()[3..] ){
	                        Some(seq) => {
	                            gene.add_to_index( seq.as_bytes(), &mut loc_idx, max_area_loc );
	                            //println!("The genes detected: {:?}", index.names_store );
	                        },
	                        None => {
	                            //missing_chr.insert( gene.chrom.to_string() );
	                            eprintln!("I do not have the sequence for the chromosome {}", gene.chrom.to_string() );
	                        }
	                    }
	                }else {
	                    match seq_records.get( &format!("chr{}", &gene.chrom.to_string()) ){
	                        Some(seq) => {
	                            gene.add_to_index( seq.as_bytes(), &mut loc_idx, max_area_loc );
	                            //println!("The genes detected: {:?}", index.names_store );
	                        },
	                        None => {
	                            //missing_chr.insert( gene.chrom.to_string() );
	                            eprintln!("I do not have the sequence for the chromosome {}", gene.chrom.to_string() );
	                        }
	                    }
	                }
	            }      
	        }
    	}


    	// Now we have an loc_idx that describes the gene family. We should now be able to identify the entry sepecific
    	// unique sequences
    	let mut entry_specific = 0;
    	let mut id = 0;
    	for mapper_entries in &loc_idx.mapper { // mapper entries is a MapperEntry object
    		if mapper_entries.has_data(){
	    		for ( secondary_key, name_entry ) in &mapper_entries.map{
	    			if name_entry.data.len() == 1 {
	    				eprintln!("I go a unique match for the famliy {} and the entry {}", self.name, loc_idx.names_store[name_entry.data[0]] );
		    			match index.incorporate_match_combo( id, name_entry.data[0], loc_idx.names_store[name_entry.data[0]].to_string() ){
		    				Err(e) => {println!("index.incorporate_match_combo hit a wall: {e:?}")}
		    				_ => (),
		    			}
		    		}
		    		else if name_entry.data.len() <= max_per_mapper{
		    			let this_name = format!("{}::{}", self.name, name_entry.data.len() );
		    			eprintln!("I go a unique low family member mapper for the famliy {} and the entry {}", self.name, this_name );
		    			match index.incorporate_match_combo( id, name_entry.data[0], this_name ){
		    				Err(e) => {println!("index.incorporate_match_combo hit a wall: {e:?}")}
		    				_ => (),
		    			}
		    		}
		    		else {
		    			let this_name = format!("{}::multi", self.name );
		    			eprintln!("I go a unique low family member mapper for the famliy {} and the entry {}", self.name, this_name );
		    			match index.incorporate_match_combo( id, name_entry.data[0], this_name ){
		    				Err(e) => {println!("index.incorporate_match_combo hit a wall: {e:?}")}
		    				_ => (),
		    			}
		    		}
	    		}
	    	}
	    	id +=1;
    	}

    	// And now get rid of all the BAD mappers - not sure of how to do that.
    	// Lets check how this works for now.
	}
}


