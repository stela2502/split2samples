use crate::fast_mapper::FastMapper;
use std::collections::HashMap;


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
	pub fn new( name: String, gene:Gene ) -> Self {
		let family = Vec::<Gene>::with_capacity( 100 ); 
		let mut chrom = Vec::<String>::with_capacity( 5 ); 
		chrom.push( gene.chrom.clone() );
		Self{
			chrom,
			family,
			name,
		}
	}

	pub fn push( &mut self, gene: Gene ) {
		let mut alread_in = false;
		for chr in &self.chrom{
			if chr == &gene.chrom {
				alread_in = true;
				break;
			}
		}
		if ! alread_in {
			self.chrom.push(  gene.chrom.to_string() );
		}
		self.family.push( gene );
	}

	/// index this gene family - This function will simply add all info into the index 
	/// Use the make_index_te_ready function
	/// I need a mut gene_id_storage where the gene information that is needed later on
	/// to 'qualify' the mapping evens. These contain all info I have about a gene.
	pub fn index( &self, index:&mut FastMapper, max_area:usize, seq_records:&HashMap<String, Vec<u8>> ) {
		let mut max_area_loc = max_area;

		let chr = Regex::new(r"^chr").unwrap();

		for gene in &self.family{
			if max_area == 0{
				max_area_loc = gene.end - gene.start;
			}
			match seq_records.get( &gene.chrom.to_string() ){
	            Some(seq) => {
	            	//println!( "Trying to add a gene to the index with a total seq length of {} and name {}", seq.len(), gene.name );
	            	//println!( "gene start {} and end {}", gene.start, gene.end );
	                gene.add_to_index( seq, index, max_area_loc, false );
	                //println!("The genes detected: {:?}", index.names_store );
	            },
	            None => {
	                if chr.is_match ( &gene.chrom.to_string() ){
	                    match seq_records.get( &gene.chrom.to_string()[3..] ){
	                        Some(seq) => {
	                            gene.add_to_index( seq, index, max_area_loc, false );
	                            //println!("The genes detected: {:?}", index.names_store );
	                        },
	                        None => {
	                            //missing_chr.insert( gene.chrom.to_string() );
	                            eprintln!("I do not have the sequence for the chromosome {}", gene.chrom );
	                        }
	                    }
	                }else {
	                    match seq_records.get( &format!("chr{}", &gene.chrom.to_string()) ){
	                        Some(_seq) => {
	                            //let gene_id = gene.add_to_index( seq, index, max_area_loc);
	                      
	                            //println!("The genes detected: {:?}", index.names_store );
	                        },
	                        None => {
	                            //missing_chr.insert( gene.chrom.to_string() );
	                            eprintln!("I do not have the sequence for the chromosome {}", gene.chrom );
	                        }
	                    }
	                }
	            }      
	        }
    	}

	}
}


