use crate::genes_mapper::MapperResult;
use crate::genes_mapper::{Cigar, CigarEndFix};

use rand::Rng;

use regex::Regex;
use core::fmt;

pub struct MultiMatch{
	data:Vec<MapperResult>
}


// Implementing Display trait for MultiMatch
impl fmt::Display for MultiMatch {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		let formatted_data: String = self.data.iter().map(|x| format!("{:?}, len: {}", x, x.get_cigar().len())).collect::<Vec<_>>().join("\n");
		write!(f, "MultiMatch has {} entries: ( \n{} )", self.data.len(), formatted_data )
	}
}



impl MultiMatch{
	pub fn new() -> Self{
		Self{
			data: Vec::<MapperResult>::with_capacity(4),
		}
	}
	pub fn push( &mut self, data:MapperResult) {
		self.data.push( data );
	}
	pub fn clear( &mut self ) {
		self.data.clear();
	}
	pub fn len( &self ) -> usize{
		self.data.len()
	}

	/// Checks if both MapperResults describe the same gene, but one maps to the unspliced and one to the spliced.
    fn get_spliced_mapper(&self, a: &MapperResult, b: &MapperResult) -> Option<MapperResult> {
        let nascent = Regex::new(r"(.*)_int$").ok()?;
        // that one matches to all gene_ids like GM12345 or EMBL000021312312 or 123123123Rik
        let contains_digit_range = Regex::new(r"\d{4}").ok()?;
        
        let a_name = a.get_name();
        let b_name = b.get_name();

        if a_name == b_name {
        	// ok that is strange - we got two matches to the same gene with the same "quality"
        	#[cfg(debug_assertions)]
        	eprintln!("We have two different matches to the same gene with the same quality!?!\n{a}\n{b}");
        	return Some(a.clone())
        }
        #[cfg(debug_assertions)]
        eprintln!("prcessing the two names {a_name} and {b_name}");
        
        let a_core = nascent.captures(a_name).and_then(|caps| caps.get(1).map(|m| m.as_str()));
        let b_core = nascent.captures(b_name).and_then(|caps| caps.get(1).map(|m| m.as_str()));
        
        match (a_core, b_core) {
            (Some(a_core_name), _) if a_core_name == b_name => Some(b.clone()),
            (_, Some(b_core_name)) if b_core_name == a_name => Some(a.clone()),
            _ => {
            	// probably mark the fact this is a multimapper for later (bam file).
            	// for now - take one at random.
            	if a_core.is_some() && b_core.is_none() {
            		Some(b.clone())
            	}else if a_core.is_none() && b_core.is_some(){
            		Some(a.clone())
            	}else {
					// Check if only one of the names ends with a stretch of 4 numbers
		            if contains_digit_range.is_match(a_name) && !contains_digit_range.is_match(b_name) {
		                Some(b.clone())
		            } else if !contains_digit_range.is_match(a_name) && contains_digit_range.is_match(b_name) {
		                Some(a.clone())
		            } else {
		                match (a_core.is_some(), b_core.is_some()) {
		                    (true, false) => Some(b.clone()),
		                    (false, true) => Some(a.clone()),
		                    _ =>  {
			                	// Compare positions on the RNA - closer to end is better
			                	if a.position_from_end() < b.position_from_end() && a.position_from_end() < 1200 {
			                		Some(a.clone())
			                	}else if b.position_from_end() < a.position_from_end() && b.position_from_end() < 1200 {
			                		Some(b.clone())
			                	}else {
			                		match a.position_from_end() < b.position_from_end() {
			                			true => Some(a.clone()),
			                			false => Some(b.clone())
			                		}
			                	}
			                }
			            }
			        }
            	}
            },
        }
    }

	pub fn get_best( &self, length:usize ) -> Result<MapperResult, &str> {
		//println!( "MultiMatch::get_best has been called! {self}");
		if self.data.len() == 0 {
			return Err("No Entry" );
		}
		if self.data.len() > 1{
			let mut best:MapperResult = MapperResult::default();

			let mut start:Option<MapperResult> = None ;
			let mut end:Option<MapperResult> = None;
			let mut rng = rand::thread_rng();
			let mut multimapper = false;

			for value in &self.data{
				//println!("I am lookjing into the matching region {value}");
				let cigar = match value.cigar(){
					Some(cig) => cig.clone(),
					None => panic!("Each match to look into needs a cigar - this one has none {value}"),
				};
				if cigar.better_as(&best.get_cigar()) {
					best = value.clone();
					#[cfg(debug_assertions)]
					println!("This is the best for now: {value}");
				}else if cigar == best.get_cigar() {
					// handle how the two matches are compared and complain if there is a discrepancy!
					// ToDo: remove the complain part
					best = match self.get_spliced_mapper( &best, value ){
						Some(entry) =>  {
							multimapper = false;
							entry
						},
						None => {
							multimapper = true;
							eprintln!("I have two matching to two different genes with the same quality!!\nold:\n{best}new:\n{value}");
							if rng.gen_bool(0.5){
								value.clone()
							}else {
								best
							}
							
						},
					};
				}
				else{
					#[cfg(debug_assertions)]
					println!("This best value \n{best}is better than the other \n{value}");
				}
				//cigar.soft_clip_start_end();

				match cigar.fixed{
					Some(CigarEndFix::Na) | Some(CigarEndFix::Both) => {
						//them should only be checked for being best
					},
					Some(CigarEndFix::End) => {
						//println!("I found a start option: {value}");
						let this = MapperResult::new( 
			    			value.gene_id(),
			    			value.start(), 
			    			true, 
			    			Some( cigar.clone() ), 
			    			cigar.mapping_quality(), 
			    			value.get_nw(),
			    			length, 
			    			cigar.edit_distance(), 
			    			value.get_name(), 
			    			value.db_length(),
						);
						if start.is_some() {
							start = match self.get_spliced_mapper( &start.unwrap(), &this ){
								Some(entry) => Some( entry ),
								None => return Err("Multiple possible end fixed"),
							}
						}else {
							start = Some(this)
						}
					},
					Some(CigarEndFix::Start) => {
						let this = MapperResult::new( 
			    			value.gene_id(),
			    			value.start(), 
			    			true, 
			    			Some( cigar.clone() ), 
			    			cigar.mapping_quality(), 
			    			value.get_nw(),
			    			length, 
			    			cigar.edit_distance(), 
			    			value.get_name(), 
			    			value.db_length(),
						);

						if end.is_some() {
							end = match self.get_spliced_mapper( &end.unwrap(), &this ){
								Some(entry) => Some( entry ),
								None => return Err("Multiple possible start fixed"),
							}
						}else {
							end = Some(this)
						}
					},
					Some(CigarEndFix::StartInsert) => {
						// this is handled during the creation of the sam strings
					},
					None=> panic!("In order to identify the best match I need the Cigar information! {}", cigar)
				}
			} // for loop ends here
			
			match (start.clone(), end.clone()) {
    			(Some(start_obj), Some(end_obj)) => {

    				//println!("I found both a start and an end fixed version!\n{start_obj}and\n{end_obj}\n");

					let start_cigar = match start_obj.cigar(){
						Some(cig) =>  cig.to_string(),
						None=> panic!("The start match has no Cigar element!? Impossible!"),
					};
					let end_cigar = match end_obj.cigar(){
						Some(cig) =>  cig.to_string(),
						None=> panic!("The end match has no Cigar element!? Impossible!"),
					};

					// Define regex patterns
					let start_length: usize;
					let end_length: usize;
					let initial_matching_pattern = Regex::new(r"^(\d+)M").unwrap();
					let trailing_matching_pattern = Regex::new(r"(\d+)M$").unwrap();
					
					// Match and extract initial matching area
					if let Some(caps) = initial_matching_pattern.captures(&start_cigar) {
						start_length = caps[1].parse().unwrap();
						//println!("Initial matching area length in start_cigar: {}", start_length);

					} else {
						//eprintln!("No initial matching area found in start_cigar");
						
						return Ok(best)
					}

				    // Match and extract trailing matching area
				    if let Some(caps) = trailing_matching_pattern.captures(&end_cigar) {
				    	end_length = caps[1].parse().unwrap();
				    	//println!("Trailing matching area length in end_cigar: {}", end_length);
				    } else {
				    	//eprintln!("No trailing matching area found in end_cigar");
				    	return Ok(best)
				    }

				    //println!("I have identified a start length: {} and end length: {}", start_length, end_length );
				    if start_obj.get_name() != end_obj.get_name() && ( start_obj.get_nw() < 0.1 &&  end_obj.get_nw() < 0.1 ) {
				    	eprintln!( "Potential translocation?!\nstart match\n{start_obj}\nend match\n{end_obj}\n" );
				    	return Err("potential translocation detected");
				    }
				    let total = start_length + end_length;
				    if total == length {
				    	eprintln!("end_obj.start(){} - (start_obj.start() {} + start_length {})", start_obj.start(), end_obj.start(), start_length);
				    	/*
			    		let cigar_str = &format!("{}M{}N{}M", start_length, 
			    			(end_obj.start() + start_length) - start_obj.start() , end_length );
			    		let cigar = Cigar::new( cigar_str );
			    		//println!("The length of start match + end match fits the expected length -> Cigar: {}", cigar );
			    		Ok(
			    			MapperResult::new( 
			    				start_obj.gene_id(),
			    				start_obj.start(), 
			    				true, 
			    				Some(cigar.clone()), 
			    				cigar.mapping_quality(), 
			    				start_obj.get_nw(),
			    				length, 
			    				cigar.edit_distance(), 
			    				start_obj.get_name(), 
			    				start_obj.db_length(),
		    				)
		    			)*/
		    			Ok(best)
			    	}else {
			    		//eprintln!("Start match and end match do not cover the whole read - additional mutations?!:\n{start:?}\n{end:?}\n");
				    	Ok(best)
			    	}
			    },
			    _ => {
		    		//eprintln!("either start or end is not detected - returning the best:\n{start:?}\n{end:?}\n");
			    	Ok(best)
			    }
			}
		}
		else { // we only have one entry
			Ok(self.data[0].clone())
		}
	}


}