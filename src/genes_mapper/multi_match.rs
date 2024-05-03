use crate::genes_mapper::MapperResult;
use crate::genes_mapper::{Cigar, CigarEndFix};


use regex::Regex;
use core::fmt;

pub struct MultiMatch{
	data:Vec<MapperResult>
}


// Implementing Display trait for MultiMatch
impl fmt::Display for MultiMatch {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		let formatted_data: String = self.data.iter().map(|x| format!("{:?}", x)).collect::<Vec<_>>().join("\n");
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

	pub fn get_best( &self, length:usize ) -> Result<MapperResult, &str> {
		//println!( "MultiMatch::get_best has been called! {self}");
		if self.data.len() == 0 {
			return Err("No Entry" );
		}
		let ret = if self.data.len() > 1{
			let mut best_mapping_qualty = 0_u8; 
			let mut best:MapperResult = MapperResult::default();

			let mut start:Option<MapperResult> = None ;
			let mut end:Option<MapperResult> = None;
			for value in &self.data{
				//println!("I am lookjing into the matching region {value}");
				let cigar = match value.cigar(){
					Some(cig) => cig.clone(),
					None => panic!("Each match to look into needs a cigar - this one has none {value}"),
				};
				if cigar.mapping_quality() > best_mapping_qualty{
					best = value.clone();
					best_mapping_qualty = cigar.mapping_quality();
					//println!("This is the best ({best_mapping_qualty}) for now: {value}");
				}
				//cigar.soft_clip_start_end();

				match cigar.fixed{
					Some(CigarEndFix::Na) | Some(CigarEndFix::Both) => {
						//them should only be checked for being best
					},
					Some(CigarEndFix::End) => {
						//println!("I found a start option: {value}");
						if start.is_some() {
							//eprintln!("Multiple possible end fixed");
							return Err("Multiple possible end fixed");
						}
						
						start = Some(
							MapperResult::new( 
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
		    				)
						);

					},
					Some(CigarEndFix::Start) => {
						if end.is_some() {
							return Err("Multiple possible start fixed");
						}
						//println!("I found a end option: {value}");
						end = Some(
							MapperResult::new( 
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
		    				)
		    			);
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
				    if start_obj.get_name() != end_obj.get_name() {
				    	eprintln!( "Potential translocation?!\nstart match\n{start_obj}\nend match\n{end_obj}\n" );
				    	return Err("potential translocation detected");
				    }
				    let total = start_length + end_length;
				    if total == length {
				    	//println!("end_obj.start(){} - (start_obj.start() {} + start_length {})", start_obj.start(), end_obj.start(), start_length);
			    		let cigar_str = &format!("{}M{}N{}M", start_length, 
			    			start_obj.start() - (end_obj.start() + start_length), end_length );
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
		    			)
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
		else {
			Ok(self.data[0].clone())
		}; // this populates the ret variable.

		// now we add a simple quality check in the end!
		match ret{
			Ok( mapped ) => {
				let cigar = mapped.get_cigar();
				if (cigar.mapped() > length as f32 * 0.7 || cigar.mapped() > mapped.db_length() as f32 * 0.8 ) && cigar.mapping_quality() > 25{
					Ok(mapped)
				}else {
					Err("No Entry" )
				}
			},
			Err(err) => Err(err),
		}
	}


}