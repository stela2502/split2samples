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
		write!(f, "MultiMatch has {} entries: ( \n{:?} )", self.data.len(), formatted_data )
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

	pub fn get_best( &self, length:usize ) -> Result<MapperResult, String> {
		if self.data.len() == 0 {
			return Err("No Entry" );
		}
		if self.data.len() > 1{
			let mut best_mapping_qualty = 0_u8; 
			let mut best:MapperResult;

			let mut start:Option<MapperResult> = None ;
			let mut end:Option<MapperResult> = None;
			for value in self.data{
				let cigar = match value.cigar(){
					Some(cig) => Cigar::new(cig),
					None => panic!("Each match to look into needs a cigar - this one has none {value}"),
				};
				match cigar.fixed{
					Some(CigarEndFix::Na) | Some(CigarEndFix::Both) => {
						if cigar.mapping_quality() > best_mapping_qualty{
							best = value;
							best_mapping_qualty = cigar.mapping_quality();
						}
					},
					Some(CigarEndFix::End) => {
						if start.is_some() {
							return Err("Multiple possible end fixed".to_string());
						}
						if cigar.mapping_quality() > best_mapping_qualty{
							best = value.clone();
							best_mapping_qualty = cigar.mapping_quality();
						}
						start = Some(value);

					},
					Some(CigarEndFix::Start) => {
						if end.is_some() {
							return Err("Multiple possible start fixed".to_string());
						}
						if cigar.mapping_quality() > best_mapping_qualty{
							best = value.clone();
							best_mapping_qualty = cigar.mapping_quality();
						}

						end = Some(value);
					},
					None=> panic!("In order to identify the best match I need the Cigar information!")
				}
			} // for loop ends here
			
			if start.is_some() && end.is_some() {
				// Define regex patterns
				let start_length: Option<usize> = None;
				let end_length: Option<usize> = None;

				let start_cigar = match start.cigar(){
					Some(cig) =>  cig,
					None=> panic!("The start match has no Cigar element!? Impossible!"),
				};
				let end_cigar = match end.cigar(){
					Some(cig) =>  cig,
					None=> panic!("The end match has no Cigar element!? Impossible!"),
				};

				let initial_matching_pattern = Regex::new(r"^(\d+)M").unwrap();
				let trailing_matching_pattern = Regex::new(r"(\d+)M$").unwrap();
				
				// Match and extract initial matching area
				if let Some(caps) = initial_matching_pattern.captures(start_cigar) {
					start_length = Some( caps[1].parse().unwrap() );
					println!("Initial matching area length in start_cigar: {:?}", start_length.ok_or(0));

				} else {
					println!("No initial matching area found in start_cigar");
				}

			    // Match and extract trailing matching area
			    if let Some(caps) = trailing_matching_pattern.captures(end_cigar) {
			    	end_length = Some( caps[1].parse().unwrap() );
			    	println!("Trailing matching area length in end_cigar: {:?}", end_length);
			    } else {
			    	println!("No trailing matching area found in end_cigar");
			    }
			    if start.gen_name() != end.get_name() {
			    	eprintln!( "Potential translocation?!\nstart match\n{start}\nend match\n{end}\n" );
			    	return Err("potential translocation detected".to_string());
			    }
			    if start_length.is_some() && end_length.is_some() {
			    	let total = start_length.ok_or(0) + end_length.ok_or(0);
			    	if total == length {
			    		println!("The length of start match + end match fits the expected length" );
			    		let cigar = Cigar::new( &format!("{}M{}N{}M", start_length.ok_or(0), 
			    			end.start - (start.start + start_length.ok_or(0)), end_length.ok_or(0) ));
			    		return Ok(
			    			MapperResult::new( 
			    				start.ok_or()?.gene_id(), 
			    				start.ok_or()?.start(), 
			    				true, 
			    				Some(format!("{}",cigar)), 
			    				cigar.mapping_quality(), 
			    				length, 
			    				cigar.edit_distance(), 
			    				start.ok_or()?.get_name() 
			    				)
			    			)
			    	}else {
			    		panic!("There seams to also be some mutations in these two matches?:\n{start:?}\n{end:?}\n")
			    	}
			    }else {
			    	// seems we only get a start or an end or none of them
			    	return Ok(best)
			    }
			}
		}
		
		else {
			Ok(self.data[0])
		}
	}


}