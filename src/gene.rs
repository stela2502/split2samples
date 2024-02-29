use crate::fast_mapper::FastMapper;

const COMPLEMENT: [Option<u8>; 256] = {
    let mut lookup = [None; 256];

    lookup[b'A' as usize] = Some(b'T');
    lookup[b'C' as usize] = Some(b'G');
    lookup[b'G' as usize] = Some(b'C');
    lookup[b'T' as usize] = Some(b'A');
    lookup[b'a' as usize] = Some(b'T');
    lookup[b'c' as usize] = Some(b'G');
    lookup[b'g' as usize] = Some(b'C');
    lookup[b't' as usize] = Some(b'A');
    lookup[b'R' as usize] = Some(b'Y');
    lookup[b'Y' as usize] = Some(b'R');
    lookup[b'S' as usize] = Some(b'W');
    lookup[b'W' as usize] = Some(b'S');
    lookup[b'K' as usize] = Some(b'M');
    lookup[b'M' as usize] = Some(b'K');
    lookup[b'B' as usize] = Some(b'V');
    lookup[b'V' as usize] = Some(b'B');
    lookup[b'D' as usize] = Some(b'H');
    lookup[b'H' as usize] = Some(b'D');
    lookup[b'N' as usize] = Some(b'N');

    lookup
};

// const CHECK: [Option<u8>; 256] = {
//     let mut lookup = [None; 256];

//     lookup[b'A' as usize] = Some(b'A');
//     lookup[b'C' as usize] = Some(b'C');
//     lookup[b'G' as usize] = Some(b'G');
//     lookup[b'T' as usize] = Some(b'T');
//     lookup[b'a' as usize] = Some(b'A');
//     lookup[b'c' as usize] = Some(b'C');
//     lookup[b'g' as usize] = Some(b'G');
//     lookup[b't' as usize] = Some(b'T');
//     lookup
// };

/// MappingInfo captures all mapping data and is a way to easily copy this data over multiple analysis runs.
pub struct Gene{
	pub chrom:String, // the cromosome id to look for the sequence
	pub start:usize, // the start position for this entry
	pub end:usize, // the end position for this entry
	exons:Vec<[usize;2]>, // a vector of start and end positions
	sense_strand:bool, // sense_strand in the genome true 1->n; false n <- 1
	pub name:String, // the gene symbol
	pub ids:Vec<String>, // e.g. ENSMBL ID and other entries like family name or class 
}

impl Gene{
	pub fn new(chrom:String, start_s:String, end_s:String, sense_strand_s:String, name:String, ids:Vec<String> ) -> Self {
		let exons = Vec::<[usize;2]>::new();
		let start = match start_s.parse::<usize>(){
			Ok(v) => v,
			Err(e) => panic!("I could not parse the start of the transcript as usize: {e:?}"),
		};
		let end = match end_s.parse::<usize>(){
			Ok(v) => v,
			Err(e) => panic!("I could not parse the end of the transcript as usize: {e:?}"),
		};

		let sense_strand = sense_strand_s == "+";

		Self{
			chrom,
			start,
			end,
			exons,
			sense_strand,
			name,
			ids,
		}
	}
	/// Return if the exon matched to the transcript
	pub fn add_exon(&mut self, start_s:String, end_s:String ) {
		let start = match start_s.parse::<usize>(){
			Ok(v) => v,
			Err(e) => panic!("I could not parse the start of the transcript as usize: {e:?}"),
		};
		let end = match end_s.parse::<usize>(){
			Ok(v) => v,
			Err(e) => panic!("I could not parse the end of the transcript as usize: {e:?}"),
		};
		self.exons.push( [start, end] );
		self.exons.sort_by(|a, b| a[0].cmp(&b[0]));
	}

	/// Select the correct regions from the gene and underlying sequences
	/// to fill in the FastMapper index.
	/// the [u8] we get here has to be utf8 encoded!
	/// not 2bit binaries!
	pub fn add_to_index(&self, seq:&[u8], index: &mut FastMapper, covered_area:usize, print: bool ){


		if let Some(mrna) = self.to_mrna(seq.to_owned()) {
		    // The mrna sequence will be unchecked for non standard nucleotides - this needs to be checked later
		    // But it will always be in sense orientation - now.
		    
		    // Stop this behaviour for now - not really necessary any more
		    //println!(">{}\n{}", self.name.to_string() + " total mRNA -> chr " + &self.chrom  , std::str::from_utf8( &mrna ).unwrap() );

		    if mrna.len() > covered_area{
				//eprintln!( "adding this mrna to the index: \n{} -> \n{}", self.name.to_string() , std::str::from_utf8(&mrna[ 0..100]).expect("Invalid UTF-8") );
				index.add( &mrna[ mrna.len()-covered_area.. ].to_owned() , self.name.to_string(), self.ids.clone() );
				if print {
					println!(">{}\n{}", self.name.to_string() + " " + &self.chrom  , std::str::from_utf8(  &mrna[ mrna.len()-covered_area.. ].to_owned()  ).unwrap() );
				}
			}
			else {
				//eprintln!( "adding this mrna to the index: \n{} -> \n{}", self.name.to_string() , std::str::from_utf8(&mrna).expect("Invalid UTF-8") );
				index.add( &mrna , self.name.to_string(), self.ids.clone() );
				if print {
					println!(">{}\n{}", self.name.to_string() + " " + &self.chrom  , std::str::from_utf8( &mrna ).unwrap() );
				}
			}
		} else {
		    eprintln!("Error in gene {} {:?} - none standard nucleotides!",self.name, self.ids );
		}
		let last_exon = match self.sense_strand{
			true => self.exons.len()-1,
			false => 0,
		};
		if self.exons[last_exon][1]- self.exons[last_exon][0] > 100 && self.exons.len() > 1{
			let addon = "_int".to_string();
			match &self.to_nascent( seq.to_owned()){
				Some( nascent ) => {
					if nascent.len() > covered_area{
						//eprintln!( "adding this nascent to the index: \n{} -> \n{}", self.name.to_string() + &addon, std::str::from_utf8(&nascent[ 0..100]).expect("Invalid UTF-8") );
						index.add( &nascent[ nascent.len()-covered_area..].to_owned() , self.name.to_string() + &addon, self.ids.clone() );
						if print {
							println!(">{}\n{}", self.name.to_string()+ &addon + " " + &self.chrom  , std::str::from_utf8(  &nascent[ nascent.len()-covered_area.. ].to_owned()  ).unwrap() );
						}
					}else{
						//eprintln!( "adding this nascent to the index: \n{} -> \n{}", self.name.to_string() + &addon, std::str::from_utf8(&nascent).expect("Invalid UTF-8") );
						index.add( nascent , self.name.to_string() + &addon, self.ids.clone()  );
						if print {
							println!(">{}\n{}", self.name.to_string()+ &addon + " " + &self.chrom  , std::str::from_utf8(  &nascent.to_owned()  ).unwrap() );
						}
					}
				},
				None=> {
					eprintln!("Error in gene {} {:?} - none standard nucleotides!",self.name, self.ids );
				}
			}
		}
	}

	/// get the mRNA sequence of the transcript in sense orientation.
	/// Fails if any other base than AGCT is in the sequence
	/// This returns the revers complement if on the opposite starnd
	fn to_mrna(&self, seq:Vec<u8> ) -> Option<Vec<u8>>{
		
		let size = self.exons.iter().map(|reg| reg[1] - reg[0] + 1).sum();
		let mut mrna = Vec::<u8>::with_capacity(size);

		let mut sorted_exons = self.exons.clone();

		sorted_exons.sort_by(|a, b| a[0].cmp(&b[0]));
		let mut lc = false;
		for reg in &sorted_exons{
			if reg[0] > seq.len() || reg[1] > seq.len() {
				eprintln!("The exon positions exeed the seq length!");
				return None;
	        }
	        if lc {
	        	let inverted_slice: Vec<u8> = seq[reg[0] - 1..reg[1]]
				    .iter()
				    .map(|&c| {
				    	c.to_ascii_lowercase()
				    })
				    .collect();
				mrna.extend_from_slice(&inverted_slice );
	        }else {
	        	let inverted_slice: Vec<u8> = seq[reg[0] - 1..reg[1]]
				    .iter()
				    .map(|&c| {
				    	c.to_ascii_uppercase()
				    })
				    .collect();
				mrna.extend_from_slice(&inverted_slice );

	        }
	        lc = ! lc;
			//println!( "gene {} exon start {} and end {}", self.name,reg[0]-1, reg[1]);
			//mrna.extend_from_slice(&seq[reg[0]-1..reg[1]]);
			//mrna.push(b'\n');
		}
		// // convert to 2bit
		// for &b in mrna.iter() {
    	// 	let _entr = match CHECK[b as usize] {
    	// 		Some(val) => val,
    	// 		None => return None,
		// 	};
		// }
		//println!(">{}\n{}", self.id.to_string() + " " + &self.name + " " + &self.chrom , std::str::from_utf8( &mrna ).unwrap() );
		if ! self.sense_strand{
			return Some ( Self::rev_compl( mrna ))
		}
		Some(mrna)
	}

	/// get the nascent RNA for this transcript (including introns).
	/// Fails if any other base than AGCT is in the sequence
	/// This returns the revers complement if on the opposite starnd
	fn to_nascent(&self, seq:Vec<u8> ) -> Option<Vec<u8>> {
		let size = self.end - self.start;
		let mut nascent = Vec::<u8>::with_capacity(size);
		nascent.extend_from_slice(&seq[self.start-1..self.end]);
		// for &b in nascent.iter() {
    	// 	let _entr = match CHECK[b as usize] {
    	// 		Some(val) => val,
    	// 		None => return None,
		// 	};
		// }
		if ! self.sense_strand{
			return Some ( Self::rev_compl( nascent ))
		}
		Some(nascent)
	}

	/// is the position (pos) after our end?
	pub fn passed( &self, pos:usize ) -> bool{
        self.end < pos
    }

    /// the reverse complement of a Vec<u8>
    fn rev_compl( seq:Vec<u8> ) -> Vec<u8>{
	    seq.iter()
	        .rev()
	        .filter_map(|&c| COMPLEMENT[c as usize])
	        .collect()
    }

    // /// the reverse complement of a &[u8]
	// fn reverse_complement(seq: &[u8]) -> Vec<u8> {
	//     let mut complement = Vec::with_capacity(seq.len());

	//     for &b in seq.iter().rev() {
	//     	let entr = match COMPLEMENT[b as usize] {
	//     		Some(val) => val,
	//     		None => panic!("Could not translate nucl {b}"),
	// 		};
	//         complement.push(entr);
	//     }

	//    complement
	// }
}