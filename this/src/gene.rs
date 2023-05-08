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

const CHECK: [Option<u8>; 256] = {
    let mut lookup = [None; 256];

    lookup[b'A' as usize] = Some(b'A');
    lookup[b'C' as usize] = Some(b'C');
    lookup[b'G' as usize] = Some(b'G');
    lookup[b'T' as usize] = Some(b'T');
    lookup[b'a' as usize] = Some(b'A');
    lookup[b'c' as usize] = Some(b'C');
    lookup[b'g' as usize] = Some(b'G');
    lookup[b't' as usize] = Some(b'T');
    lookup
};

/// MappingInfo captures all mapping data and is a way to easily copy this data over multiple analysis runs.
pub struct Gene{
	pub chrom:String, // the cromosome id to look for the sequence
	start:usize, // the start position for this entry
	end:usize, // the end position for this entry
	exons:Vec<[usize;2]>, // a vector of start and end positions
	sense_strand:bool, // sense_strand in the genome true 1->n; false n <- 1
	pub name:String, // the gene symbol
	pub id:String, // e.g. ENSMBL ID
}

impl Gene{
	pub fn new(chrom:String, start_s:String, end_s:String, sense_strand_s:String, name:String, id:String ) -> Self {
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
			id,
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
	}

	/// Select the correct regions from the gene and underlying sequences
	/// to fill in the FastMapper index.
	pub fn add_to_index(&self, seq:&[u8], index: &mut FastMapper ){

		if ! self.sense_strand{ // I need the reverse sequence in the index!
			// assume that the first exon would be the one that we need to care about.
			// 8bp initial mapper and 32bp additional - does the exon boundary lie in that area?
			match  &self.to_mrna( seq.to_owned()){
				Some( mrna ) => {
					if mrna.len() > 100{
						index.add( &mrna[ 0..100 ].to_owned() , self.name.to_string() );
					}
					else {
						index.add( &mrna , self.name.to_string() );
					}
				},
				None=> {
					eprintln!("Error in gene {} - none standard nucleotides!",self.name.to_string());
					return
				}
			};

			//if self.exons[ self.exons.len()-1 ][1] - self.exons[ self.exons.len()-1 ][0] < 100 {
			if self.exons[ 0 ][1] - self.exons[ 0 ][0] < 100 {
				// we could reach the intron!
				let addon = "_int".to_string();
				match  &self.to_nascent( seq.to_owned()){
					Some( nascent ) => {
						if nascent.len() > 100{
							index.add( &nascent[ 0..100].to_owned() , self.name.to_string() + &addon );
						}else{
							index.add( &nascent , self.name.to_string() + &addon );
						}
					},
					None=> {
						eprintln!("Error in gene {} - none standard nucleotides!",self.name.to_string());
						return
					}
				};
			}
		}else {
			match  self.to_mrna( seq.to_owned()){
				Some( mrna ) => {
					let compl_mrna = Self::rev_compl ( mrna );
					if compl_mrna.len() > 100{
						index.add( &compl_mrna[ 0..100 ].to_owned() , self.name.to_string() );
					}
					else {
						index.add( &compl_mrna , self.name.to_string() );
					}
				},
				None=> {
					eprintln!("Error in gene {} - none standard nucleotides!",self.name.to_string());
					return
				}
			};

			if self.exons[ self.exons.len()-1 ][1] - self.exons[ self.exons.len()-1 ][0] < 100{
				// we could reach the intron!
				let addon = "_int".to_string();
				match  self.to_nascent( seq.to_owned()){
					Some( nascent ) => {
						let compl = Self::rev_compl ( nascent );

						if compl.len() > 100{
							index.add( &compl[ 0..100 ].to_owned() , self.name.to_string() + &addon );
						}else{
							index.add( &compl , self.name.to_string() + &addon );
						}
					},
					None=> {
						eprintln!("Error in gene {} - none standard nucleotides!",self.name.to_string());
						return
					}
				};
			}
		}
	}

	/// get the mRNA sequence of the transcript in sense orientation.
	/// Fails if any other base than AGCT is in the sequence
	fn to_mrna(&self, seq:Vec<u8> ) -> Option<Vec<u8>>{
		let mut size = 0;
		for reg in &self.exons{
			size += reg[1] - reg[0];
		}
		let mut mrna = Vec::<u8>::with_capacity(size);
		for reg in &self.exons{
			println!( "gene {} exon start {} and end {}", self.name,reg[0]-1, reg[1]);
			mrna.extend_from_slice(&seq[reg[0]-1..reg[1]]); 
		}
		for &b in mrna.iter() {
    		let _entr = match CHECK[b as usize] {
    			Some(val) => val,
    			None => return None,
			};
		}
		println!(">{} antisense\n{}", self.id.to_string() + " " + &self.name + " " + &self.chrom , std::str::from_utf8( &mrna ).unwrap() );
		Some(mrna)
	}

	/// get the nascent RNA for this transcript (including introns).
	/// Fails if any other base than AGCT is in the sequence
	fn to_nascent(&self, seq:Vec<u8> ) -> Option<Vec<u8>> {
		let size = self.end - self.start;
		let mut nascent = Vec::<u8>::with_capacity(size);
		nascent.extend_from_slice(&seq[self.start-1..self.end]);
		for &b in nascent.iter() {
    		let _entr = match CHECK[b as usize] {
    			Some(val) => val,
    			None => return None,
			};
		}
		Some(nascent)
	}

	/// is the position (pos) after our end?
	pub fn passed( &self, pos:usize ) -> bool{
        self.end < pos
    }

    /// the reverse complement of a Vec<u8>
    fn rev_compl( seq:Vec<u8> ) -> Vec<u8>{
    	let mut complement = Vec::with_capacity(seq.len());

	    for &b in seq.iter().rev() {
	    	let entr = match COMPLEMENT[b as usize] {
	    		Some(val) => val,
	    		None => panic!("Could not translate nucl {b}"),
			};
	        complement.push(entr);
	    }

	    complement
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