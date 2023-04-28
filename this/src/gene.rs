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

    lookup
};

/// MappingInfo captures all mapping data and is a way to easily copy this data over multiple analysis runs.
pub struct Gene{
	pub chrom:String, // the cromosome id to look for the sequence
	start:usize, // the start position for this entry
	end:usize, // the end position for this entry
	exons:Vec<[usize;2]>, // a vector of start and end positions
	sense_strand:bool, // sense_strand in the genome true 1->n; false n <- 1
	name:String, // the gene symbol
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
		let mut addon="".to_string();
		if self.sense_strand{
			// assume that the first exon would be the one that we need to care about.
			// 8bp initial mapper and 32bp additional - does the exon boundary lie in that area?
			index.add( &self.to_mrna( seq.to_owned()), self.name.to_string() );
			if self.exons[0][1] - self.start > 38{
				// we could reach the intron!
				addon = "_int".to_string();
				index.add( &seq.to_owned(), self.name.to_string() + &addon );
			}
		}else {
			let compl = Self::reverse_complement( seq );
			let compl_mrna = Self::rev_compl ( self.to_mrna( seq.to_owned()));

			index.add( &compl_mrna, self.name.to_string() );
			if self.end - self.exons[self.exons.len()-1][0] > 38{
				// we could reach the intron!
				addon = "_int".to_string();
				index.add( &compl, self.name.to_string() + &addon );
			}
		}
	}

	fn to_mrna(&self, seq:Vec<u8> ) -> Vec<u8>{
		let mut size = 0;
		for reg in &self.exons{
			size += reg[1] - reg[0];
		}
		let mut mrna = Vec::<u8>::with_capacity(size);
		for reg in &self.exons{
			println!( "exon start {} and end {}", reg[0]-1, reg[1]);
			mrna.extend_from_slice(&seq[reg[0]-1..reg[1]]); 
		}
		println!(">{}\n{}\n", self.id, std::str::from_utf8( &mrna ).unwrap() );
		mrna
	}	

	pub fn passed( &self, pos:usize ) -> bool{
        self.end < pos
    }

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
	fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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
}