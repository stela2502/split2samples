use crate::fast_mapper::FastMapper;
use crate::genes_mapper::GenesMapper;
use core::fmt;

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

// Implementing Display trait for Gene
impl fmt::Display for Gene {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Gene named {} with {} exons", self.name, self.exons.len() )
    }
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

	pub fn exon_length(&self) -> usize{
		self.exons.len()
	}

	

	/// Select the correct regions from the gene and underlying sequences
	/// to fill in the FastMapper index.
	/// the [u8] we get here has to be utf8 encoded!
	/// not 2bit binaries!
	pub fn add_to_index(&self, seq:&[u8], index: &mut FastMapper, covered_area:usize, print: bool ){

		let ( m_rna, raw) = self.generate_rna_and_nascent_strings( seq, covered_area );

		if let Some(mrna) = m_rna {
			index.add( &mrna.to_owned() , self.name.to_string(), self.ids.clone() );
			if print {
				println!(">{}\n{}", self.name.to_string() + " " + &self.chrom  , std::str::from_utf8(  &mrna[ mrna.len()-covered_area.. ].to_owned()  ).unwrap() );
			}
			if let Some(nascent) = raw{
				index.add( &nascent.to_owned() , self.name.to_string() +"_int", self.ids.clone() );
				if print {
					println!(">{}\n{}", self.name.to_string() + "_int " + &self.chrom  , std::str::from_utf8(  &mrna[ mrna.len()-covered_area.. ].to_owned()  ).unwrap() );
				}
			}
		} else {
		    eprintln!("Error in gene {} {:?} - none standard nucleotides!",self.name, self.ids );
		}

		

	}

	/// Select the correct regions from the gene and underlying sequences
	/// to fill in the FastMapper index.
	/// the [u8] we get here has to be utf8 encoded!
	/// not 2bit binaries!
	pub fn add_to_genes_mapper(&self, seq:&[u8], index: &mut GenesMapper, covered_area:usize, print: bool ){

		let ( m_rna, raw) = self.generate_rna_and_nascent_strings( seq, covered_area );

		if let Some(mrna) = m_rna {
			index.add( &mrna.to_owned() , self.name.to_string(), self.chrom.to_string(), 0   );
			if print {
				println!(">{}\n{}", self.name.to_string() + " " + &self.chrom  , std::str::from_utf8(  &mrna.to_owned()  ).unwrap() );
			}
			if let Some(nascent) = raw{
				index.add( &nascent.to_owned() , self.name.to_string() +"_int", self.chrom.to_string(), 0  );
				if print {
					println!(">{}\n{}", self.name.to_string() + "_int " + &self.chrom  , std::str::from_utf8(  &nascent.to_owned()  ).unwrap() );
				}
			}
		}
		else {
		    eprintln!("Error in gene {} {:?} - none standard nucleotides!",self.name, self.ids );
		}

		
	}

	/// generate mRNA and nacent RNA (if the last exon is small enough)
	pub fn generate_rna_and_nascent_strings(&self, seq: &[u8], covered_area: usize) -> (Option<Vec<u8>>, Option<Vec<u8>>) {
	    let mrna = self.to_mrna(seq.to_owned(), covered_area);
	    let nascent = self.to_nascent(seq.to_owned(), covered_area);
	    (mrna, nascent)
	}

	/// get the mRNA sequence of the transcript in sense orientation.
	/// Fails if any other base than AGCT is in the sequence
	/// This returns the revers complement if on the opposite starnd
	pub fn to_mrna(&self, seq: Vec<u8>, covered_area:usize) -> Option<Vec<u8>> {

	    let mut mrna = Vec::<u8>::with_capacity(seq.len() ); // Allocate more space for potential additions

	    let mut sorted_exons = self.exons.clone();
	    sorted_exons.sort_by(|a, b| a[0].cmp(&b[0]));

	    // exons upper/lower case iterations to see the breaks
	    let mut lc = false;
	    for reg in &sorted_exons {
	        if reg[0] > seq.len() || reg[1] > seq.len() {
	            eprintln!("The exon positions exceed the sequence length!");
	            return None;
	        }
	        let exon_slice = &seq[reg[0] - 1..reg[1]];

	        // Convert characters to uppercase or lowercase based on `lc`
	        let exon_mrna: Vec<u8> = exon_slice
	            .iter()
	            .map(|&c| if lc { c.to_ascii_lowercase() } else { c.to_ascii_uppercase() })
	            .collect();

	        mrna.extend_from_slice(&exon_mrna);

	        lc = !lc;
	    }

	    self.cut_to_size(mrna, covered_area )

	}

	/// cut the RNA to the right size
	fn cut_to_size( &self, seq:Vec<u8>, covered_area:usize) -> Option<Vec<u8>> {
		let start = seq.len().saturating_sub(covered_area);
		if ! self.sense_strand{
			Some ( Self::rev_compl( seq)[ start.. ].to_vec() )
		}else {
			Some (  seq[ start.. ].to_vec() )
		}
	}
	
	

	/// get the nascent RNA for this transcript (including introns).
	/// Fails if any other base than AGCT is in the sequence
	/// This returns the revers complement if on the opposite starnd
	fn to_nascent(&self, seq:Vec<u8> , covered_area:usize) -> Option<Vec<u8>> {
		if self.exons.len() == 0{
			eprintln!("I have no exons?! - something is wrong here! {self}");
			return None
		}
		let last_exon = match self.sense_strand{
			true => self.exons.len()-1,
			false => 0,
		};
		if self.exons[last_exon][0] >= self.exons[last_exon][1]{
			eprintln!("What is this exon {last_exon} has a strange size: {} to {}bp?", self.exons[last_exon][0], self.exons[last_exon][1] );
			None
		}
		else if self.exons[last_exon][1]- self.exons[last_exon][0] < covered_area && self.exons.len() > 1{
			//println!("Worked! - last exon length = {} in {self}", self.exons[last_exon][1]- self.exons[last_exon][0] );
			//let size = self.end - self.start;
			let start_index = self.start.saturating_sub(1);
			let end_index = self.end.min(seq.len());
			let nascent = seq.get(start_index..end_index).unwrap_or_default().to_vec();

			self.cut_to_size(nascent, covered_area )
		}
		else {
			None
		}
	}

	/// the reverse complement of a Vec<u8>
    fn rev_compl( seq:Vec<u8> ) -> Vec<u8>{
	    seq.iter()
	        .rev()
	        .filter_map(|&c| COMPLEMENT[c as usize])
	        .collect()
    }

	/// is the position (pos) after our end?
	pub fn passed( &self, pos:usize ) -> bool{
        self.end < pos
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