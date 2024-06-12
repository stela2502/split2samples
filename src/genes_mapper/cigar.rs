//gene_mapper::cigar.rs

use crate::traits::Cell;
use crate::traits::BinaryMatcher;
use crate::genes_mapper::gene_data::GeneData;

//use crate::genes_mapper::gene_data::GeneData;
use core::cmp::Ordering;


use regex::Regex;

use core::cmp::max;
use core::fmt;

//use std::fs::File;
//use std::io::Write;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum CigarEnum{
	Match,
	Mismatch,
	Insertion,
	Deletion,
	Empty,
}

impl CigarEnum{
	pub fn opposite(&self, other: &Self ) ->bool {
		match &self{
			CigarEnum::Match => other == &CigarEnum::Mismatch,
			CigarEnum::Mismatch => other == &CigarEnum::Match,
			CigarEnum::Insertion => other == &CigarEnum::Deletion,
			CigarEnum::Deletion => other == &CigarEnum::Insertion,
			CigarEnum::Empty => panic!("You can not compare CigarEnum::Empty to anything"),
		}
	}
	pub fn to_id(&self) -> usize{
		match &self{
			CigarEnum::Match => 0,
			CigarEnum::Mismatch => 1,
			CigarEnum::Insertion => 2,
			CigarEnum::Deletion => 3,
			CigarEnum::Empty => panic!("You can not compare CigarEnum::Empty to anything"),
		}
	}
	pub fn from_str(c: &str) -> Option<CigarEnum> {
        match c {
            "I" => Some(CigarEnum::Insertion),
            "D" => Some(CigarEnum::Deletion),
            "M" => Some(CigarEnum::Match),
            "X" => Some(CigarEnum::Mismatch),
            // Add more cases as needed
            _ => None,
        }
    }
}

impl fmt::Display for CigarEnum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let direction_str = match self {
            CigarEnum::Match => "M",
			CigarEnum::Mismatch => "X",
			CigarEnum::Insertion => "I",
			CigarEnum::Deletion => "D",
			CigarEnum::Empty => panic!("There is an empty cigar entry in your vector!"),
        };
        write!(f, "{}", direction_str)
    }
}

/// This enum denotes the way a Cigar's end areas have been modified:
/// NA - no modification; Start - a longer caotic area at the start has been fixed -> end match;
/// End - a longer fix at the end of the Cigar has been fixed - a start match;
/// Both - and internal match?! Let's see if that even happens.
/// StartInsert - this cigar indicated that teh match was longer on the starting end - Likely the db too short?
///               Anyhow - this need to be fixed in the analysis scripts and therefore a CigarEndFix is needed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarEndFix{
	Na,
	Start,
	End,
	Both,
	StartInsert,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cigar{
	// the cigar string
	pub cigar:String,
	pub fixed:Option<CigarEndFix>,
	debug:bool,
	/// contains whether <CigarEnum>.to_id() is part of this cigar
	contains: Vec<bool>,
	/// state_changes is a measure of quality for a Cigar - the more changes the less likely that this is real.
	state_changes: usize,

}

// Implementing Display trait for SecondSeq
impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} - {:?}", self.cigar, self.fixed )
    }
}

impl Default for Cigar {
    fn default() -> Self {
        Cigar {
            cigar: "".to_string(),
            fixed: None,
            debug: false,
            contains: vec![false;4],
            state_changes: 1000,
        }
    }
}

impl PartialOrd for Cigar {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Cigar {
    fn cmp(&self, other: &Self) -> Ordering {
        // Define your comparison logic here. For example:
        // 1. Compare by `state_changes` (ascending)
        // 2. Then by length of `cigar` string (descending)
        // 3. Finally by `cigar` string lexicographically (ascending)

        self.mapping_quality().cmp(&other.mapping_quality()) // more is better
        .then_with( || self.len().cmp(&other.len())) // more is better
        .then_with( || other.state_changes.cmp(&self.state_changes)) // less is better
    }
}

impl Cigar{

	pub fn new( cigar: &str ) -> Self{
		Self{
			cigar: cigar.to_owned(),
			fixed:None,
			debug:false,
			contains: vec![false;4],
			state_changes: 0,
		}
	}

	pub fn better_as( &self, other: &Self ) -> bool{
		self > other
	}

	pub fn set_debug(&mut self, state:bool) {
		self.debug = state;
	}

	pub fn state_changes(&self) -> usize{
		self.state_changes
	}
	
	pub fn clear(&mut self) {
		self.cigar ="".to_string();
		self.fixed = None;
		self.contains = vec![false;4];
		self.state_changes = 0;
	}

	fn max3<T: Ord>(a: T, b: T, c: T) -> T {
        max(a, max(b, c))
    }

    pub fn to_string(&self) -> String {
    	self.cigar.to_string()
    }

    // this function is used in combination with to_sam_string to
    // fix the issues with internal insertions/dseletions in the read sequence.
    pub fn length_del_start(&self) -> usize {
		if ! self.contains[CigarEnum::Deletion.to_id()] {
			// not necessary to process this!
			return 0;
		}
		let re_start = Regex::new(r"^(\d+)D").unwrap();
		for cap in re_start.captures_iter( &self.to_string() ) {
			let count: usize = cap[1].parse().unwrap();
			return count
		}
		return 0;
    }

    pub fn to_sam_string(&self) -> String{
    	let mut ret = self.cigar.to_string();
		if ! self.contains[CigarEnum::Deletion.to_id()] {
			// not necessary to process this!
			return ret;
		}
		let re_start = Regex::new(r"^(\d+D)").unwrap();
		if let Some(mat) =re_start.captures(&ret) {
			if let Some(clippable) = mat.get(1) {
				ret.replace_range(0..clippable.end(), "")
			}
		}
		let re_end = Regex::new(r"(\d+D)$").unwrap();
		if let Some(mat) =re_end.captures(&ret) {
			if let Some(clippable) = mat.get(1) {
				ret.replace_range(clippable.start()..clippable.end(), "")
			}
		}
		return ret
    }

    pub fn fix_border_insertion( &mut self, mapping_start: usize, read:&GeneData, database:&GeneData )->usize{
    	let mut ret =0 ;
    	if self.contains[CigarEnum::Insertion.to_id()] {
			let re_start = Regex::new(r"^(\d+)I").unwrap();
			if let Some(mat) =re_start.captures(&self.cigar) {
				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
				println!("The cigar has a start I stretch - checking if the database would have the same sequence: {}",&self.cigar);

				let mut cigar_vec = self.string_to_vec( &self.cigar );
				// this means likely that my read does reach into the < start area of the gene
				let count: usize = mat[1].parse().unwrap();

				for i in 0..count  { //  runs once for count==1
					if read.get_nucleotide_2bit( i ) == database.get_nucleotide_2bit( mapping_start - (count -i) ){
						#[cfg(all(debug_assertions, feature = "mapping_debug"))]
						println!("The sequence at database position {} is the same as the sequence on read position {}", 
							mapping_start - (count -i), i );
						cigar_vec[i] = CigarEnum::Match;
					}else {
						#[cfg(all(debug_assertions, feature = "mapping_debug"))]
						println!("The sequence at database position {} is NOT the same as the sequence on read position {}",
							 mapping_start - (count -i),
							 i
						);
						cigar_vec[i] = CigarEnum::Mismatch;
					}
				}
				ret = count;
				self.clear();
				self.convert_to_cigar( &cigar_vec );
				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
				println!("The updated Cigar looks like that: {self}");
			}
			// we now also need to check if the end would also contain I's
			let re_end = Regex::new(r"(\d+)I$").unwrap();
			if let Some(mat) =re_end.captures(&self.cigar) {

				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
				println!("The cigar has a end I stretch - checking if the database would have the same sequence: {}",&self.cigar);

				let mut cigar_vec = self.string_to_vec( &self.cigar );
				// this means likely that my read does reach into the < start area of the gene
				let count: usize = mat[1].parse().unwrap();
				let (mine, other) = self.calculate_covered_nucleotides(&self.cigar );
				let cigar_vec_len = cigar_vec.len();
				for i in 0..count  { //  runs once for count==1
					if read.get_nucleotide_2bit( read.len() - count + i ) == database.get_nucleotide_2bit( mapping_start + other +i  ) {
						#[cfg(all(debug_assertions, feature = "mapping_debug"))]
						println!("The sequence at database position {} is the same as the sequence on read position {}",
							mapping_start + other +i  , read.len() - count + i +1
						);
						cigar_vec[ cigar_vec_len - i -1 ] = CigarEnum::Match;
					}else {
						#[cfg(all(debug_assertions, feature = "mapping_debug"))]
						println!("The sequence at database position {} ({:?}) is NOT the same as the sequence on read position {} ({:?})", 
							mapping_start + other +i +1, database.get_nucleotide_2bit( mapping_start + other +i +1  ),
							read.len() - count + i , read.get_nucleotide_2bit( read.len() - count + i +1 )
						);
						cigar_vec[ cigar_vec_len - i -1 ] = CigarEnum::Mismatch;
					}
				}
				self.clear();
				self.convert_to_cigar( &cigar_vec );
				#[cfg(all(debug_assertions, feature = "mapping_debug"))]
				println!("The updated Cigar looks like that: {self}");
			}
		}
		ret
    }

    fn string_to_vec(&self, cigar_string:&str ) -> Vec<CigarEnum> {
    	let mut result = Vec::new();
    	let re = Regex::new(r"(\d+)(\w)").unwrap();
    	for cap in re.captures_iter(cigar_string) {
	        let count: usize = cap[1].parse().unwrap();
	        let operation = &cap[2];
	        let cigar_type = match CigarEnum::from_str(operation) {
	        	None => { 
	        		panic!("You try to revert the Cigar String {cigar_string} back to a vector of CigarEnum - but we do not recognize the {operation} as a CigarEnum!")
	        	},
	            Some (cig) => cig,
	        };

	        // Add the operation to the result vector `count` times
	        for _ in 0..count {
	            result.push(cigar_type);
	        }
   		}
   		result
    }

    pub fn restart_from_cigar(&mut self, cigar_string:&str ) {
    	self.clear();

    	let result = self.string_to_vec(cigar_string);

	    self.convert_to_cigar( &result );
    }


    /// cleanup in this regards is meant to fix alignment errors.
    pub fn clean_up_cigar<T>(&mut self, _seq1:&T, _seq2:&T)
    where
    	T:BinaryMatcher + std::fmt::Display,
    {
    	if self.fixed.is_some() {
    		//eprintln!("This cigar was already fixed!");
   			return ();
    	}

        //println!("Before the soft clip I have this result {self} for these sequences:{seq1}\n {seq2}\n");
        self.soft_clip_start_end();

	    //println!("After the soft clip I have this result {self} for these sequences:{seq1}\n {seq2}\n");
    }

    pub fn mapping_quality(&self ) -> u8 {

    	if &self.cigar == ""{
    		return 0
    	}
    	let re = Regex::new(r"(\d+)(\w)").unwrap();
    
	    // Initialize counts for 'M' and other operations
	    let mut m_count = 0;
	    let mut other_count = 0;
	    
	    // Iterate through matches of the regular expression
	    for cap in re.captures_iter(&self.cigar) {
	        let count: usize = cap[1].parse().unwrap();
	        let operation = &cap[2];
	        match operation {
	            "M" => m_count += count,
	            //"S" => (),// soft clipped is ignored here
	            //"H" => (),// hard clipped is ignored here
	            "N" => (),// not matched (intron) is ignored, too
	            _ => other_count += count,
	        }
	    }
	    if other_count == 0 {
	    	return 40_u8
	    }
	    let ratio = (m_count as f32) / ((other_count + m_count) as f32 );
	    (40.0  * ratio ) as u8
    }

    pub fn edit_distance(&self ) -> f32 {
    	let re = Regex::new(r"(\d+)(\w)").unwrap();
    
	    // Initialize counts for 'M' and other operations
	    let mut other_count = 0.0;
	    let mut total = 0.0;	    
	    // Iterate through matches of the regular expression
	    for cap in re.captures_iter(&self.cigar) {
	        let count: f32 = cap[1].parse().unwrap();
	        total += count;
	        let operation = &cap[2];
	        match operation {
	            "M" => (), // match
	            //"S" => total -= count,// soft clipped is igniored here
	            //"H" => total -= count,// hard clipped is igniored here
	            "N" => total -= count,// not matched (intron) also ignored
	            "=" => (), // not quite sure, but should likely be good too - or?
	            _ => other_count += count,
	        }
	    }
	    other_count  / total
    }

    pub fn mapped(&self) -> f32{
    	let re = Regex::new(r"(\d+)(\w)").unwrap();
    	let mut total = 0.0;
    	for cap in re.captures_iter(&self.cigar) {
    		let count: f32 = cap[1].parse().unwrap();
    		match &cap[2]{
    			"M" => total += count,
    			_ => (),//only count MATCH
    		}
    	}
    	total
    }


    /// This will soft clip (or better replace with X's) unstable alignements that shift between states frequently and have a low matching count
    /// These areas are identified by those RegExp constructs:
    /// let start = r"^((?:[1-7]M|[1-9][0-9]*[IXD]){4,})";
    /// let end = r"((?:[1-7]M|[1-9][0-9]*[IXD]){4,})$";
    pub fn soft_clip_start_end( &mut self) {
    	if self.state_changes < 6 {
    		self.fixed = Some(CigarEndFix::Na);
    		return;
    	}
    	let start = r"^((?:[1-7]M|[1-9][0-9]*[IXD]){4,})";
    	//let start = r"^((?:[123456789][IXMD]){5,})";
    	let re_start = regex::Regex::new(start).unwrap();

		let end = r"((?:[1-7]M|[1-9][0-9]*[IXD]){4,})$";
    	//let end = r"((?:[123456789][IXMD]){5,})$";
    	let re_end = regex::Regex::new(end).unwrap();

    	let old_cigar = &self.cigar.clone();
		self.fixed = Some(CigarEndFix::Na);

	    if let Some(mat) = re_end.captures(&old_cigar) {
	        if let Some(clippable) = mat.get(1) {
	        	if clippable.start() > 0{
	        		let prev_char:&char = &self.cigar[(clippable.start()-1)..clippable.start()].chars().next().unwrap();
				    if prev_char.is_digit(10) {
				        // The character is a digit - overmatched
				        let clipped_part = &self.cigar[(clippable.start()+2)..clippable.end()];

						#[cfg(debug_assertions)]	        
				        let (mine, other) = self.calculate_covered_nucleotides(clipped_part);
				        #[cfg(not(debug_assertions))]
				        let (mine, _other) = self.calculate_covered_nucleotides(clipped_part);
				        #[cfg(debug_assertions)]
				        {
				        	println!("Detected an over match!");
				        	println!("	clipped part = {clipped_part} with a mine size of {mine} and an other size of {other}");
				        }
				        self.cigar.replace_range((clippable.start()+2)..clippable.end(), &format!("{}X", mine));
				        self.fixed = Some(CigarEndFix::End);
				        #[cfg(debug_assertions)]
				        {
				        	println!("updated cigar!\n{self}");
				        }
				    } else {
				        // The character is not a digit - great
				        let clipped_part = &self.cigar[clippable.start()..clippable.end()];
		        		let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
		        		//println!("I'll clip this part: {clipped_part} with a mine length of {mine}");
		            	self.cigar.replace_range(clippable.start()..clippable.end(), &format!("{}X", mine));
		            	self.fixed = Some(CigarEndFix::End);
				    }
					
				}else {
			        // The character is not a digit - great
			        let clipped_part = &self.cigar[clippable.start()..clippable.end()];
			        let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
			        //println!("not an digit - I'll clip this part: {clipped_part} with a mine length of {mine}");
			        self.cigar.replace_range(clippable.start()..clippable.end(), &format!("{}X", mine));
			        self.fixed = Some(CigarEndFix::End);
			    }
	        }
	    }

	    let problem = r"^\d\d*S$";
		let re_problem = regex::Regex::new(problem).unwrap();

	   	if let Some(_mat) = re_problem.captures(&self.cigar.clone()) {
	   		//panic!("This should not happen in the tests!");
	   		self.state_changes = 1;
	   		return;
	   	}
		#[cfg(debug_assertions)]
	   	println!("This is the Cigar after the end fix: {self}");
	    
	    if let Some(mat) = re_start.captures(&old_cigar) {
	        if let Some(clippable) = mat.get(1) {
	            if clippable.start() <= clippable.end() && clippable.end() <= self.cigar.len() {
	        		let clipped_part = &self.cigar[clippable.start()..clippable.end()];
	        		let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
	            	self.cigar.replace_range(clippable.start()..clippable.end(), &format!("{}X", mine));
	            	match self.fixed{
	            		Some(CigarEndFix::End) => self.fixed = Some( CigarEndFix::Both),
	            		Some(CigarEndFix::Na) => self.fixed =Some( CigarEndFix::Start),
	            		_ => unreachable!() ,
	            	}
	        	}else {
	        		// This can be totally normal for really really crappy matches.
	        		// panic!("With {self} I found a crappy match {} {}: {} old: {}", clippable.start(), clippable.end(), self.cigar.len(), old_cigar);
	        	}
	            
	        }
	    }

	    self.state_changes = 0;
	    for c in self.cigar.chars() {
	        if ! c.is_digit(10) {
	        	self.state_changes +=1;
	        }
	    }


    }

    /// Soft clip the sequences at start or end of a read that map with a horrible mapping like 1X1M2X1M4X3M1X1D1X1M3X1M2X1M2X1M1X1I3X59M,
    /// instead of that string I would like to see a 23S59M 

	/// calculates the nucleotides on both mine and the other sequence that has passed at the end of the cigar string
	pub fn calculate_covered_nucleotides(&self, cigar_string: &str) -> (usize, usize) {
	    let mut mine = 0;
	    let mut other = 0;
	    let mut current_number = String::new();
	    let mut inserts = 0;
	    let mut deletions = 0;
	    
	    for c in cigar_string.chars() {
	        if c.is_digit(10) {
	            // If the character is a digit, append it to the current number
	            current_number.push(c);
	        } else {
	            // If the character is not a digit, process the operation
	            let count = current_number.parse::<usize>().unwrap_or(1); // Parse the count, default to 1 if parsing fails

	            match c {
	                'M' | '=' | 'X' => {
	                	mine  += count; // Match, mismatch, or sequence match
	                	other += count; 
	                },
	                'I' => {
	                	inserts += count;
	                	//mine += count; // Insertion
	                },
	                'D' => {
	                	deletions += count;
	                	//other += count;// Deletion or intron
	                },
	                'S' => {
	                	mine += count;// Deletion or intron
	                },
	                'H' => {
	                	mine += count;// hard klipped
	                },
	                _ => {}, // Other CIGAR operations (e.g., P)
	            }
	            current_number.clear(); // Clear the current number for the next operation
	        }
	    }   
	    mine = mine + inserts;//.saturating_sub(deletions);
	    other = other + deletions;//.saturating_sub(inserts) ;
	    (mine, other)
	}

	pub fn len (&self) -> usize{
		let mut mine = 0;
		let mut current_number = String::new();
		let mut inserts = 0;
	    let mut deletions = 0;

		for c in self.cigar.chars() {
	        if c.is_digit(10) {
	            // If the character is a digit, append it to the current number
	            current_number.push(c);
	        } else {
	            // If the character is not a digit, process the operation
	            let count = current_number.parse::<usize>().unwrap_or(1); // Parse the count, default to 1 if parsing fails

	            match c {
	                'M' | '=' | 'X' => {
	                	mine  += count; // Match, mismatch, or sequence match
	                },
	                'I' => {
	                	inserts += count;
	                	//mine += count; // Insertion
	                },
	                'D' => {
	                	deletions += count;
	                	//other += count;// Deletion or intron
	                },
	                'S' => {
	                	mine += count;// soft klipped
	                },
	                'H' => {
	                	mine += count;// hard klipped
	                }
	                _ => {}, // Other CIGAR operations (e.g., S, H, P)
	            }
	            current_number.clear(); // Clear the current number for the next operation
	        }
	    }   
	    //mine + deletions.saturating_sub(inserts)
	    mine + inserts//.saturating_sub(deletions)
	}

	pub fn convert_to_cigar(&mut self, path: &[CigarEnum] ){
	    let mut count = 0;
	    let mut last_direction = None;

		self.cigar.clear();
		self.state_changes = 0;

		//let mut loc_path = path.to_vec();

	    for &direction in path.iter() {
	        if Some(direction) == last_direction {
	            count += 1;
	        } else {
	        	self.state_changes += 1;
	            if count > 0 {
	                self.cigar.push_str(&format!("{}{}",count, last_direction.unwrap()));
	                self.contains[last_direction.unwrap().to_id()] =  true ;
	            }
	            count = 1;
	            last_direction = Some(direction);
	        }
	    }

	    if count > 0 {
	        self.cigar.push_str(&format!("{}{}",count, last_direction.unwrap()));
	        self.contains[last_direction.unwrap().to_id()] =  true ;
	    }    
	}

	fn populate_contains( &mut self, cigar: &[CigarEnum] ) {
		for entry in cigar{
			self.contains[entry.to_id()] = true;
		}
	}

	/// this function literally checks for 1D1J or vice versa. An artifact from earlier fixes.
	pub fn fix_1d1i_1i1d(&mut self, cigar: &mut Vec<CigarEnum>, start_pos: Option<usize>) {

		if ! self.contains.iter().any(|&x| x){
			self.populate_contains(&cigar);
		}
		if ! self.contains[CigarEnum::Deletion.to_id()] &&  ! self.contains[CigarEnum::Insertion.to_id()] {
			// not necessary to process this!
			return ();
		}

	    let mut start_index = start_pos.clone().unwrap_or(cigar.len() - 1);

	    if start_index >= 3 {
	        for i in 0..=start_index - 3 {
	        	if i > cigar.len() -4 {
	        		break
	        	}
	            if ((cigar[i] == CigarEnum::Match || cigar[i] == CigarEnum::Mismatch) &&
	               (cigar[i + 1] == CigarEnum::Insertion || cigar[i + 1] == CigarEnum::Deletion) &&
	               (cigar[i + 2] == CigarEnum::Match || cigar[i + 2] == CigarEnum::Mismatch) &&
	               (cigar[i + 3] == CigarEnum::Deletion || cigar[i + 3] == CigarEnum::Insertion)) &&
	               cigar[i + 1] != cigar[i + 3] {
	                   
	                   //println!("Found a pattern - {} to {} {}{}{}{}", i + 3, i, cigar[i], cigar[i + 1], cigar[i + 2], cigar[i + 3]);
	                   
	                   // Resolve the pattern
	                   //cigar[i] = CigarEnum::Mismatch;
	                   cigar[i+2] = CigarEnum::Mismatch;
					   cigar[i + 1] = CigarEnum::Mismatch; // each DI ID element represents one X
	                   //println!("Removing {}: {}", i + 3, cigar[i + 3]);
	                   cigar.remove(i + 3);
	                   
	                   //println!("pattern changed to {}{}{}",  cigar[i], cigar[i + 1], cigar[i + 2] );
	                   
	                   // Adjust start_index
	                   start_index += 1;
	            }
	        }
	    }

	    // this sometimes can create 1D1J or 1J1D entries - which in reality are a simple X.
		let mut start_index = start_pos.unwrap_or(cigar.len() - 1);
		if start_index >= 2 {
	        for i in 0..=start_index - 2 {
	        	if i > cigar.len() -2 {
	        		break
	        	}
	        	if 
	               ((cigar[i] == CigarEnum::Insertion || cigar[i] == CigarEnum::Deletion) &&
	               (cigar[i + 1] == CigarEnum::Deletion || cigar[i + 1] == CigarEnum::Insertion)) &&
	               cigar[i] != cigar[i + 1] {
	                   
	                   //println!("Found a pattern - {} to {} {}{}", i , i+1, cigar[i], cigar[i + 1]);
	                   
	                   // Resolve the pattern
	                   cigar[i] = CigarEnum::Mismatch;

	                   //println!("Removing {}: {}", i + 1, cigar[i + 1]);
	                   cigar.remove(i + 1);
	                   
	                   //println!("pattern changed to {}",  cigar[i] );
	                   
	                   // Adjust start_index
	                   start_index += 1;
	            }
	        }
	    }
	}



	pub fn calculate_cigar(&mut self, matrix : &Vec<Vec<Cell>>, last_match:bool )  {
			// Trace back the alignment path
		
		    let mut path = Vec::with_capacity( matrix.len().max( matrix[0].len() ));

		    let rows = matrix.len();
		    let cols =  matrix[0].len();

		    let mut i = rows - 1;
		    let mut j = cols - 1;

		    // fix the final entry
		    if  last_match {
				path.push(CigarEnum::Match);
		    }else {
		    	path.push(CigarEnum::Mismatch);
		    }

		    while i > 0 && j > 0 {
		        let current_score = matrix[i][j].score;
		        //println!("Current score = {current_score}");
		        let diagonal_score = matrix[i - 1][j - 1].score;
		        let up_score = matrix[i - 1][j].score;
		        let left_score = matrix[i][j - 1].score;

		        let max = Self::max3(diagonal_score, up_score, left_score);
		        if max == diagonal_score {
		        	if diagonal_score > current_score {
		        		// mismatch!!!
		        		//println!("adding Mismatch");
		        		path.push(CigarEnum::Mismatch);
		        	}else {
		        		// match
		        		//println!("adding Match");
		        		path.push(CigarEnum::Match);
		        	}
		            i -= 1;
		            j -= 1;
		        } else if max == up_score {
		        	//println!("adding Deletion");
		        	path.push(CigarEnum::Insertion);
		            i -= 1;
		        } else {
		        	//println!("adding Insertion");
		            path.push(CigarEnum::Deletion);
		            j -= 1;
		        }
		    }

		    // If there are remaining gaps at the beginning of the sequences, fill them with the corresponding directions
		    while i > 0 {
		    	//println!("End not reached by main - adding Deletion");
		        path.push(CigarEnum::Insertion);
		        i -= 1;
		    }

		    while j > 0 {
		    	//println!("End not reached by main - adding Insertion");
		        path.push(CigarEnum::Deletion);
		        j -= 1;
		    }

		    path.reverse();

			// Convert the alignment path to CIGAR string
		    
		    //TTCATATCGACAATTAGGGTTTACGACCTCGATGTTGGATCAGGACATCCCA
		    //TTCATATCGACAATTAGGGTTTACGACCTCGATGTTTCAGGACTAGATAGTA
		    // 37M3D7M3I4M ???! 37M OK 
			self.convert_to_cigar(&path);
		    
		    /*
		    let mut csv_table = String::new();

		    for i in 0..rows {
		    	let mut line = "".to_string(); 
		        for j in 0..cols {
		        	line += &format!("{}//{}\t", matrix[i][j].score, matrix[i][j].direction );   
		        }
		        csv_table.push_str(&format!("{}\n", line));
		    }

		    // Write CSV table to a file
		    let mut file = File::create(&format!("{}.csv", &self.cigar)).unwrap();
		    file.write_all(csv_table.as_bytes()).unwrap();

		    println!("You can check my alignement table : '{}'", format!("{}.csv",  &self.cigar ) );
			*/
	}
}