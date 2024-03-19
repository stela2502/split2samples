//gene_mapper::cigar.rs

use crate::traits::Cell;
use crate::traits::BinaryMatcher;
use crate::genes_mapper::gene_data::GeneData;


use regex::Regex;

use core::cmp::max;
use core::fmt;

use std::fs::File;
use std::io::Write;

#[derive(Clone, Copy, PartialEq)]
pub enum CigarEnum{
	Match,
	Mismatch,
	Insertion,
	Deletion,
}

impl fmt::Display for CigarEnum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let direction_str = match self {
            CigarEnum::Match => "M",
			CigarEnum::Mismatch => "X",
			CigarEnum::Insertion => "I",
			CigarEnum::Deletion => "D",
        };
        write!(f, "{}", direction_str)
    }
}



#[derive(Debug, Clone)]
pub struct Cigar{
	// the cigar string
	pub cigar:String,
}

// Implementing Display trait for SecondSeq
impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.cigar )
    }
}


impl Cigar{

	pub fn new( cigar: &str ) -> Self{
		Self{
			cigar: cigar.to_owned()
		}
	}

	fn max3<T: Ord>(a: T, b: T, c: T) -> T {
        max(a, max(b, c))
    }

    /// cleanup in this regards meand to fix alignment errors.
    pub fn clean_up_cigar<T>(&mut self, seq1:&T, seq2:&T)
    where
    	T:BinaryMatcher,
    {

		let missmatch_followd_by_deletion = r"(\d+)X(\d+)([ID])(\d+)M";
		let re_msbd = regex::Regex::new(missmatch_followd_by_deletion).unwrap();

		let old_cigar = &self.cigar.clone();
		let matches: Vec<regex::Captures> = re_msbd.captures_iter( &old_cigar ).collect::<Vec<regex::Captures>>();
		
		for mat in matches.into_iter().rev() {
			if let (Some(mismatch_length), Some(deletion_length), Some(direction), Some(following_match_length)) = (mat.get(1), mat.get(2), mat.get(3), mat.get(4)) {

				//panic!( "With {old_cigar} I got the end values for the three matches as {}, {} and {}", mismatch_length.end(), deletion_length.end(), following_match_length.end()); 
                let mismatch_end = mismatch_length.end();
                let deletion_end = deletion_length.end();

                let mismatch_part = &self.cigar[..(mismatch_end+1)];
           		let deletion_part = &self.cigar[..(deletion_end+1)];
           		//let match_part = &self.cigar[..(following_match_length.end()+1)];
           		let direction_str = &self.cigar[direction.start()..direction.end()];

           		let mismatch_usize = mismatch_length.as_str().parse::<usize>().unwrap();
            	let deletion_usize = deletion_length.as_str().parse::<usize>().unwrap();
            	let old_match_usize = following_match_length.as_str().parse::<usize>().unwrap();

                let ( mm_mine, mm_other) = &self.calculate_covered_nucleotides( mismatch_part );
                //let ( del_mine, del_other) = &self.calculate_covered_nucleotides( deletion_part );
                let mut new_following_match = old_match_usize;

                let mut new_mismatch = mismatch_usize;
                //panic!("Do I get the correct parts? whole: {}, missmatch {}, deletion {} and finally with the match {}",old_cigar, mismatch_part, deletion_part, match_part );

                for add in 0..mismatch_usize.min(deletion_usize){
                	//eprintln!("{self} and\nA: {seq1}\nB: {seq2}: I am trying to access the position {} in seq1 and {} in seq2", mm_mine-add, del_other-add );
                	if direction_str == "D"{
	                	if seq2.get_nucleotide_2bit(mm_mine-add-1) == seq1.get_nucleotide_2bit(mm_mine-add-1+deletion_usize){
	                		//eprintln!("And I got a match");
	                		new_mismatch -= 1;
	                		new_following_match +=1;
	                	}	
                	}else{ // "I"
                		if seq2.get_nucleotide_2bit(mm_mine-add-1+deletion_usize) == seq1.get_nucleotide_2bit(mm_mine-add-1){
	                		//eprintln!("And I got a match");
	                		new_mismatch -= 1;
	                		new_following_match +=1;
	                	}
                	}
                	
                }
                // fix the following match part
                self.cigar.replace_range((deletion_end+1)..following_match_length.end(), &format!("{}", new_following_match)) ;
                // fix the mismatch part
                if new_mismatch == 0{
                	// remove that part
                	self.cigar.replace_range(mismatch_length.start()..(mismatch_length.end()+1), "");
                }else {
                	// replace the value
                	self.cigar.replace_range(mismatch_length.start()..(mismatch_length.end()), &format!("{}", new_mismatch) );
                }

            }
        }
        self.soft_clip_start_end();
    }

    pub fn mapping_quality(&self ) -> u8 {
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
	            _ => other_count += count,
	        }
	    }
	    if other_count == 0 {
	    	return 40_u8
	    }
	    let ratio = (count as i32) / ((other_count + count) as i32 );
	    (40.0  * (1.0 - ratio)) as u8
    }

    pub fn edit_distance(&self ) -> usize {
    	let re = Regex::new(r"(\d+)(\w)").unwrap();
    
	    // Initialize counts for 'M' and other operations
	    let mut other_count = 0;
	    
	    // Iterate through matches of the regular expression
	    for cap in re.captures_iter(&self.cigar) {
	        let count: usize = cap[1].parse().unwrap();
	        let operation = &cap[2];
	        match operation {
	            "M" => (),
	            _ => other_count += count,
	        }
	    }
	    other_count 
    }

    pub fn soft_clip_start_end( &mut self) {
    	let start = r"^((?:[123456][IXMD]){5,})";
    	let re_start = regex::Regex::new(start).unwrap();

    	let end = r"((?:[123456][IXMD]){5,})$";
    	let re_end = regex::Regex::new(end).unwrap();

    	let old_cigar = &self.cigar.clone();
		
	    if let Some(mat) = re_end.captures(&old_cigar) {
	        if let Some(clippable) = mat.get(1) {
	        	if clippable.start() > 1{
	        		let prev_char:&char = &self.cigar[(clippable.start()-1)..clippable.start()].chars().next().unwrap();
				    if prev_char.is_digit(10) {
				        // The character is a digit - overmatched
				        //eprintln!("Detected an over match!");
				        let clipped_part = &self.cigar[(clippable.start()+2)..clippable.end()];
				        //eprintln!("Detected an over match! - clipped part = {clipped_part}");
				        let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
				        self.cigar.replace_range((clippable.start()+2)..clippable.end(), &format!("{}S", mine));
				    } else {
				        // The character is not a digit - great
				        let clipped_part = &self.cigar[clippable.start()..clippable.end()];
		        		let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
		        		//println!("I'll clip this part: {clipped_part} with a mine length of {mine}");
		            	self.cigar.replace_range(clippable.start()..clippable.end(), &format!("{}S", mine));
				    }
					
				}else {
			        // The character is not a digit - great
			        let clipped_part = &self.cigar[clippable.start()..clippable.end()];
			        let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
			        //println!("I'll clip this part: {clipped_part} with a mine length of {mine}");
			        self.cigar.replace_range(clippable.start()..clippable.end(), &format!("{}S", mine));
			    }
	        }
	    }

	    let problem = r"^\d\d*S$";
		let re_problem = regex::Regex::new(end).unwrap();
	   	if let Some(mat) = re_problem.captures(&self.cigar.clone()) {
	   		//panic!("This should not happen in the etsts!");
	   		return;
	   	}

	   	//println!("This is the Cigar after the end fix: {self}");
	    
	    if let Some(mat) = re_start.captures(&old_cigar) {
	        if let Some(clippable) = mat.get(1) {
	            if clippable.start() <= clippable.end() && clippable.end() <= self.cigar.len() {
	        		let clipped_part = &self.cigar[clippable.start()..clippable.end()];
	        		let (mine, _) = self.calculate_covered_nucleotides(clipped_part);
	            	self.cigar.replace_range(clippable.start()..clippable.end(), &format!("{}S", mine));
	        	}else {
	        		// This can be totally normal for really really crappy matches.
	        		// panic!("With {self} I found a crappy match {} {}: {} old: {}", clippable.start(), clippable.end(), self.cigar.len(), old_cigar);
	        	}
	            
	        }
	    }

    }

    /// Soft clip the sequences at start or end of a read that map with a horrible mapping like 1X1M2X1M4X3M1X1D1X1M3X1M2X1M2X1M1X1I3X59M,
    /// instead of that string I would like to see a 23S59M 

	/// calculates the nucleotides on both mine and the other sequence that has passed at the end of the cigar string
	fn calculate_covered_nucleotides(&self, cigar_string: &str) -> (usize, usize) {
	    let mut mine = 0;
	    let mut other = 0;
	    let mut current_number = String::new();
	    
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
	                	mine += count; // Insertion
	                },
	                'D' => {
	                	other += count;// Deletion or intron
	                }, 
	                _ => {}, // Other CIGAR operations (e.g., S, H, P)
	            }
	            current_number.clear(); // Clear the current number for the next operation
	        }
	    }   
	    (mine, other)
	}

	fn convert_to_cigar(&mut self, path: &[CigarEnum] ){
	    let mut count = 0;
	    let mut last_direction = None;

		self.cigar.clear();

	    for &direction in path.iter() {
	        if Some(direction) == last_direction {
	            count += 1;
	        } else {
	            if count > 0 {
	                self.cigar.push_str(&format!("{}{}",count, last_direction.unwrap()));
	            }
	            count = 1;
	            last_direction = Some(direction);
	        }
	    }

	    if count > 0 {
	        self.cigar.push_str(&format!("{}{}",count, last_direction.unwrap()));
	    }    
	}

	pub fn calculate_cigar(&mut self, matrix : &Vec<Vec<Cell>> )  {
			// Trace back the alignment path
		
		    let mut path = Vec::with_capacity( matrix.len().max( matrix[0].len() ));

		    let rows = matrix.len();
		    let cols =  matrix[0].len();

		    let mut i = rows - 1;
		    let mut j = cols - 1;

		    while i > 0 && j > 0 {
		        let current_score = matrix[i][j].score;
		        let diagonal_score = matrix[i - 1][j - 1].score;
		        let up_score = matrix[i - 1][j].score;
		        let left_score = matrix[i][j - 1].score;

		        let max = Self::max3(diagonal_score, up_score, left_score);
		        if max == diagonal_score {
		        	if diagonal_score > current_score {
		        		// mismatch!!!
		        		 path.push(CigarEnum::Mismatch);
		        	}else {
		        		// match
		        		path.push(CigarEnum::Match);
		        	}
		            i -= 1;
		            j -= 1;
		        } else if max == up_score {
		            path.push(CigarEnum::Deletion);
		            i -= 1;
		        } else {
		            path.push(CigarEnum::Insertion);
		            j -= 1;
		        }
		    }

		    // If there are remaining gaps at the beginning of the sequences, fill them with the corresponding directions
		    while i > 0 {
		        path.push(CigarEnum::Deletion);
		        i -= 1;
		    }

		    while j > 0 {
		        path.push(CigarEnum::Insertion);
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
		    let mut file = File::create(&format!("{}.csv", &me.cigar)).unwrap();
		    file.write_all(csv_table.as_bytes()).unwrap();
		    

		    println!("You can check my alignement table : '{}'", format!("{}.csv",  &me.cigar ) );
			*/
	}
}