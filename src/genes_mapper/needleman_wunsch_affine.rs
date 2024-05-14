const MATCH_SCORE: i32 = 1;
const MISMATCH_SCORE: i32 = -2;
const GAP_OPEN_PENALTY: i32 = -4; // Penalty for opening a gap
const GAP_EXTENSION_PENALTY: i32 = -3; // Penalty for extending a gap

use crate::genes_mapper::cigar::{Cigar, CigarEnum};
use crate::traits::BinaryMatcher;

use std::fs::File;
use std::io::{self,Write};

pub struct NeedlemanWunschAffine{
	dp: Vec<Vec<i32>>,
	n: usize,
	m: usize,
	cigar_vec: Option<Vec<CigarEnum>>,
	circles:usize,
	debug:bool,
}



impl <'a> NeedlemanWunschAffine {
	pub fn new( ) ->Self {
		// Initialize the DP matrix

		let me = Self {
	    	dp : Vec::<Vec::<i32>>::new(),
	    	n: 0,
	    	m: 0,
	    	cigar_vec: None,
	    	circles:0,
	    	debug: false,
	    };
	    me
	}

	pub fn debug( &self ) -> bool{
		self.debug
	}
	pub fn set_debug(&mut self, state:bool){
		self.debug = state
	}

	/// initialize does exactly that - initializes the internal storage - if necessary.
	pub fn initialize( &mut self, rows:usize, cols:usize ){

		let size = rows.max(cols);

		if self.dp.len() < size +2 || self.dp[0].len() < size +2 {
			self.dp = vec![vec![0; size + 2]; size + 2];
			for i in 1..=size {
		        self.dp[i][0] = self.dp[i - 1][0] + GAP_OPEN_PENALTY + (i - 1) as i32 * GAP_EXTENSION_PENALTY;
		        self.dp[0][i] = self.dp[0][i - 1] + GAP_OPEN_PENALTY + (i - 1) as i32 * GAP_EXTENSION_PENALTY;
		    }
		}

	}

	pub fn to_string<T>( &mut self, seq1: &T, seq2: &T, humming_cut: f32 ) -> String 
	where
    T: BinaryMatcher + std::fmt::Display{

    	let cigar_vec = self.to_cigar_vec(seq1, seq2, humming_cut);
		let ( aligned1, aligned2 ) = self.needleman_wunsch_affine_backtrack(seq1, seq2, &cigar_vec );
		let cig_str = Self::cigar_to_string( &cigar_vec );
		format!("read:\n{}\ndatabse\n{}\nalignement:\n{:?}\n{:?}\n{cig_str:?}", 
			&seq1, &seq2, 
			&aligned1 , 
			&aligned2 
		)
	}

	pub fn int_state_to_string<T>( &mut self, seq1: &T, seq2: &T, cigar_vec: &[CigarEnum] ) -> String 
	where
    T: BinaryMatcher + std::fmt::Display{
		let ( aligned1, aligned2 ) = self.needleman_wunsch_affine_backtrack(seq1, seq2, cigar_vec );
		let cig_str = Self::cigar_to_string( &cigar_vec );
		format!("alignement:\n{:?}\n{:?}\n{cig_str:?}", 
			//&seq1.as_dna_string(), &seq2.as_dna_string(), 
			&aligned1 , 
			&aligned2 
		)
	}

	/// calculates a needleman_wunsch_affine matrix and returns the final value scaled to the shorter sequence length.
	pub fn needleman_wunsch_affine<T>(&mut self, seq1: &T, seq2: &T, humming_cut: f32 ) -> f32 
	where
    T: BinaryMatcher {

	    self.n = seq1.len();
	    self.m = seq2.len();
	    let n = self.n;
	    let m = self.m;
	    if n.min(m) < 15  || seq1.tri_nuc_abs_diff(seq2) > humming_cut{
            return 100.0
        }

        self.initialize( n, m );

		// Fill in the DP matrix
		for i in 1..=n {
		    for j in 1..=m {
		        if let (Some(nuc1), Some(nuc2)) = (seq1.get_nucleotide_2bit(i - 1), seq2.get_nucleotide_2bit(j - 1)) {
		            // Calculate match/mismatch score
		            let match_mismatch_score = if nuc1 == nuc2 { MATCH_SCORE } else { MISMATCH_SCORE };

		            // Compute the maximum score directly
		            self.dp[i][j] = [
		            	// Score for diagonal (match/mismatch)
		                self.dp[i - 1][j - 1] + match_mismatch_score,
		                // Score for gap in seq1 (insertion in read)
		                self.dp[i - 1][j] + GAP_OPEN_PENALTY,
		                // Score for gap in seq2 (deletion in read)
		                self.dp[i][j - 1] + GAP_OPEN_PENALTY, 
		                // Score for gap extension in seq1 (insertion in read)
		                (1..=i - 1).map(|k| self.dp[i - k][j] + GAP_OPEN_PENALTY + k as i32 * GAP_EXTENSION_PENALTY).max().unwrap_or(std::i32::MIN),
		                // Score for gap extension in seq2 (deletion in read)
		                (1..=j - 1).map(|k| self.dp[i][j - k] + GAP_OPEN_PENALTY + k as i32 * GAP_EXTENSION_PENALTY).max().unwrap_or(std::i32::MIN),
		            ].iter().copied().max().unwrap();
		        } else {
		            // Panic if unable to retrieve nucleotide sequences
		            panic!("Could not get sequence for nuc1 or nuc2 needleman_wunsch_affine")
		        }
		    }
		}
	    let size = n.min(m) as f32;
	    self.cigar_vec = Some(self.to_cigar_vec( seq1, seq2, humming_cut ));
	    // this call might have swapped but the self.n and self.m have been updated!

	    (size - self.dp[self.n][self.m] as f32) / size

	}

	/// This will give two vectors - the aligned sequence 1 and sequence 2
	pub fn needleman_wunsch_affine_backtrack<T>(&mut self, seq1: &T, seq2: &T, cig_vec: &[CigarEnum] ) -> ( String, String ) 
	where
    T: BinaryMatcher {
		//let (mut i, mut j) = (seq1.len(), seq2.len());
		let (mut i, mut j);
	    let mut align_seq1 = Vec::<u8>::with_capacity( cig_vec.len() );
	    let mut align_seq2 = Vec::<u8>::with_capacity( cig_vec.len() );
	    i = 0;
	    j = 0;

	    for val in cig_vec {
	    	match val {
	    		CigarEnum::Insertion => {
	    			if let Some(nuc1) = seq1.get_nucleotide_2bit(i){
	    				align_seq1.push( Self::to_utf8_lc(nuc1) );
	            		align_seq2.push( b'-' );
	            		i += 1;
	    			}else {
	    				panic!("Could not get sequence for nuc1 #{i}")
	    			}
	    		},
	    		CigarEnum::Deletion => {
	    			if let Some(nuc2) = seq2.get_nucleotide_2bit(j){
	    				align_seq1.push( b'-' ); 
						align_seq2.push( Self::to_utf8_lc(nuc2) );
						j +=1;
	    			}else {
	    				eprintln!("Could not get sequence for nuc2 #{j} Something is seriousely wrong here!!!");
	    				break;
	    			}
	            	
	    		},
	    		CigarEnum::Mismatch => {
	    			if let Some(nuc1) =  seq1.get_nucleotide_2bit(i){
	    				align_seq1.push( Self::to_utf8_lc( nuc1 ) );
	    				i += 1;
	    			}else {
	    				align_seq1.push( b'-' ); 
	    			}

	    			if let Some(nuc2) =  seq2.get_nucleotide_2bit(j){
	    				align_seq2.push( Self::to_utf8_lc( nuc2 ) );
	    				j += 1;
	    			}else {
	    				align_seq2.push( b'-' ); 
	    			}
	    		},
	    		_ => {

	    			if let Some(nuc1) =  seq1.get_nucleotide_2bit(i){
	    				align_seq1.push( Self::to_utf8( nuc1 ) );
	    				i += 1;
	    			}else {
	    				align_seq1.push( b'-' ); 
	    			}
	    			if let Some(nuc2) =  seq2.get_nucleotide_2bit(j){
	    				align_seq2.push( Self::to_utf8( nuc2 ) );
	    				j += 1;
	    			}else {
	    				align_seq2.push( b'-' ); 
	    			}
	    		},
	    	}

	    }

	    ( String::from_utf8_lossy( &align_seq1).to_string(), String::from_utf8_lossy( &align_seq2).to_string() )
	}

	fn to_utf8( twobit:u8 ) -> u8 {
		match twobit{
			0 => b'A',
			1 => b'C',
			2 => b'G',
			3 => b'T',
			_ => b'N',
		}
	}

	fn to_utf8_lc( twobit:u8 ) -> u8 {
		match twobit{
			0 => b'a',
			1 => b'c',
			2 => b'g',
			3 => b't',
			_ => b'n',
		}
	}

	pub fn to_cigar_vec_inverted<T>( &mut self, seq1: &T, seq2: &T, humming_cut: f32 )-> Vec<CigarEnum> 
	where
    T: BinaryMatcher {

    	self.circles +=1;
    	if  self.circles < 2 {
    		let _ =self.needleman_wunsch_affine( seq2, seq1, humming_cut );
    	}
	    let results = self.cigar_vec();
	    self.circles = 0;

	    results.iter()
    		.map(|cigar_enum| if *cigar_enum == CigarEnum::Deletion { CigarEnum::Insertion } else if *cigar_enum == CigarEnum::Insertion { CigarEnum::Deletion } else { *cigar_enum })
    		.collect()
	}

	pub fn cigar_vec( &self ) -> Vec<CigarEnum>{
		match &self.cigar_vec{
			Some(vec) => vec.to_vec(),
			None => Vec::new(),
		}
	}


	/// This will create the Cigar states vector from the internal matrix.
	pub fn to_cigar_vec<T>(&mut self, seq1: &T, seq2: &T, _humming_cut: f32) -> Vec<CigarEnum> 
	where
    T: BinaryMatcher {

		let (mut i, mut j) = (seq1.len(), seq2.len());

	    let mut cigar = vec![CigarEnum::Empty;i.max(j)];
	    
	    let mut rev_id = cigar.len().saturating_sub(1);


	    while i > 0 || j > 0 {

	    	//if let (Some(nuc1), Some(nuc2)) = (seq1.get_nucleotide_2bit(i.saturating_sub(1)), seq2.get_nucleotide_2bit(j.saturating_sub(1)) ){
	    		let this_value = self.dp[i][j];
	    		if self.debug {
	    			println!("I am testing this part of the dp matrix with the lower right corner i:{i};j:{j}:\n{}\t{}\n{}\t{}",
		    			self.dp[i.saturating_sub(1)][j.saturating_sub(1)], self.dp[i][j.saturating_sub(1)],
		    			self.dp[i.saturating_sub(1)][j],self.dp[i][j]
		    			);
	    		}
		    	/*
		    	// this checks if the next diagonal would be a match - if iot is this is the wanted outcome - finish here!
		    	if nuc1 == nuc2 {
					cigar[rev_id] = CigarEnum::Match; // most likely anyhow?!
					i = i.saturating_sub(1);
					j =j.saturating_sub(1);
		    	}
		    	else 
		    	*/if let Some((max_index, max_value)) = vec![
		    	self.dp[i.saturating_sub(1)][j.saturating_sub(1)], 
		    	self.dp[i.saturating_sub(1)][j],
		    	self.dp[i][j.saturating_sub(1)]
		    	].iter().enumerate().max_by_key(|(_, &val)| val) {
		    		cigar[rev_id] = match max_index{
		    			0 => {
		    				// a match or a mismatch self.dp[i.saturating_sub(1)][j.saturating_sub(1)]
		    				i = i.saturating_sub(1);
		    				j =j.saturating_sub(1);
		    				if  *max_value == this_value - MATCH_SCORE {
		    					CigarEnum::Match
		    				}else {
		    					CigarEnum::Mismatch
		    				}
		    			},
		    			1=>{
		    				// an instertion is most likely
		    				i = i.saturating_sub( 1 );
		    				CigarEnum::Insertion
		    			},
		    			2=> {
		    				// a deletion is most likely
		    				j = j.saturating_sub(1);
		    				CigarEnum::Deletion
		    			},
		    			_ => unreachable!()
		    		};
		    		if self.debug{
		    			println!("inserted the value {} at position {rev_id}", &cigar[rev_id] );
		    		}
		    		//println!("Here ({i};{j} I had a max index of {max_index} and a max value of {max_value} and decided on a {}",cigar[rev_id]);
		    		
		    	}else {
		    		panic!("I can not decode the df matrix to aCigar state at :{i}; j{j}:\n{}\t{}\n{}\t{}",
		    			self.dp[i.saturating_sub(1)][j.saturating_sub(1)], self.dp[i][j.saturating_sub(1)],
		    			self.dp[i.saturating_sub(1)][j],self.dp[i][j]
		    			);
		    	}
		    //}
		    if rev_id > 0 {
		    	rev_id -= 1;
		    } else if i > 0 || j > 0  {
		    	// we have etimated the wrongpath length!
		    	cigar.insert(0, CigarEnum::Empty);
		    	//println!("Inserting a new CigarEnum::Empty at position 0");
		    	//self.debug = true;
		    }   
	    }

	    if rev_id != 0 {
	    	panic!("We have not filled in all the valiues here!?!?");
	    }

	    // I need that for the debug:
	    let mut cig = Cigar::new("");
	    cig.convert_to_cigar( &cigar );
	    if self.debug {
	    	let (alng1, alng2 ) = self.needleman_wunsch_affine_backtrack( seq1, seq2, &cigar );
	    	println!("Initial cigar string:\n{cig}\n{}\n{}", alng1, alng2 );
	    	println!("The initial alignement:\n{}", self.int_state_to_string( seq1, seq2, &cigar ));
	    }

	    //println!("Before fixup I have:\n{cig}");

	    //println!("Fixing aroung deletion/insertion mapping");

	    // The alignement is consitently bad at mapping bp around a Deletion/Insertion.
	    // It is more likely that bp that would match somewhere in the gap are scattered over the gap,
	    // even if the bp would 100% match the gap start.
	    // This functionality fixed that issue.

	    
	    let mut gap_start: Option<(usize, CigarEnum, usize ) > = None;
		let mut matching: usize;
		let mut seq1_id = seq1.len();
		let mut seq2_id = seq2.len();
		let mut i = cigar.len();
		// Find the position of the last non-deletion cigar operation
		while i > 0 && cigar[i - 1] == CigarEnum::Deletion {
		    i -= 1;
		    seq2_id -=1;
		}
		//let mut restarts = 0;
		while i > 0 {
		    i -= 1;
		    // this is where I try to match the unmatched bits from the opposit side.
		    match gap_start {
		    	// count how many nucleotides might be wrongly placed in the area 
		        None => {
		            if cigar[i] == CigarEnum::Deletion || cigar[i] == CigarEnum::Insertion{
		                matching = 0;
		                let replace_with = cigar[i];
		                let mut drop_replaces = 0;
		                while let (Some(nuc1), Some(nuc2)) = (seq1.get_nucleotide_2bit(seq1_id.saturating_sub(1)), seq2.get_nucleotide_2bit(seq2_id.saturating_sub(1) )) {
		                    if nuc1 != nuc2 {
		                    	if self.debug{
		                    		println!("I detected a {replace_with} at position {i} (seq1_id = {seq1_id}; seq2_id = {seq2_id}) - but {} does not match {}", 
		                    			Self::to_utf8(nuc1) as char, Self::to_utf8(nuc2) as char)
		                    	}
		                    	//seq1_id -= 1;
		                    	//seq2_id -= 1;
		                    	break;
		                    }
		                    if self.debug{
		                    	println!("I detected a {replace_with} at position {i} - and {} does match {}", 
		                    		Self::to_utf8(nuc1) as char, Self::to_utf8(nuc2) as char)
		                    }
		                    matching += 1;
		                    // both sequences had a match so we need to check the next seqence position in both sequences:
		                    seq1_id -= 1;
		                    seq2_id -= 1;
		                   
		                    // this is going to be the match that we just detected
		                    cigar[i] = CigarEnum::Match;
		                    // now let's look in the sequence before the one we just made.
		                    i -= 1;
		                    if cigar[i].opposite(&replace_with) {
		                    	// if we just replace that we change the alignement length!
		                    	// actually this is a wrong alignement here: just drop it and check if that works
		                    	if self.debug{
		                    		println!("#1 position: We drop the {i}th entry - a {} as it is the opposite of our search {}", cigar[i], replace_with);
		                    	}
		                    	cigar.remove(i);
		                    	i-=1;
		                    	drop_replaces+=1;
		                    }
		                    // make sure we overwrite the old value to make the required move of the deletion / insert
		                    else if cigar[i] != replace_with{
		                    	cigar[i] = replace_with;
		                    	if self.debug {
		                    		cig.convert_to_cigar( &cigar );
		                    		println!("#1 at position {i}+1 (now looking into {i}) - I will replace the {} with {} here {cig}",  cigar[i], replace_with);
		                    		println!("The updated alignement:\n{}", self.int_state_to_string( seq1, seq2, &cigar ));
		                    	}

		                        matching -= 1;
		                    }
		                    if i == 0 {
		                        break;
		                    }
		                }
		                if matching > 0 {

		                	gap_start = Some( (matching, replace_with, drop_replaces) );
		                	if self.debug{
		                		println!("#1 I created a gap_start {gap_start:?}");
		                	}
		                }
		                /*
		                println!("We decided that this needs to be changed: {gap_start:?}");
		                cig.convert_to_cigar( &cigar );
		                println!("Before this change the intermediate cigar looks like that:\n{cig}\nremember we will change {matching} 'not {replace_with}' to '{replace_with}'");
		                */
		            }
		        }
		        // move the remaining deletion / insert 
		        Some( (to_ignore, replace_with, drop_replaces) ) => {
		        	// whatever the system thought it would do here - do not do that
		        	if cigar[i] != replace_with {
			        	// at least if it consumed a nucleotiode from both DNA strings - not a good idea!
			            //if cigar[i] == CigarEnum::Match || cigar[i] == CigarEnum::Mismatch {
			            if to_ignore > 0 {
			            	if cigar[i].opposite(&replace_with) {
		                    	// if we just replace that we change the alignement length!
		                    	// actually this is a wrong alignement here: just drop it and check if that works!
		                    	
		                    	cigar.remove(i);
		                    	i-=1;
		                    	gap_start = Some((to_ignore, replace_with, drop_replaces + 1) );
		                    	if self.debug {
		                    		println!("#2 position: We drop the {i}th entry - a {} as it is the opposite of our search {}", cigar[i], replace_with);
		                    		println!("The updated alignement: {}", self.int_state_to_string( seq1, seq2, &cigar ));
		                    	}
		                    }else if cigar[i] != replace_with {
		                    	cigar[i] = replace_with;
		                    	if self.debug{
		                    		cig.convert_to_cigar( &cigar );
		                    		println!("#2 position {i} - replcaing a {} with {} ({cig})", cigar[i], replace_with);
		                    		println!("The updated alignement: {}", self.int_state_to_string( seq1, seq2, &cigar ));
		                    	}
		                    	gap_start = Some((to_ignore - 1, replace_with, drop_replaces) );
		                    }
		                } else {
		                    gap_start = None;
		                }
		            }else if drop_replaces > 0 {
		            	cigar[i] = replace_with;
		                gap_start = Some((to_ignore - 1, replace_with, drop_replaces -1 ) );
		            }
		        }
		    };

		    if cigar[i] != CigarEnum::Deletion {
		    	seq1_id = seq1_id.saturating_sub(1);
		    }
		    if cigar[i] != CigarEnum::Insertion {
		    	seq2_id = seq2_id.saturating_sub(1);
		    }
		}

		if self.debug {
			cig.convert_to_cigar( &cigar );
			let (alng1, alng2 ) = self.needleman_wunsch_affine_backtrack( seq1, seq2, &cigar );
	    	println!("del/ins remapped cigar string:\n{cig}\n{}\n{}", alng1, alng2);
	    }
		//cig.convert_to_cigar( &cigar );
		//println!("{cig}\nIs what we have after and before the fix_1d1i_1i1d() call");

		Cigar::fix_1d1i_1i1d( &mut cigar, None );


		if self.debug {
			cig.convert_to_cigar( &cigar );
			let (alng1, alng2 ) = self.needleman_wunsch_affine_backtrack( seq1, seq2, &cigar );
	    	println!("fix_1d1i_1i1d'ed cigar string:\n{cig}\n{}\n{}", alng1, alng2);
	    }
	    //cig.convert_to_cigar( &cigar );
	    //println!("After fix_1d1i_1i1d() I have:\n{cig}");

	    cigar
	}

	/// a rather useless function converting the Cigar vector to a String - No Cigar string - but the exptended version.
	pub fn cigar_to_string(cigar:&[CigarEnum]) -> String{
		let mut ret = "".to_string();
		for i in cigar{
			ret += &format!("{}",i)
		}
		ret
	}



// Function to export DP matrix to a file
pub fn export_dp_matrix(&self, file_path: &str) -> io::Result<()> {
    let mut file = match File::create(file_path) {
    	Ok(file) => file,
    	Err(err) => {
    		eprintln!( "Trying to write th dp matrix of a needleman wunsch affine mapping I hit this file system error:\n{err:?}");
    		// this is not too crucial here and it would be a shame to kill the whole process over this!
    		return Ok(()) 
    	},
    };

    for row in &self.dp {
        let row_str: Vec<String> = row.iter().map(|&val| val.to_string()).collect();
        let row_line = row_str.join("\t");
        writeln!(file, "{}", row_line)?;
    }
    println!("I have exported the needleman wunsch affine primary matrix to {file_path}");
    Ok(())
}


}


