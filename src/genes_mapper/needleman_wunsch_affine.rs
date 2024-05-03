const MATCH_SCORE: i32 = 1;
const MISMATCH_SCORE: i32 = -1;
const GAP_OPEN_PENALTY: i32 = -6; // Penalty for opening a gap
const GAP_EXTENSION_PENALTY: i32 = -2; // Penalty for extending a gap

use crate::genes_mapper::cigar::{CigarEnum};
use crate::traits::BinaryMatcher;


pub struct NeedlemanWunschAffine{
	dp: Vec<Vec<i32>>,
	n: usize,
	m: usize,
	cigar_vec: Option<Vec<CigarEnum>>,
	circles:usize,
}

impl NeedlemanWunschAffine {
	pub fn new( size:usize ) ->Self {
		// Initialize the DP matrix
		let mut dp = vec![vec![0; size + 1]; size + 1];
		// Initialize the first row and column with gap penalties
    	for i in 1..=size {
	        dp[i][0] = dp[i - 1][0] + GAP_OPEN_PENALTY + (i - 1) as i32 * GAP_EXTENSION_PENALTY;
	    }
	    for j in 1..=size {
	        dp[0][j] = dp[0][j - 1] + GAP_OPEN_PENALTY + (j - 1) as i32 * GAP_EXTENSION_PENALTY;
	    }
	    Self {
	    	dp,
	    	n: size,
	    	m:size,
	    	cigar_vec: None,
	    	circles:0,
	    }
	}

	pub fn to_string<T>( &mut self, seq1: &T, seq2: &T, humming_cut: f32 ) -> String 
	where
    T: BinaryMatcher + std::fmt::Display{

		let ( aligned1, aligned2, cigar_vec ) = self.needleman_wunsch_affine_backtrack(seq1, seq2, humming_cut );
		let cig_str = Self::cigar_to_string( &cigar_vec );
		format!("read:\n{}\ndatabse\n{}\nalignement:\n{:?}\n{:?}\n{cig_str:?}", 
			&seq1, &seq2, 
			String::from_utf8_lossy( &aligned1 ), 
			String::from_utf8_lossy( &aligned2 )
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
	    if self.dp.len() < n +10 || self.dp[0].len() < m +10 {
	    	eprintln!("dp dimensions need to be adjusted to ({n} +10 ; {m} +10 )");
	    	// Initialize the DP matrix
	    	self.dp = vec![vec![0; m + 10]; n + 10];
	    	// Initialize the first row and column with gap penalties
	    	for i in 1..=n {
		        self.dp[i][0] = self.dp[i - 1][0] + GAP_OPEN_PENALTY + (i - 1) as i32 * GAP_EXTENSION_PENALTY;
		    }
		    for j in 1..=m {
		        self.dp[0][j] = self.dp[0][j - 1] + GAP_OPEN_PENALTY + (j - 1) as i32 * GAP_EXTENSION_PENALTY;
		    }
	    }

	    // Fill in the DP matrix
	    for i in 1..=n {
	        for j in 1..=m {
	        	if let (Some(nuc1), Some(nuc2)) = (seq1.get_nucleotide_2bit(i - 1), seq2.get_nucleotide_2bit(j - 1) ){
		            let match_mismatch_score = if nuc1 == nuc2 { MATCH_SCORE } else { MISMATCH_SCORE };
		            let mut scores = vec![
		                self.dp[i - 1][j - 1] + match_mismatch_score, // Diagonal (match/mismatch)
		                self.dp[i - 1][j] + GAP_OPEN_PENALTY + GAP_EXTENSION_PENALTY, // Gap in seq1
		                self.dp[i][j - 1] + GAP_OPEN_PENALTY + GAP_EXTENSION_PENALTY, // Gap in seq2
		            ];
		            // Check if an extension of an existing gap is more favorable
		            for k in 1..=i - 1 {
		                scores[1] = scores[1].max(self.dp[i - k][j] + GAP_OPEN_PENALTY + k as i32 * GAP_EXTENSION_PENALTY);
		            }
		            for k in 1..=j - 1 {
		                scores[2] = scores[2].max(self.dp[i][j - k] + GAP_OPEN_PENALTY + k as i32 * GAP_EXTENSION_PENALTY);
		            }
		            self.dp[i][j] = scores.iter().copied().max().unwrap();
		        }else {
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
	pub fn needleman_wunsch_affine_backtrack<T>(&mut self, seq1: &T, seq2: &T, humming_cut: f32 ) -> (Vec<u8>, Vec<u8>, Vec<CigarEnum>) 
	where
    T: BinaryMatcher {
		let (mut i, mut j) = (seq1.len(), seq2.len());
		let max= i.max(j);
	    let mut align_seq1 = Vec::<u8>::with_capacity( max );
	    let mut align_seq2 = Vec::<u8>::with_capacity( max );
	    i = 0;
	    j = 0;
	    let cig_vec = self.to_cigar_vec( seq1, seq2, humming_cut );

	    for val in &cig_vec {
	    	match val {
	    		CigarEnum::Insertion => {
	    			if let Some(nuc1) = seq1.get_nucleotide_2bit(i){
	    				align_seq1.push( Self::to_utf8(nuc1) );
	            		align_seq2.push( b'-' );
	            		i += 1;
	    			}else {
	    				panic!("Could not get sequence for nuc1")
	    			}
	    		},
	    		CigarEnum::Deletion => {
	    			if let Some(nuc2) = seq2.get_nucleotide_2bit(j){
	    				align_seq1.push( b'-' ); 
						align_seq2.push( Self::to_utf8(nuc2) );
						j +=1;
	    			}else {
	    				panic!("Could not get sequence for nuc2")
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

	    (align_seq1, align_seq2, cig_vec)
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
	pub fn to_cigar_vec<T>(&mut self, seq1: &T, seq2: &T, humming_cut: f32) -> Vec<CigarEnum> 
	where
    T: BinaryMatcher {

		let (mut i, mut j) = (seq1.len(), seq2.len());
	    let mut cigar = vec![CigarEnum::Deletion;i.max(j)];
	    
	    let mut rev_id = cigar.len().saturating_sub(1);

	    while i > 0 || j > 0 {
	    	if let (Some(nuc1), Some(nuc2)) = (seq1.get_nucleotide_2bit(i.saturating_sub(1)), seq2.get_nucleotide_2bit(j.saturating_sub(1)) ){
		        if i > 0 && j > 0 && self.dp[i][j] == self.dp[i.saturating_sub(1)][j.saturating_sub(1)] + (if nuc1 == nuc2 { MATCH_SCORE } else { MISMATCH_SCORE }) {
		        	cigar[rev_id] = if nuc1 == nuc2 { CigarEnum::Match } else { CigarEnum::Mismatch } ;
		            i = i.saturating_sub(1);
		            j =j.saturating_sub(1);
		        }else if i > 0 && (self.dp[i][j] == self.dp[i.saturating_sub(1)][j] + GAP_OPEN_PENALTY || self.dp[i][j] == self.dp[i.saturating_sub(1)][j] + GAP_EXTENSION_PENALTY) {
		        	// CRAP - this is not good!
		        	// this algorithm SUCKS with insertions
		        	// better invert the whole issue here
		        	return self.to_cigar_vec_inverted( seq1, seq2, humming_cut);
		        	/*
		        	// not sure how to handle this here!
		        	cigar[rev_id] = CigarEnum::Insertion;
		            i = i.saturating_sub( 1 );
		            */
		        } else {
		        	cigar[rev_id] = CigarEnum::Deletion ;
		            j = j.saturating_sub(1);
		        }
		        if rev_id == 0 {
		        	break
		        }
		        rev_id = rev_id.saturating_sub( 1 );
		       } else {
		       		panic!("Could not get sequence for nuc1 or nuc2 to_cigar_vec ( {i}, {j} )")
		       }
	    }

	    // The alignement is consitently bad at mapping bp around a deletion.
	    // It is more likely that bp that would match somewhere in the gap are scattered over the gap,
	    // even if the bp would 100% match the gap start.
	    // This functionality fixed that issue.

	    
	    let mut gap_start: Option<usize> = None;
	    let mut matching:usize;
	    let mut seq1_id = 0;
	    let mut seq2_id = 0;
	    let mut i = 0;
	    while i < cigar.len(){
	    	//eprintln!("(round #{i})" );
	    	match gap_start{
	    		None => {
	    			if cigar[i] == CigarEnum::Deletion{
	    				//eprintln!("(round #{i}) Starting the gap at {}; {} (seq1#{seq1_id}, seq2#{seq2_id})", String::from_utf8_lossy(&vec![seq1[seq1_id]]), String::from_utf8_lossy(&vec![seq2[seq2_id]]));
	    				matching = 0;
	    				while let (Some(nuc1), Some(nuc2)) = (seq1.get_nucleotide_2bit(seq1_id), seq2.get_nucleotide_2bit(seq2_id) ) {
	    					if nuc1 != nuc2 {
	    						break
	    					}
							matching +=1;
							seq1_id +=1;
							seq2_id +=1;
							cigar[i] = CigarEnum::Match;
							i +=1;
							if cigar[i] != CigarEnum::Deletion{
								cigar[i] = CigarEnum::Deletion;
								matching -=1;
							}
	    				}
	    				//eprintln!("I have found a total of {matching} matching bp at the start of the gap!");
	    				gap_start = Some( matching );
	    			}
	    		},
	    		Some(to_ignore) => {
	    			if cigar[i] != CigarEnum::Deletion{
	    				if to_ignore > 0 {
	    					gap_start = Some(to_ignore -1);
	    					cigar[i] = CigarEnum::Deletion;
	    				}else {
	    					// the logics has broken. Could still be a mismatch, but this is not straight forward to check here.
	    					//eprintln!("Locis has boken seq2[{seq1_id}] ({}) != seq1[{seq2_id}] ({})", String::from_utf8_lossy(&vec![seq2[seq2_id]]), String::from_utf8_lossy(&vec![seq1[seq1_id]]));
	    					gap_start = None;
	    				}
	    			}
	    		}

	    	};

			if cigar[i] != CigarEnum::Deletion{
				seq1_id +=1;
				if seq1_id == seq1.len(){
					break;
				}
			}
			if cigar[i] != CigarEnum::Insertion{
				seq2_id +=1;
				if seq2_id == seq2.len(){
					break;
				}
			}
			i+=1;
	    }

	    cigar
	}

	/// a rather useless function converting the Cigar vector to a String - No Cigar string - but the exptended version.
	pub fn cigar_to_string(cigar:&[CigarEnum]) -> String{
		let mut ret = "".to_string();

		for value in cigar{
			ret += match value{
				CigarEnum::Match => "M",
				CigarEnum::Mismatch => "X",
				CigarEnum::Insertion => "I",
				CigarEnum::Deletion => "D",
			}
		}
		ret
	}
}


