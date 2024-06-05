// /tests/create_genmome_index.rs

// tests/quantify_executable_gene_mapper.rs

#[cfg(test)]
mod tests {
    use std::process::Command;
    use std::fs;
    use rustody::genes_mapper::{GenesMapper, NeedlemanWunschAffine};

    fn check_sequence(index: &GenesMapper, seq:&[u8], gname:Option<&str>, unique_name: Option<&str>, nwa: &mut NeedlemanWunschAffine ) {

    	match index.get_strict( seq, 1, nwa) {
            Ok(matched) => {
                //panic!("I got this result: {}", matched[0]);
                match gname {
                	Some(gname) => {
                		assert!(true, "Matched sequence found");
                		assert_eq!( matched[0].get_name(), gname, "Found the right gene back");
                        match unique_name {
                            Some(name) => {
                                match index.get_gene( matched[0].gene_id() ){
                                    Some(gene_data) => {
                                        assert_eq!(gene_data.get_unique_name(), name, "The transcript name did not match the expectations {}", String::from_utf8_lossy(seq) );
                                    },
                                    None => {
                                        assert!(false, "The expected gene_data object could not be found")
                                    }
                                }
                            },
                            None => {
                                // that is fine as we do not need to check anything if we do not have the request!
                            }
                        }
                	},
                	None => {
                		// we expected a faliure
                		assert!( matched.is_empty(), "No result found - as expected {matched:?}");
                	}
                }
            }
            Err(e) => {
            	match gname {
                	Some(gname) => {
                		println!("I got an unexpected error: {:?} as I had expected to find {gname}", e);
                		assert!(false, "An unexpected error occurred");
                	},
                	None => {
                		assert!( true, "Yes there is no match!" )
                	}
                }
            }
        };
    }
    #[test]
    fn test_genome_index_creation() {
        let is_release_mode = !cfg!(debug_assertions);

        let command = if is_release_mode {
            "./target/release/create_gene_mapper_index"
        } else {
            "./target/debug/create_gene_mapper_index"
        };

        let outpath = "testData/test_index/KI270728.1";

        let args = &[
            "-g", "testData/KI270728.1.gtf.gz",
            "-o", outpath,
            "-f", "testData/KI270728.1.fa.gz",
            "--genename", "gene_id",
            "--transcript", "transcript_id",
        ];

        if let Ok(metadata) = fs::metadata(outpath) {
            if metadata.is_dir() {
                if let Err(err) = fs::remove_dir_all(outpath) {
                    eprintln!("Failed to remove directory tree: {}", err);
                }
            }
        }

        let output = Command::new(command).args(args)
            .output()
            .expect("Failed to execute command");

        if !output.status.success() {
            panic!("Command failed: {}", format!("{} {}", command, args.join(" ")));
        }

        let output_str = String::from_utf8_lossy(&output.stdout);
        println!("The output from the program:\n{output_str}");

        let test_path = format!("{}/index.bin", outpath);

        assert!(fs::metadata(&test_path).is_ok(), "The index outfile has not been created!");

        let fasta = format!("{}/indexed_sequences.fa.gz", outpath);

        assert!(fs::metadata(&fasta).is_ok(), "The fasta.gz outfile has not been created!");

        let index = match GenesMapper::load_index( outpath ){
        	Ok(ind) => ind,
        	Err(e) => panic!("The index could not be loaded! {e}"),
        };

        assert_eq!(index.len(), 13, "In total 14 genes/transcripts were indexed");

        let mut nwa = NeedlemanWunschAffine::new();

        // two different transcription ends - do we get them both?
        // matches to both of them - I can not check for the transcript name as it is a randoim process which one gets reported - they are both equal at this test
        check_sequence(&index, b"GCCAGGTGGTCCTTACCATGACCAACATGGACCCTGTGGACACAGCCACGTATTACTGT", Some("ENSG00000278510"), None, &mut nwa );
        // matches to only ENST00020619806 (artificially created!!)
		check_sequence(&index, b"ATGCCTCCTGTACAAGAACCCAGGCTGCGTCTCAGTGGTGCTCCCTCCCTACCTCTGCAGAACAGGAAAGT", Some("ENSG00000278510"),  Some("ENST00020619806"), &mut nwa );


        // Now the transcript in the opposite direction
        // The transcript
        check_sequence(&index, b"GGACATTCTTACTGTGCTAAAAAGCCACTGCAAACATAGCAATAAAAACCTGTCATTTTCCAAAG", Some("ENSG00000273554"), None, &mut nwa );


        // the revers of that but not in the right location?!
        check_sequence(&index, b"ACGGCCGTGTATTACTGTACCACAGGGGGAGGGGTC", Some("ENSG00000277761"), None,  &mut nwa );


        /*match index.get_strict(b"ACTGTGCTAAAAAGCCACTGCAAACATAGCAATAAAAACCTGTCATTTTCCAAAGCAGG", 1, &mut nwa) {
            Ok(matched) => {
                //panic!("I got this result: {}", matched[0]);
                assert!(true, "Matched sequence found");
                assert_eq!( matched[0].get_name(), "ENSG00000277761", "Found the right gene back")
            }
            Err(e) => {
                println!("I got an unexpected error: {:?}", e);
                assert!(false, "An unexpected error occurred");
            }
        };*/

        //GTCTTCACGCATCCCTTGAATTGGAAATTGTGCCCTGGAGACTGTATACAAGACTGGATTAAAAAGACTGGATTCTGTTT



    }
}
