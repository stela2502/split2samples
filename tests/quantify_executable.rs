// tests/quantify_executable.rs

// Import necessary modules
use std::process::Command;
use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::fs::File;
use std::fs;
use std::io::BufReader;
use std::io::BufRead;
use std::collections::HashMap;


#[test]
fn test_quantify_rhapsody_multi() {
	
	let is_release_mode = !cfg!(debug_assertions);

    // Execute the command with the provided arguments
    let output = Command::new(
    	if is_release_mode { "./target/release/quantify_rhapsody_multi" } 
    	else { "./target/debug/quantify_rhapsody_multi" }
    	).args(&[
            "-r", "testData/1e5_mRNA_S1_R1_001.fastq.gz",
            "-f", "testData/1e5_mRNA_S1_R2_001.fastq.gz",
            "-o", "testData/output_1e5",
            "-s", "mouse",
            "-e", "testData/genes.fasta",
            "-a", "testData/MyAbSeqPanel.fasta",
            "-m", "1",
            "-v", "v1"
        ])
        .output()
        .expect("Failed to execute command");

	// Check if the command was successful (exit code 0)
    assert!(output.status.success());
    // Convert output to string
    let output_str = String::from_utf8_lossy(&output.stdout);

	// Check if output contains the expected lines
	let file_path ="testData/output_1e5/SampleCounts.tsv";
    assert!(output_str.contains(file_path));

    assert!(fs::metadata(file_path).is_ok(), "expected outfile does exists");

    let file = File::open(file_path).unwrap();
	let reader = BufReader::new(file);

	let mut data = HashMap::<String, usize>::new();

    // Process each line in the file
    for line in reader.lines() {
        //let line = line?;
        // Split the line by tabs
        match line{
        	Ok(line) => {
        		let columns: Vec<&str> = line.split('\t').collect();
        		if columns.len() >= 13 {
        			*data.entry(columns[13].to_string()).or_insert(0) += 1;
        		}
        	},
        	Err(e) => {eprintln!("Oops an error occured? {e:?}")},
        }
    }

    let mut exp =HashMap::<String, usize>::new();
    exp.insert( "na".to_string(), 34292 );
    exp.insert( "Sample1".to_string(), 98 );
    exp.insert( "Sample2".to_string(), 185 );
    exp.insert( "Sample3".to_string(), 236 );
    exp.insert( "Sample4".to_string(), 143 );
    exp.insert( "Sample5".to_string(), 19 );
    exp.insert( "Sample6".to_string(), 115 );
    exp.insert( "AsignedSampleName".to_string(), 1 );

    // Iterate over the actual hashmap and assert each key-value pair separately
    for (key, actual_value) in &data {
        assert_eq!(
            exp.get(key),
            Some(actual_value),
            "Unexpected value for key {}",
            key
        );
    }


    
}
