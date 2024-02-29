// tests/quantify_executable.rs

// Import necessary modules
use std::process::Command;
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
    exp.insert( "na".to_string(), 34413 );
    exp.insert( "SampleTag01_mm".to_string(), 140 );
    exp.insert( "SampleTag02_mm".to_string(), 216 );
    exp.insert( "SampleTag03_mm".to_string(), 245 );
    exp.insert( "SampleTag04_mm".to_string(), 169 );
    exp.insert( "SampleTag05_mm".to_string(), 21 );
    exp.insert( "SampleTag06_mm".to_string(), 132 );
    exp.insert( "AssignedSampleName".to_string(), 1 );

    // Iterate over the actual hashmap and assert each key-value pair separately
    let mut failed = false;
    for (key, actual_value) in &data {
        match exp.get(key){
            Some(expected_value) => {
                if actual_value != expected_value {
                    eprintln!("{key} should be {expected_value} but was {actual_value}");
                    failed = true;
                }
            }
            None => {
                eprintln!("{key} should not be found at all but was {actual_value}");
                failed = true;
            }
        }
    }

    if failed{
        panic!("Sample detcetion is faulty!");
    }
    
}
