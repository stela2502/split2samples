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

    let command = if is_release_mode {
        "./target/release/quantify_rhapsody_multi"
    } else {
        "./target/debug/quantify_rhapsody_multi"
    };

    let args = &[
        "-r", "testData/1e5_mRNA_S1_R1_001.fastq.gz",
        "-f", "testData/1e5_mRNA_S1_R2_001.fastq.gz",
        "-o", "testData/output_1e5",
        "-s", "mouse",
        "-a", "testData/MyAbSeqPanel.fasta",
        "-e", "testData/genes.fasta",
        "-m", "1",
        "-v", "v1"
    ];  

    // Execute the command with the provided arguments
    let output = Command::new( command ).args( args )
        .output()
        .map_err(|e| {
            eprintln!("Failed to execute command: {}", e);
            e
        }).unwrap();

    let cmd = format!("{} {}", command, args.join(" "));
    if !output.status.success() {
        eprintln!("Command failed: {}", cmd);
        // Handle failure accordingly
    }else {
        println!("{}", cmd );
    }


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

	let mut data = HashMap::<String, isize>::new();

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
    /*
    The data of that part needs to be collected, too!
collected read counts:
expression reads  : 45553 reads (66.45% of cellular)
antibody reads    : 19991 reads (29.16% of cellular)
sample reads      : 993 reads (1.45% of cellular)

reported UMI counts:
expression reads  : 593 UMIs (0.86% of cellular)
antibody reads    : 352 UMIs (0.51% of cellular)
sample reads      : 12 UMIs (0.02% of cellular)
    */

    for line in output_str.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 8 {
            // Extract the name and value
            let name = format!("{} {}", parts[0], parts[4]); // Combine the first CellType and the CountsType
            if name == "-> reads"{
                continue
            }
            if let Ok(value) = parts[3].parse::<isize>() { // Parse the fourth part as the value
                //println!("{}", name); // This prints "expression reads"
                //println!("{}", value); // This prints 593
                data.insert(name, value);
            }
        }
    }


    /*
    SampleTag03_mm should be 245 but was 251
    SampleTag05_mm should be 21 but was 23
    SampleTag02_mm should be 216 but was 220
    SampleTag01_mm should be 140 but was 147
    SampleTag04_mm should be 169 but was 175
    SampleTag06_mm should be 132 but was 139
    na should be 34413 but was 34577
    */

    /*
    SampleTag03_mm should be 250 but was 254
    SampleTag04_mm should be 172 but was 176
    SampleTag02_mm should be 220 but was 224
    SampleTag06_mm should be 138 but was 140
    SampleTag01_mm should be 147 but was 150
    SampleTag05_mm should be 23 but was 24
    na should be 33305 but was 35592
    */

    /*
    SampleTag06_mm should be 140 but was 134
    SampleTag01_mm should be 150 but was 142
    SampleTag04_mm should be 176 but was 171
    na should be 35592 but was 35609
    SampleTag03_mm should be 254 but was 252
    SampleTag05_mm should be 24 but was 22
    SampleTag02_mm should be 224 but was 218
    */

    /*
    SampleTag06_mm should be 141 but was 134
    SampleTag02_mm should be 225 but was 218
    na should be 35558 but was 35533
    SampleTag01_mm should be 151 but was 142
    SampleTag04_mm should be 177 but was 171
    */

    let mut exp =HashMap::<String, isize>::new();
    exp.insert( "na".to_string(), 35533 );
    exp.insert( "SampleTag01_mm".to_string(), 142 );
    exp.insert( "SampleTag02_mm".to_string(), 218 );
    exp.insert( "SampleTag03_mm".to_string(), 254 );
    exp.insert( "SampleTag04_mm".to_string(), 171 );
    exp.insert( "SampleTag05_mm".to_string(), 24 );
    exp.insert( "SampleTag06_mm".to_string(), 134 );
    exp.insert( "AssignedSampleName".to_string(), 1 );

    // overall stats:
    exp.insert( "cellular reads".to_string(), 68556 );
    exp.insert( "filtered reads".to_string(), 14078 );

    //collected read counts:
    exp.insert( "expression reads".to_string(), 45429 );
    exp.insert( "antibody reads".to_string(),   19983 );
    exp.insert( "sample reads".to_string(),     959 );

    //reported UMI counts:
    exp.insert( "expression UMIs".to_string(),  45348 );
    exp.insert( "antibody UMIs".to_string(),    19940 );
    exp.insert( "sample UMIs".to_string(),      959 );

    // Iterate over the actual hashmap and assert each key-value pair separately
    let mut failed = false;
    for (key, actual_value) in &data {
        match exp.get(key){
            Some(expected_value) => {
                if (actual_value - expected_value).abs() > 3 {
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



