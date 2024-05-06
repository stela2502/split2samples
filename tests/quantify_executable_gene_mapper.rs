// tests/quantify_executable_gene_mapper.rs

// Import necessary modules
use std::process::Command;
use std::fs::File;
use std::fs;
use std::io::BufReader;
use std::io::BufRead;
use std::collections::HashMap;

#[test]
fn test_quantify_gene_mapper() {
    
    let is_release_mode = !cfg!(debug_assertions);

    let command = if is_release_mode {
        "./target/release/quantify_gene_mapper"
    } else {
        "./target/debug/quantify_gene_mapper"
    };

    let args = &[
        "-r", "testData/1e5_mRNA_S1_R1_001.fastq.gz",
        "-f", "testData/1e5_mRNA_S1_R2_001.fastq.gz",
        "-o", "testData/output_1e5_gm",
        "-s", "mouse",
        "-e", "testData/genes.fasta",
        "-a", "testData/MyAbSeqPanel.fasta",
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
        // print the command here as all errors now also show this print out!
        println!("{}", cmd );
    }

    // Check if the command was successful (exit code 0)
    assert!(output.status.success());
    // Convert output to string
    let output_str = String::from_utf8_lossy(&output.stdout);

    // Check if output contains the expected lines
    let file_path ="testData/output_1e5_gm/SampleCounts.tsv";
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
                //println!("adding {name} -> {value}");
                data.insert(name, value);
            }
        }
    }

    /*
    ---- test_quantify_gene_mapper stdout ----
    SampleTag01_mm should be 142 but was 151
    SampleTag03_mm should be 252 but was 254
    SampleTag05_mm should be 22 but was 24
    SampleTag06_mm should be 134 but was 141
    SampleTag02_mm should be 218 but was 225
    na should be 35533 but was 35559
    SampleTag04_mm should be 171 but was 177
    thread 'test_quantify_gene_mapper' panicked at tests/quantify_executable.rs:240:9:
    Sample detcetion is faulty!

    */

    /*
    na should be 35558 but was 35559
    expression UMIs should be 45463 but was 45464
    expression reads should be 45554 but was 45555
    */
    let mut exp =HashMap::<String, isize>::new();
    /*
    exp.insert( "na".to_string(), 35533 );
    exp.insert( "SampleTag01_mm".to_string(), 142 );
    exp.insert( "SampleTag02_mm".to_string(), 218 );
    exp.insert( "SampleTag03_mm".to_string(), 254 );
    exp.insert( "SampleTag04_mm".to_string(), 171 );
    exp.insert( "SampleTag05_mm".to_string(), 24 );
    exp.insert( "SampleTag06_mm".to_string(), 134 );
    exp.insert( "AssignedSampleName".to_string(), 1 );
    */
    exp.insert( "na".to_string(), 35149 );
    exp.insert( "SampleTag01_mm".to_string(), 129 );
    exp.insert( "SampleTag02_mm".to_string(), 209 );
    exp.insert( "SampleTag03_mm".to_string(), 242 );
    exp.insert( "SampleTag04_mm".to_string(), 163 );
    exp.insert( "SampleTag05_mm".to_string(), 20 );
    exp.insert( "SampleTag06_mm".to_string(), 132 );
    exp.insert( "AssignedSampleName".to_string(), 1 );


    // overall stats:
    exp.insert( "cellular reads".to_string(), 68556 );
    exp.insert( "filtered reads".to_string(), 14078 );

    //collected read counts:
    exp.insert( "expression reads".to_string(), 44506 );
    exp.insert( "antibody reads".to_string(),   19648 );
    exp.insert( "sample reads".to_string(),     911 );

    //reported UMI counts:
    exp.insert( "expression UMIs".to_string(),  44431 );
    exp.insert( "antibody UMIs".to_string(),    19606 );
    exp.insert( "sample UMIs".to_string(),      911 );

    // Iterate over the actual hashmap and assert each key-value pair separately
    let mut failed = false;
    for (key, actual_value) in &data {
        match exp.get(key){
            Some(expected_value) => {
                if (actual_value - expected_value).abs() > 3_isize {
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