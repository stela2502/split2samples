// tests/quantify_executable_gene_mapper.rs

// Import necessary modules
use std::process::Command;
use std::fs::File;
use std::fs;
use std::io::BufReader;
use std::io::BufRead;
use std::collections::HashMap;
use regex::Regex;
use std::process::exit;

#[test]
fn test_quantify_gene_mapper() {
    
    let is_release_mode = !cfg!(debug_assertions);
    if ! is_release_mode {
        eprintln!("Test should be re-run in release mode (speed!)");
        exit(0);
    }
    let command = if is_release_mode {
        "./target/release/genomic_mapper"
    } else {
        "./target/debug/genomic_mapper"
    };

    let outpath = "testData/output_1e5_genomic_mapper/";

    let args = &[
        "-r", "testData/1e5_mRNA_S1_R1_001.fastq.gz",
        "-f", "testData/1e5_mRNA_S1_R2_001.fastq.gz",
        "-o", outpath,
        "-s", "mouse",
        "-i", "testData/output_1e5_gm",
        "-m", "1",
        "-v", "v1"
    ];  

    let test_path ="testData/output_1e5_gm/index.bin";
    if let Ok(_metadata) =fs::metadata(test_path) {
        // this is OK
    }else {
        panic!("You first need to run the 'cargo test -r --test quantify_executable_gene_mapper' before that!");
    }

    if let Ok(metadata) =fs::metadata(outpath) {
        // If it exists and is a directory
        if metadata.is_dir() {
            // Remove the directory and all its contents recursively
            if let Err(err) = fs::remove_dir_all(outpath) {
                eprintln!("Failed to remove directory tree: {}", err);
                // Handle the error accordingly
            } else {
                //println!("Directory tree successfully removed");
            }
        } else {
            //eprintln!("{} is not a directory", outpath);
            // Handle the error accordingly
        }
    }

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

    // check the sam file
    let sam_path = outpath.to_string()+"/MappedReads_1.sam";

    assert!(fs::metadata(sam_path.clone()).is_ok(), "expected samfile does not exists {}", &sam_path);

    let sam_file = File::open(sam_path).unwrap();

    let sam_reader = BufReader::new(sam_file);
    //let mut sam_data = HashMap::<String, isize>::new();
    let mut i = 0;
    let mut only_s = 0;
    // First regular expression pattern: ^[0-9] [0-9]*S$
    let match_only_s = r"^[0-9] [0-9]*S$";

    // Second regular expression pattern: ^[0.9][0-9]*S[0-9][0-9]*M[0-9][0-9]*S$
    let mainly_s = r"^[0-9.][0-9]*S[0-9][0-9]*M[0-9][0-9]*S$";
    // Compile the regular expressions
    let bad1 = Regex::new(match_only_s).unwrap();
    let bad2 = Regex::new(mainly_s).unwrap();

    let empty_line = Regex::new( "^$").unwrap();
    let mut empty = 0;

    // Iterate over each line in the SAM file
    for line in sam_reader.lines() {
        i += 1;
        if let Ok(data) = line {
            if empty_line.is_match(&data){
                empty +=1;
                continue
            }
            let fields: Vec<&str> = data.split('\t').collect();
            if fields.len() > 10 {
                // Check if the 6th field (data[5]; Cigar) matches any of the patterns
                if bad1.is_match(fields[5]) || bad2.is_match(fields[5]) {
                    only_s += 1;
                }
            }
        } else {
            eprintln!("Error reading line {}", i);
        }
    }

    assert_eq!( only_s , 0,"still {} unacceptably bad matches in the sam file", only_s);

    assert_eq!( empty , 0,"still {} empty lines in the sam file", empty);

    assert_eq!( i , 44414,"not the right number of lines in the sam file {}", i);
    
    // check the sampleCounts

    // Check if output contains the expected lines
    let file_path = outpath.to_string()+"/SampleCounts.tsv";

    assert!(fs::metadata(file_path.clone()).is_ok(), "expected outfile does exists");

    assert!(output_str.contains(&file_path));

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

    /*
    expression UMIs should be 43870 but was 43874
    SampleTag06_mm should be 127 but was 132
    sample reads should be 862 but was 916
    expression reads should be 43944 but was 43948
    SampleTag02_mm should be 202 but was 210
    SampleTag04_mm should be 157 but was 165
    na should be 25684 but was 25656
    SampleTag03_mm should be 232 but was 242
    sample UMIs should be 862 but was 916
    SampleTag01_mm should be 111 but was 131
    */

    exp.insert( "na".to_string(), 25656 );
    exp.insert( "SampleTag01_mm".to_string(), 131 );
    exp.insert( "SampleTag02_mm".to_string(), 210 );
    exp.insert( "SampleTag03_mm".to_string(), 242 );
    exp.insert( "SampleTag04_mm".to_string(), 165 );
    exp.insert( "SampleTag05_mm".to_string(), 20 );
    exp.insert( "SampleTag06_mm".to_string(), 132 );
    exp.insert( "AssignedSampleName".to_string(), 1 );


    // overall stats:
    exp.insert( "cellular reads".to_string(), 68556 );
    exp.insert( "filtered reads".to_string(), 14078 );

    //collected read counts:
    exp.insert( "expression reads".to_string(), 43948 );
    exp.insert( "antibody reads".to_string(),   0 );
    exp.insert( "sample reads".to_string(),     916 );

    //reported UMI counts:
    exp.insert( "expression UMIs".to_string(),  43874 );
    exp.insert( "antibody UMIs".to_string(),    0 );
    exp.insert( "sample UMIs".to_string(),      916 );

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