// tests/test_mapping_efficiency.rs

// Import necessary modules
use std::process::Command;


#[test]
fn test_quantify_rhapsody_multi_ok_seq() {
	
	let is_release_mode = !cfg!(debug_assertions);

	// this is mapping to Sry using BLASTN
    // but the sequence is really really bad - better not match?!
	let seq1 = "AGAAGCAGCAGTTCCATGACCACCACCAAAAGAAGAGAAGAATAAAACAACAACAACAAAAACAACAACAAAAC";

    // Execute the command with the provided arguments
    let output = Command::new(
    	if is_release_mode { "./target/release/map_one_sequence" } 
    	else { "./target/debug/map_one_sequence" }
    	).args(&[
    		"-s", "mouse",
            "-o", "testData/output_1e5",
            "-e", "testData/2276_20220531_chang_to_rpl36a_amplicons.fasta",
            "-a", "testData/MyAbSeqPanel.fasta",
            "-d", seq1
        ])
        .output()
        .expect("Failed to execute command");

    let stdout = std::str::from_utf8( &output.stdout).unwrap();

    assert!( stdout.contains( "I have collected these genes: {453: 1}"), "No gene detected!");

}

#[test]
fn test_quantify_rhapsody_multi_failed_strict() {
    
    let is_release_mode = !cfg!(debug_assertions);

    // this is mapping to Sry using BLASTN
    // but the sequence is really really bad - better not match?!
    let seq1 = "AGAAGCAGCAGTTCCATGACCACCACCAAAAGAAGAGAAGAATAAAACAACAACAACAAAAACAACAACAAAAC";

    // Execute the command with the provided arguments
    let output = Command::new(
        if is_release_mode { "./target/release/map_one_sequence" } 
        else { "./target/debug/map_one_sequence" }
        ).args(&[
            "-s", "mouse",
            "-o", "testData/output_1e5",
            "-e", "testData/2276_20220531_chang_to_rpl36a_amplicons.fasta",
            "-a", "testData/MyAbSeqPanel.fasta",
            "-d", seq1,
            "--highest-nw-val", "0.2",
        ])
        .output()
        .expect("Failed to execute command");

    let stdout = std::str::from_utf8( &output.stdout).unwrap();

    assert!( stdout.contains( "I have collected these genes: {}"), "No gene detected!");

}

#[test]
fn test_cxcr4_mapping() {
    
    let is_release_mode = !cfg!(debug_assertions);

    // this is mapping to Sry using BLASTN
    // but the sequence is really really bad - better not match?!
    let seq1 = "GCTGCATAATCTCTTCATTCCGAGGAGCACCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGCGGGGTTGTGGTTGT";

    // Execute the command with the provided arguments
    let output = Command::new(
        if is_release_mode { "./target/release/map_one_sequence" } 
        else { "./target/debug/map_one_sequence" }
        ).args(&[
            "-s", "mouse",
            "-o", "testData/output_1e5",
            "-e", "testData/2276_20220531_chang_to_rpl36a_amplicons.fasta",
            "-a", "testData/MyAbSeqPanel.fasta",
            "-d", seq1,
            "--highest-nw-val", "0.2",
        ])
        .output()
        .expect("Failed to execute command");

    let stdout = std::str::from_utf8( &output.stdout).unwrap();

    assert!( stdout.contains( "I have collected these genes: {49: 1}"), "Cxcr4 49 not found!");

}


#[test]
fn test_quantify_rhapsody_multi_bad_seq() {
	
	let is_release_mode = !cfg!(debug_assertions);

	// this is mapping to Washc1 using BLASTN
	let seq1 = "AAGAAGCAGCAGTTCCATGACCACCACCACAGCAGCCGGCAGGAGACGAGGATGAGGAGGACTGGGAGTCCTA";

    // Execute the command with the provided arguments
    let output = Command::new(
    	if is_release_mode { "./target/release/map_one_sequence" } 
    	else { "./target/debug/map_one_sequence" }
    	).args(&[
    		"-s", "mouse",
            "-o", "testData/output_1e5",
            "-e", "testData/2276_20220531_chang_to_rpl36a_amplicons.fasta",
            "-a", "testData/MyAbSeqPanel.fasta",
            "-d", seq1
        ])
        .output()
        .expect("Failed to execute command");

    let stdout = std::str::from_utf8( &output.stdout).unwrap();

    assert!( stdout.contains( "I have collected these genes: {}"), "Gene detected!");

}