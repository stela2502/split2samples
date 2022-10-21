use clap::Parser;
//use flate2;
//use flate2::read;
use flate2::write;
use flate2::Compression;
//use needletail::bitkmer::BitNuclKmer;
use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord;

use serde::Deserialize;
//use std::collections::HashMap;
use std::io::Write;
use std::io::BufWriter;
use std::str;
mod utils;
//use crate::utils::get_all_snps;
use std::path::Path;
use std::ffi::OsStr;
use std::fs::File;

use std::collections::BTreeMap;



// first, reproduce the appproach from
// https://github.com/jeremymsimon/SPLITseq/blob/main/Preprocess_SPLITseq_collapse_bcSharing.pl

/// This doc string acts as a help message when the user runs '--help'
/// as do all doc strings on fields
#[derive(Parser)]
#[clap(version = "0.1.0", author = "Rob P. <rob@cs.umd.edu>, Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the input R1 file
    #[clap(short, long)]
    read_file1: String,
    /// the input R2 file
    #[clap(short, long)]
    read_file2: String,
    /// the specie of the library
    #[clap(short, long)]
    specie: String,
    /// the input outpath
    #[clap(short, long)]
    outpath: String,
    /// consider 1-hamming distance neighbors of random hexamers
    #[clap(short, long)]
    one_hamming: bool,
}

#[derive(Debug, Deserialize)]
struct Sample {
    oligo:  &str,
    id: i16,
    search:  BTreeMap< &str, &str>,
    file1:  Box<dyn Write>, 
    file2:  Box<dyn Write>,
}


/// from https://users.rust-lang.org/t/solved-how-to-split-string-into-multiple-sub-strings-with-given-length/10542/7
fn sub_strings(string: &str, sub_len: usize) -> Vec<&str> {
    let mut subs = Vec::with_capacity(string.len() / sub_len);
    let mut iter = string.chars();
    let mut pos = 0;

    while pos < string.len() {
        let mut len = 0;
        for ch in iter.by_ref().take(sub_len) {
            len += ch.len_utf8();
        }
        subs.push(&string[pos..pos + len]);
        pos += len;
    }
    subs
}

/// from https://stackoverflow.com/questions/65976432/how-to-remove-first-and-last-character-of-a-string-in-rust
fn rem_first(value: &str) -> &str {
    let mut chars = value.chars();
    chars.next();
    chars.as_str()
}


fn get_sample(sample: &Sample, primer: &str, id: usize, sub_len: usize, outpath: &str, read1: &str, read2: &str) {
    sample.oligo = primer;
    sample.id = id;
    // and here comes the fun. Get me all possible sub_len strings in that oligo into the BTreeMap
    for n in 1..sub_len {
	  for subS in sub_strings( primer, sub_len){
		sample.search.insert( subS, subS);
	  }
    }
    // start the two files buffers
    let f = Path::new( outpath );
    let mut File2 = f.join( Path::new( read1 ).file_name().unwrap());
    sample.file1 = writer( &File2 );
    File2 = f.join( Path::new( read2 ).file_name().unwrap());
    sample.file2 = writer( &File2 );

}	

/// from https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561
/// Write normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
// Attempting to have a file writer too
pub fn writer(filename: &Path ) -> &Box<dyn Write> {

    let file = match File::create(&filename) {
        Err(why) => panic!("couldn't open {}: {}", filename.display(), why.description()),
        Ok(file) => file,
    };

    if filename.extension() == Some(OsStr::new("gz")) {
        // Error is here: Created file isn't gzip-compressed
        return &Box::new(BufWriter::with_capacity(
            128 * 1024,
            write::GzEncoder::new(file, Compression::default()),
        ))
    } else {
        return &Box::new(BufWriter::with_capacity(128 * 1024, file))
    }

}	

pub fn writeFastq( buf: Box<dyn Write>, res: SequenceRecord ) {
    res = res.unwrap();
    buf.write_all( res.all() ).unwrap();
}



fn main() {
    // parse the options
    let opts: Opts = Opts::parse();

    let mut samples: Vec<&Sample> = Vec::with_capacity(12);
    //let File1 = Path::new(outpath).join( Path::new(read_file1).file_name());
    let mut File1 = Path::new(&opts.outpath);

    let mut File2 = File1.join( Path::new(&opts.read_file1).file_name().unwrap());

    let file1 = writer( &File2 );
    File2 = File1.join( Path::new(&opts.read_file2).file_name().unwrap() );
    let file2 = writer( &File2);

    if  &opts.specie.eq("human") {
        // get all the human sample IDs into this.
        let mut sam1: Sample;
        samples[0] = sam1;
        get_sample(&sam1, "ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG", 1, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam2: Sample;
        samples[1] = sam2;
        get_sample(&sam2, "TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG", 2, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam3: Sample;
        samples[2] = sam3;
        get_sample(&sam3, "CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT", 3, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam4: Sample;
        samples[3] = sam4;
        get_sample(&sam4, "ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT", 4, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam5: Sample;
        samples[4] = sam5;
        get_sample(&sam5, "CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG", 5, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam6: Sample;
        samples[5] = sam6;
        get_sample(&sam6, "TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC", 6, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam7: Sample;
        samples[6] = sam7;
        get_sample(&sam7, "TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT", 7, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );    
        let mut sam8: Sample;
        samples[7] = sam8;
        get_sample(&sam8, "CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT", 8, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam9: Sample;
        samples[8] = sam9;
        get_sample(&sam9, "GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG", 9, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam10: Sample;
        samples[9] = sam10;
        get_sample(&sam10, "GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC", 10, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam11: Sample;
        samples[10] = sam11;
        get_sample(&sam11, "CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC", 11, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam12: Sample;
        samples[11] = sam12;
        get_sample(&sam12, "GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG", 12, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );

    }
    else if &opts.specie.eq("mouse") {
        // and the mouse ones
        let mut sam1: Sample;
        samples[0] = sam1;
        get_sample(&sam1, "AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", 1, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam2: Sample;
        samples[1] = sam2;
        get_sample(&sam2, "ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 2, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam3: Sample;
        samples[2] = sam3;
        get_sample(&sam3, "AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", 3, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam4: Sample;
        samples[3] = sam4;
        get_sample(&sam4, "TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT", 4, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam5: Sample;
        samples[4] = sam5;
        get_sample(&sam5, "GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", 5, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam6: Sample;
        samples[5] = sam6;
        get_sample(&sam6, "GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC", 6, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam7: Sample;
        samples[6] = sam7;
        get_sample(&sam7, "ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", 7, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );    
        let mut sam8: Sample;
        samples[7] = sam8;
        get_sample(&sam8, "TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG", 8, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam9: Sample;
        samples[8] = sam9;
        get_sample(&sam9, "GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", 9, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam10: Sample;
        samples[9] = sam10;
        get_sample(&sam10, "TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC", 10, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam11: Sample;
        samples[10] = sam11;
        get_sample(&sam11, "GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", 11, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
        let mut sam12: Sample;
        samples[11] = sam12;
        get_sample(&sam12, "CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC", 12, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 );
    }
    else {
        println!("Sorry, but I have no primers for specie {}", &opts.specie);
        std::process::exit(1)
    }


    // for now, we're assuming FASTQ and not FASTA.
    let mut reader1 = parse_fastx_file(&opts.read_file1).expect("valid path/file");
    let mut reader2 = parse_fastx_file(&opts.read_file2).expect("valid path/file");

    while let Some(record2) = reader2.next() {
        let Some(record1) = reader1.next();	
        let seqrec = record2.expect("invalid record");
        let seq = seqrec.seq().into_owned();

	
	// create 9mers of the R2 read and check which of teh 12 sample ids matches best:
	let mut subs = sub_strings( &seq, 9 );
    let mut res = Vec::with_capacity(12);
 
	for i in 0..11{
	    res[i] = 0;
	    for s in subs{
		  if samples[i].search.contains_key( s ){
		    res[i] +=1;
		  }
	    }
	}
	let maxValue = res.iter().max();
	let mut z = 0;
    let mut id = -1;
    if maxValue > 2{
        for i in 0..res.len(){
    	    if res[i] == maxValue {
	           id = i;
            z += z;
            }
	    }
 	}
	if z == 1{
	    writeFastq( samples[id].file1, record1 );
	    writeFastq( samples[id].file2, record2 );	
	}
	else {
	    writeFastq( file1, record1 );
	    writeFastq( file2, record1 );
	}

    }

    for i in 0..11{
	   samples[i].file1.flush().unwrap();
	   samples[i].file2.flush().unwrap();
    }
    file1.flush().unwrap();
    file2.flush().unwrap();
}
