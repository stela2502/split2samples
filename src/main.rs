use clap::Parser;
use needletail::bitkmer::BitNuclKmer;
use needletail::parse_fastx_file;
use serde::Deserialize;
use std::collections::HashMap;
use std::io::Write;
use std::str;
mod utils;
use crate::utils::get_all_snps;
use std::path::Path;
use std::ffi::OsStr;

use std::collections::BTreeMap;



// first, reproduce the appproach from
// https://github.com/jeremymsimon/SPLITseq/blob/main/Preprocess_SPLITseq_collapse_bcSharing.pl

/// This doc string acts as a help message when the user runs '--help'
/// as do all doc strings on fields
#[derive(Parser)]
#[clap(version = "0.1.0", author = "Rob P. <rob@cs.umd.edu>")]
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
struct sample {
    oligo: &str,
    id: short,
    search: &BTreeMap<&str, &str>,
    file1: &Box<dyn Write>, 
    file2: &Box<dyn Write>,
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


fn get_sample(primer: &str, id: usize, sub_len: usize, outpath: &str) -> &sample {
    let mut sam: sample;
    sam.oligo = primer;
    sam.id = id;
    /// and here comes the fun. Get me all possible sub_len strings in that oligo into the BTreeMap
    for n in 1..sub_len {
	for subS in sub_strings( primer, sub_len){
		sam.search.insert( subS, subS);
	}
    }
    /// start the two files buffers
    let f = Path::new(&opts.outpath);
    let mut File2 = f.join( Path::new(&opts.read_file1).file_name().unwrap());
    sam.file1 = writer( &File2 );
    File2 = f.join( Path::new(&opts.read_file2).file_name().unwrap());
    sam.file2 = writer( &File2 );

    return &sam;
}	

/// from https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561
/// Write normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
// Attempting to have a file writer too
pub fn writer(filename: &Path ) -> &Box<dyn Write> {

    let file = match File::create(&filename) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        // Error is here: Created file isn't gzip-compressed
        return &Box::new(BufWriter::with_capacity(
            128 * 1024,
            write::GzEncoder::new(file, Compression::default()),
        ))
    } else {
        return &Box::new(BufWriter::with_capacity(128 * 1024, file))
    }

}	

pub fn writeFastq( buf Box<dyn Write>, res FastqRecord ) {
    rec = res.unwrap();
    buf.write_all( rec.all() ).unwrap();
}



fn main() {
    // parse the options
    let opts: Opts = Opts::parse();

    let mut samples: Vec<&sample> = Vec::with_capacity(12);

    //let File1 = Path::new(outpath).join( Path::new(read_file1).file_name());
    let mut File1 = Path::new(&opts.outpath);

    let mut File2 = File1.join( Path::new(&opts.read_file1).file_name().unwrap());

    let file1 = writer( &File2 );
    File2 = File1.join( Path::new(&opts.read_file2).file_name().unwrap() );
    let file2 = writer( &File2);

    if ( opts.specie.eq("human") ){
        /// get all the human sample IDs into this.
        samples[0] = get_sample(primer="ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG", 1, 9, outpath);
        samples[1] = get_sample("TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG", 2, 9, outpath);
        samples[2] = get_sample("CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT", 3, 9, outpath);
        samples[3] = get_sample("ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT", 4, 9, outpath);
        samples[4] = get_sample("CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG", 5, 9, outpath);
        samples[5] = get_sample("TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC", 6, 9, outpath);
        samples[6] = get_sample("TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT", 7, 9, outpath);
        samples[7] = get_sample("CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT", 8, 9, outpath);
        samples[8] = get_sample("GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG", 9, 9, outpath);
        samples[9] = get_sample("GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC", 10, 9, outpath);
        samples[10] = get_sample("CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC", 11, 9, outpath);
        samples[11] = get_sample("GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG", 12, 9, outpath);

    }
    else if ( opts.specie.eq("mouse") ){
        // and the mouse ones
        samples[0] = get_sample("AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", 1, 9, outpath);
        samples[1] = get_sample("ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 2, 9, outpath);
        samples[2] = get_sample("AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", 3, 9, outpath);
        samples[3] = get_sample("TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT", 4, 9, outpath);
        samples[4] = get_sample("GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", 5, 9, outpath);
        samples[5] = get_sample("GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC", 6, 9, outpath);
        samples[6] = get_sample("ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", 7, 9, outpath);
        samples[7] = get_sample("TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG", 8, 9, outpath);
        samples[8] = get_sample("GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", 9, 9, outpath);
        samples[9] = get_sample("TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC", 10, 9, outpath);
        samples[10] = get_sample("GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", 11, 9, outpath);
        samples[11] = get_sample(primer="CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC", 12, 9, outpath);
    }
    else {
        println!("Sorry, but I have no primers for specie {}", specie);
        std::process::exit(1)
    }


    // for now, we're assuming FASTQ and not FASTA.
    let mut reader1 = parse_fastx_file(&opts.read_file1).expect("valid path/file");
    let mut reader2 = parse_fastx_file(&opts.read_file2).expect("valid path/file");

    while let Some(record2) = reader2.next() {
        let Some(record1) = reader1.next();	
        let seqrec = record2.expect("invalid record");
        let seq = format("{}",seqrec.seq().into_owned());

	
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
	let maxValue: Integer = res.iter().max();
	let mut z = 0;
        let ID: int;
        if maxValue > 2{
	for i in 0..res.len(){
    	    if res[i] == max {
	        ID = i;
                z += z;
	    }
 	}
	if (z == 1){
	    writeFastq( samples[i].file1, record1 );
	    writeFastq( samples[i].file2, record2 );	
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
