use clap::Parser;
//use flate2;
//use flate2::read;
use flate2::Compression;
use flate2::write::GzEncoder;

//use needletail::bitkmer::BitNuclKmer;
use needletail::{parse_fastx_file, Sequence, FastxReader};

//use serde::Deserialize;
//use std::collections::HashMap;
use std::io::Write;
use std::io::BufWriter;
use std::str;
use anyhow::{Context, Result};

mod utils;
//use crate::utils::get_all_snps;
use std::path::Path;
use std::path::PathBuf;
use std::ffi::OsStr;
use std::fs::File;
use std::fs;

use kmers::naive_impl::Kmer;

use std::collections::BTreeSet;



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

struct Sample {
    oligo:  String,
    id: usize,
    search: BTreeSet<u64>,
    file1:  GzEncoder<File>, 
    file2:  GzEncoder<File>,
}

impl Sample {
    fn from_description(primer: &[u8], id: usize, sub_len: usize, outpath: &str, read1: &str, read2: &str) -> Self {

        let mut search = BTreeSet::<u64>::new();
        // and here comes the fun. Get me all possible sub_len strings in that oligo into the BTreeMap
        for kmer in needletail::kmer::Kmers::new(primer, 9) {
            let km = Kmer::from(kmer).into_u64();
            search.insert(km);
        }

        let file1_path = PathBuf::from(outpath).join(read1).join(format!("{}.fq.gz", id));
        let file2_path = PathBuf::from(outpath).join(read2).join(format!("{}.fq.gz", id));
        
        // need better error handling here too
        let file1 = GzEncoder::new(File::create(file1_path).unwrap(), Compression::default());
        let file2 = GzEncoder::new(File::create(file2_path).unwrap(), Compression::default());

        Self {
            oligo: String::from_utf8_lossy(primer).into_owned(), // handle errors better
            id,
            search,
            file1,
            file2
        }
    }
}

fn fill_kmer_vec<'a>(seq: needletail::kmer::Kmers<'a>, kmer_vec: &mut Vec<u64>) {
    kmer_vec.clear();
   for km in seq {
        kmer_vec.push(Kmer::from(km).into_u64());
   } 
}

fn main() -> anyhow::Result<()> {
    // parse the options
    let opts: Opts = Opts::parse();

    // for now just write some tests to the stdout!

    let stdout = std::io::stdout();
    let lock = stdout.lock();
    let mut buf = std::io::BufWriter::with_capacity(32 * 1024, lock);

    let mut samples: Vec<Sample>;// = Vec::with_capacity(12);
    // //let File1 = Path::new(outpath).join( Path::new(read_file1).file_name());
    // let mut File1 = Path::new(&opts.outpath);

    // let mut File2 = File1.join( Path::new(&opts.read_file1).file_name().unwrap());

    // let file1 = get_writer( &File2 );
    // File2 = File1.join( Path::new(&opts.read_file2).file_name().unwrap() );
    // let file2 = get_writer( &File2);
 
    fs::create_dir_all(&opts.outpath)?;

    if  opts.specie.eq("human") {
        // get all the human sample IDs into this.
        samples = vec![
            Sample::from_description( b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG", 1, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG", 2, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT", 3, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT", 4, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG", 5, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC", 6, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT", 7, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),    
            Sample::from_description( b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT", 8, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG", 9, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC", 10, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC", 11, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG", 12, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 )
        ];

    }
    else if opts.specie.eq("mouse") {
        // and the mouse ones
        samples = vec![
            Sample::from_description( b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", 1, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 2, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", 3, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT", 4, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", 5, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC", 6, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", 7, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),    
            Sample::from_description( b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG", 8, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", 9, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC", 10, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", 11, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
            Sample::from_description( b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC", 12, 9, &opts.outpath, &opts.read_file1, &opts.read_file2 ),
        ];
    } else {
        println!("Sorry, but I have no primers for species {}", &opts.specie);
        std::process::exit(1)
    }

    let file1_path = PathBuf::from(&opts.outpath).join(&opts.read_file1).join("ambig.fq.gz");
    let file2_path = PathBuf::from(&opts.outpath).join(&opts.read_file2).join("ambig.fq.gz");
        
    // need better error handling here too
    let mut file1_ambig_out = GzEncoder::new(File::create(file1_path).unwrap(), Compression::default());
    let mut file2_ambig_out = GzEncoder::new(File::create(file2_path).unwrap(), Compression::default());

    // for now, we're assuming FASTQ and not FASTA.
    let mut reader1 = parse_fastx_file(&opts.read_file1).expect("valid path/file");
    let mut reader2 = parse_fastx_file(&opts.read_file2).expect("valid path/file");

    let mut kmer_vec = Vec::<u64>::with_capacity(12);

    while let Some(record2) = reader2.next() {
        if let Some(record1) = reader1.next() {
            let seqrec = record2.expect("invalid record");
            let seqrec1 = record1.expect("invalid record");
            //let seq = seqrec.seq().into_owned();

            let norm_seq = seqrec.normalize(false);

            // create 9mers of the R2 read and check which of the 12 sample ids matches best:
            let kmers = norm_seq.kmers(9); 
            fill_kmer_vec(kmers, &mut kmer_vec);

            let mut res = vec![0; 12]; //Vec::with_capacity(12);
            let mut maxValue = 0;

            for i in 0..11{
                for s in &kmer_vec {
                    if samples[i].search.contains(&s){
                        res[i] +=1;
                    }
                }
                if res[i] > maxValue{
                    maxValue = res[i];
                }
            }

            let mut z = 0;
            let mut id = 0;
            if maxValue > 2 {
                for i in 0..res.len(){
                    if res[i] == maxValue {
                        id = i;
                        z += z;
                    }
                }
            }

            if z == 1 {
                seqrec1.write(&mut samples[id].file1, None)?;
                seqrec.write(&mut samples[id].file2, None)?;
            } else {
                seqrec1.write(&mut file1_ambig_out, None)?;
                seqrec.write(&mut file2_ambig_out, None)?;
            }

        } else {
            anyhow::bail!("file 2 had reads remaining, but file 1 ran out of reads!");
        }
    }
    Ok(())
}
