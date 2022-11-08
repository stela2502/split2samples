use clap::Parser;
//use flate2;
//use flate2::read;
use flate2::Compression;
use flate2::write::GzEncoder;

//use needletail::bitkmer::BitNuclKmer;
//use needletail::{parse_fastx_file, Sequence, FastxReader};
use needletail::{parse_fastx_file, Sequence};
use std::convert::TryInto;

//use serde::Deserialize;
use std::collections::HashSet;
//use std::io::Write;
use std::io::BufWriter;
use std::str;
//use anyhow::{Context, Result};

mod cellids;
mod sampleids;

use crate::cellids::CellIds;
//mod utils;
//use crate::utils::get_all_snps;
//use std::path::Path;
use std::path::PathBuf;
//use std::ffi::OsStr;
use std::fs::File;
use std::fs;

use kmers::naive_impl::Kmer;

use std::collections::BTreeSet;
use std::io::prelude::*;

use byteorder::{BigEndian, ReadBytesExt};


// first, reproduce the appproach from
// https://github.com/jeremymsimon/SPLITseq/blob/main/Preprocess_SPLITseq_collapse_bcSharing.pl

/// Split a pair of BD rhapsody fastq files (R1 and R2) into sample specific fastq pairs
#[derive(Parser)]
#[clap(version = "0.1.0", author = "Stefan L. <stefan.lang@med.lu.se>, Rob P. <rob@cs.umd.edu>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    reads: String,
    /// the input R2 samples file
    #[clap(short, long)]
    file: String,
    /// the specie of the library [mouse, human]
    #[clap(short, long)]
    specie: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
}


fn main() -> anyhow::Result<()> {
    // parse the options
    let opts: Opts = Opts::parse();

    // for now just write some tests to the stdout!
    // thanks to Rob I can now skip the stdout part!
    //let stdout = std::io::stdout();
    //let lock = stdout.lock();
    //let buf = std::io::BufWriter::with_capacity(32 * 1024, lock);

    let sub_len = 9;
    let mut samples: SampleIds::new( sub_len );// = Vec::with_capacity(12);
    // //let File1 = Path::new(outpath).join( Path::new(reads).file_name());
    // let mut File1 = Path::new(&opts.outpath);

    // let mut File2 = File1.join( Path::new(&opts.reads).file_name().unwrap());

    // let file1 = get_writer( &File2 );
    // File2 = File1.join( Path::new(&opts.file).file_name().unwrap() );
    // let file2 = get_writer( &File2);
 
    fs::create_dir_all(&opts.outpath)?;

    let mut cell2sample: Vec<HashSet<u32>>;
    cell2sample = vec![
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
    ];

    if  opts.specie.eq("human") {
        // get all the human sample IDs into this.
        samples.add( b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG" );
        samples.add( b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG" );
        samples.add( b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT" );
        samples.add( b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT" );
        samples.add( b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG" );
        samples.add( b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC" );
        samples.add( b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT" );
        samples.add( b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT" );
        samples.add( b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG" );
        samples.add( b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC" );
        samples.add( b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC" );
        samples.add( b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG" );

    }
    else if opts.specie.eq("mouse") {
        // and the mouse ones
        samples.add( b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG" );
        samples.add( b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC" );
        samples.add( b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC" );
        samples.add( b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT" );
        samples.add( b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT" );
        samples.add( b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC" );
        samples.add( b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC" );
        samples.add( b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG" );
        samples.add( b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC" );
        samples.add( b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC" );
        samples.add( b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT" );
        samples.add( b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC" );

    } else {
        println!("Sorry, but I have no primers for species {}", &opts.specie);
        std::process::exit(1)
    }

    let file1_path = PathBuf::from(&opts.outpath).join("ambig.R1.fq.gz");
    let file2_path = PathBuf::from(&opts.outpath).join("ambig.R2.fq.gz");
        
    // need better error handling here too
    //println!( "why does this file break? {}", file1_path.display() );
    let mut file1_ambig_out = GzEncoder::new(File::create(file1_path).unwrap(), Compression::default());
    let mut file2_ambig_out = GzEncoder::new(File::create(file2_path).unwrap(), Compression::default());

    // for now, we're assuming FASTQ and not FASTA.
    let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
    let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");


    //  now we need to get a CellIDs object, too
    let cells = CellIds::new();

    let mut kmer_vec = Vec::<u64>::with_capacity(60);
    let mut unknown = 0;

    while let Some(record2) = readefile.next() {
        if let Some(record1) = readereads.next() {
            let seqrec = record2.expect("invalid record");
            let seqrec1 = record1.expect("invalid record");
            //let seq = seqrec.seq().into_owned();

            match samples.get( &seqrec.seq() ){
                Ok(id) => {
                    match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                        Ok(val) => cell2sample[id].insert( val ), // will never insert one element twice. Great!
                        Err(err) => {
                            println!("{}",err);
                            continue
                        }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
                    };
                },
                Err(err) => {
                    println!("{}",err);
                    seqrec1.write(&mut file1_ambig_out, None)?;
                    seqrec.write(&mut file2_ambig_out, None)?;
                }
            }
        } else {
            anyhow::bail!("file 2 had reads remaining, but file 1 ran out of reads!");
        }
    }
    
    println!( "collected sample info:");
    for i in 0..cell2sample.len(){
        if cell2sample[i].hash_set.len() > 0 {
            println!( "    sample {}: {} reads and {} cells", i+1, samples.reads[i as u32], cell2sample[i].hash_set.len() );

            let file_path = PathBuf::from(&opts.outpath).join(format!("sample{}.ints.txt",i+1));
            let mut file = File::create( file_path )?;
            let mut writer = BufWriter::new(&file);
            for int in cell2sample[i].drain(){
                //file.write( )
                //BigEndian::write_u32(file, int).unwrap();
                writeln!(writer, "{}", int);
                //println!( "        with sample {}",int );
            }
        }
        // this is not working - I need to first make sure the i32 does work for me here....
        //let file_path = PathBuf::from(&opts.outpath).join("sample{i}.ints.dat");
        //let mut file = File::create( file_path )?;
        //for int in samples[i].hash_set.drain(){
        //    file.write( int )?;
        //}
        
        //samples[i].file1.try_finish()?;
        //samples[i].file2.try_finish()?;
    }
    println!(     "genomic reads: {} reads", unknown );
    file1_ambig_out.try_finish()?;
    file2_ambig_out.try_finish()?;
    
    // now lets rewind the fastq files and actually process the 

    println!("There is more to come here!");


    Ok(())
}



// macro_rules! to_cellid {
//    ($r1: expr) => {
//       //1-9,22-30,44-52
//       to_cellid($r1, [0,8], [21,29], [43,51]  )
//    };
// }
