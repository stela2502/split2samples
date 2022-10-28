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
//use std::io::BufWriter;
use std::str;
//use anyhow::{Context, Result};
mod cellids;
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

struct Sample {
//    oligo:  String,
    id: usize,
    search: BTreeSet<u64>,
    file1:  GzEncoder<File>, 
    file2:  GzEncoder<File>,
    hashSet: HashSet<i32>,
}

impl Sample {
    fn from_description(primer: &[u8], id: usize, sub_len: usize, outpath: &str) -> Self {

        let mut search = BTreeSet::<u64>::new();
        // println!( "split this into {}bp kmers", sub_len );
        // and here comes the fun. Get me all possible sub_len strings in that oligo into the BTreeMap
        for kmer in needletail::kmer::Kmers::new(primer, sub_len.try_into().unwrap() ) {
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            search.insert(km);
        }

        let file1_path = PathBuf::from(outpath).join(format!("{}.R1.fq.gz", id));
        let file2_path = PathBuf::from(outpath).join(format!("{}.R2.fq.gz", id));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let file1 = GzEncoder::new(File::create(file1_path).unwrap(), Compression::default());
        let file2 = GzEncoder::new(File::create(file2_path).unwrap(), Compression::default());

        let hashSet HashSet<i32>= HashSet::with_capacity(10000);

        let id = 0;

        Self {
  //          oligo: String::from_utf8_lossy(primer).into_owned(), // handle errors better
            id,
            search,
            file1,
            file2,
            hashSet,
        }
    }
}


fn fill_kmer_vec<'a>(seq: needletail::kmer::Kmers<'a>, kmer_vec: &mut Vec<u64>) {
   kmer_vec.clear();
   let mut bad = 0;
   for km in seq {
        // I would like to add a try catch here - possibly the '?' works? 
        // if this can not be converted it does not even make sense to keep this info
        for nuc in km{
            if *nuc ==b'N'{
                bad = 1;
            }
        }
        if bad == 0{
            // let s = str::from_utf8(km);
            // println!( "this is the lib: {:?}",  s );
            kmer_vec.push(Kmer::from(km).into_u64());
        }
   } 
}

fn main() -> anyhow::Result<()> {
    // parse the options
    let opts: Opts = Opts::parse();

    // for now just write some tests to the stdout!
    // thanks to Rob I can now skip the stdout part!
    //let stdout = std::io::stdout();
    //let lock = stdout.lock();
    //let buf = std::io::BufWriter::with_capacity(32 * 1024, lock);

    let mut samples: Vec<Sample>;// = Vec::with_capacity(12);
    let sub_len = 9;
    // //let File1 = Path::new(outpath).join( Path::new(reads).file_name());
    // let mut File1 = Path::new(&opts.outpath);

    // let mut File2 = File1.join( Path::new(&opts.reads).file_name().unwrap());

    // let file1 = get_writer( &File2 );
    // File2 = File1.join( Path::new(&opts.file).file_name().unwrap() );
    // let file2 = get_writer( &File2);
 
    fs::create_dir_all(&opts.outpath)?;

    if  opts.specie.eq("human") {
        // get all the human sample IDs into this.
        samples = vec![
            Sample::from_description( b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG", 1, sub_len, &opts.outpath ),
            Sample::from_description( b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG", 2, sub_len, &opts.outpath ),
            Sample::from_description( b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT", 3, sub_len, &opts.outpath ),
            Sample::from_description( b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT", 4, sub_len, &opts.outpath ),
            Sample::from_description( b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG", 5, sub_len, &opts.outpath ),
            Sample::from_description( b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC", 6, sub_len, &opts.outpath ),
            Sample::from_description( b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT", 7, sub_len, &opts.outpath ),    
            Sample::from_description( b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT", 8, sub_len, &opts.outpath ),
            Sample::from_description( b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG", 9, sub_len, &opts.outpath ),
            Sample::from_description( b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC", 10, sub_len, &opts.outpath ),
            Sample::from_description( b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC", 11, sub_len, &opts.outpath ),
            Sample::from_description( b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG", 12, sub_len, &opts.outpath )
        ];

    }
    else if opts.specie.eq("mouse") {
        // and the mouse ones
        samples = vec![
            Sample::from_description( b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", 1, sub_len, &opts.outpath ),
            Sample::from_description( b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 2, sub_len, &opts.outpath ),
            Sample::from_description( b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", 3, sub_len, &opts.outpath ),
            Sample::from_description( b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT", 4, sub_len, &opts.outpath ),
            Sample::from_description( b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", 5, sub_len, &opts.outpath ),
            Sample::from_description( b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC", 6, sub_len, &opts.outpath ),
            Sample::from_description( b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", 7, sub_len, &opts.outpath ),    
            Sample::from_description( b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG", 8, sub_len, &opts.outpath ),
            Sample::from_description( b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", 9, sub_len, &opts.outpath ),
            Sample::from_description( b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC", 10, sub_len, &opts.outpath ),
            Sample::from_description( b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", 11, sub_len, &opts.outpath ),
            Sample::from_description( b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC", 12, sub_len, &opts.outpath ),
        ];
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

            let norm_seq = seqrec.normalize(true);

            // create 9mers of the R2 read and check which of the 12 sample ids matches best:
            // println!( "split R2 into {}bp kmers", sub_len );

            let kmers = norm_seq.kmers(sub_len.try_into().unwrap());
            //let kmers = seqrec.kmers(sub_len.try_into().unwrap());
            fill_kmer_vec(kmers, &mut kmer_vec);

            let mut res = vec![0; 12]; //Vec::with_capacity(12);
            let mut max_value = 0;

            for i in 0..samples.len(){
                let mut id = 0;
                // reduce the number of checks by factor 9. Hope this still leads to ome usable
                // results.
                while id < kmer_vec.len() {
                    if samples[i].search.contains(&kmer_vec[id]){
                        //println!( "kmer matches to sample {}", samples[i].id );
                        res[i] +=1;
                    }
                    id += 9;
                }

                //for s in &kmer_vec {
                //    if samples[i].search.contains(&s){
                //        //println!( "kmer matches to sample {}", samples[i].id );
                //        res[i] +=1;
                //    }
                //}
                if res[i] > max_value{
                    max_value = res[i];
                }
            }

            let mut z = 0;
            let mut id = 0;
            if max_value > 2 {
                for i in 0..res.len(){
                    if res[i] == max_value {
                        id = i;
                        z += 1;
                    }
                }
            }
            //println!( "we have a match to sample {} with a max value of {} and {} samples reaching this value ",id, max_value,z );
            if z == 1 {
                let cellid:i32 = match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                    Ok(val) => val,
                    Err(err) => continue, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
                }
                samples[id].hashSet.insert( cellid ); // will never insert one element twice. Great!

                //println!("MATCH - A read matching to one sample z={} id = {} and max_value={} for this sequence: {}\n   sampleID= {}",
                // z, id, max_value, str::from_utf8(&seqrec.seq())?, cellid );

                samples[id].id += 1;

                //seqrec1.write(&mut samples[id].file1, None)?;
                //seqrec.write(&mut samples[id].file2, None)?;
            } else if z > 1 { 
                println!("REALLY?! A read matching to more than one samples z={} and max_value={}\n{}", z, max_value, str::from_utf8(&seqrec.seq())? );
            }
            else {
                unknown += 1;
                //seqrec1.write(&mut file1_ambig_out, None)?;
                //seqrec.write(&mut file2_ambig_out, None)?;
            }

        } else {
            anyhow::bail!("file 2 had reads remaining, but file 1 ran out of reads!");
        }
    }
    
    println!( "collected sample info:");
    for i in 0..samples.len(){
        println!( "    sample {}: {} reads and {} cells", i+1, samples[i].id, samples[i].hashSet.len() );
        //samples[i].file1.try_finish()?;
        //samples[i].file2.try_finish()?;
    }
    println!(     "genomic reads: {} reads", unknown );
    //file1_ambig_out.try_finish()?;
    //file2_ambig_out.try_finish()?;
    
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
