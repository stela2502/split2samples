use clap::Parser;

//use this::sampleids::SampleIds;
//use this::analysis::
use bio::io::gff;
use needletail::parse_fastx_file;

use std::path::PathBuf;
use std::fs;

use std::time::SystemTime;

use std::fs::File;
use this::ofiles::Ofiles;
use std::io::BufRead;


// use std::collections::HashSet;
// use std::convert::TryInto;

/// Just run a test case for speed reproducability and simplicity

#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the gff gene information table
    #[clap(default_value= "testData/",short, long)]
    gtf: String,
    /// the fasta genome data
    #[clap( default_value=  "testData/",short, long)]
    file: String,
    /// the outpath
    #[clap(default_value=  "testData/mapperTest",short, long)]
    outpath: String,
    #[clap(default_value= "testData/MyAbSeqPanel.fasta", short, long)]
    antibody: String,
    /// thje mapping kmer length
    #[clap(default_value_t=32, long)]
    gene_kmers: usize,
}

/*
    /// UMI min count - use every umi (per gene; 1) or only reoccuring ones (>1)
    #[clap(default_value_t=1,short, long)]
    umi_count: u8,
*/





// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

    let mut kmer_size = opts.kmer_size;
    if opts.kmer_size > 32{
        println!("Sorry the max size of the kmers is 32 bp");
        kmer_size = 32;
    }

    fs::create_dir_all(&opts.outpath).expect("AlreadyExists");

    let log_file_str = PathBuf::from(&opts.outpath).join(
        "index_log.txt"
    );

    println!( "the log file: {}", log_file_str.file_name().unwrap().to_str().unwrap() );
    
    let log_file = match File::create( log_file_str ){
        Ok(file) => file,
        Err(err) => {
            panic!("Error: {err:#?}" );
        }
    };

    // now we need to read the gtf info and get all genes out of that
    // we specificly need the last 64bp of the end of the gene as that is what we are going to map
    // if there is a splice entry in that area we need to create a spliced and an unspliced entry
    let fp1 = PathBuf::from(opts.gtf);
    let f1 = match File::open(fp1){
        Ok(file) => file,
        Err(err) => panic!("The file {} not exists: {err}", fp1.file_name().unwrap().to_str().unwrap() )
    };
    let file1 = GzDecoder::new( f1 );
    let buff1 = BufReader::new( file1 );

    let mut reader = gff::Reader::new( buff1, gff::GffType::GTF);

    for record in reader.records() {
        let rec = record.ok().expect("Error reading record.");
        println!("{}", rec.seqname());

    }
}