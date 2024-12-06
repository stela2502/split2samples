use clap::Parser;

//use regex::Regex;

use rustody::genes_mapper::GenesMapper;
//use rustody::fast_mapper::FastMapper;
//use rustody::gene::Gene;

use needletail::parse_fastx_file;

//use this::sampleids::SampleIds;
//use this::analysis::

use std::path::PathBuf;
use std::fs;
//use std::path::Path;

use std::time::SystemTime;

use std::fs::File;
//use std::io::BufReader;
//use std::io::BufRead;

//use flate2::read::GzDecoder;
//use std::collections::HashMap;
//use std::collections::HashSet;

use rustody::mapping_info::MappingInfo;
use rustody::ofiles::Ofiles;

//use std::thread;
//use rayon::slice::ParallelSlice;
//use rayon::iter::ParallelIterator;

//use indicatif::{ProgressStyle, ProgressBar, MultiProgress};

//static EMPTY_VEC: Vec<String> = Vec::new();

// use std::convert::TryInto;

/// In contrast to the create_index programs this will index the complete fastq entries
/// The results of this would more be like a genomic index - and this is untested!

#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the fasta file (should be gzipped!)
    #[clap(short, long)]
    file: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// create text outfile in addition to the binary [default: false]
    #[clap(default_value_t=false, long)]
    text: bool,
    /// how many threads to use to analyze this (default 1)
    #[clap( short, long)]
    num_threads: Option<usize>,
}

/*
    /// UMI min count - use every umi (per gene; 1) or only reoccuring ones (>1)
    #[clap(default_value_t=1,short, long)]
    umi_count: u8,
*/



// Function to check if a file exists
fn check_file_existence(file_path: &str, option: &str, errors: &mut Vec<String>) {
    if fs::metadata(file_path).is_err() {
        errors.push(format!("Option {option} - File not found: {file_path}"));
    }
}


// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    let now = SystemTime::now();
    
    //// create the report object /////////////////////////////////////
    let opts: Opts = Opts::parse();

    // Check if each file exists
    let mut errors = Vec::new();

    check_file_existence(&opts.file, "file", &mut errors);

    // If there are errors, print them and exit
    if !errors.is_empty() {
        eprintln!("Error: Some files do not exist:");
        for error in &errors {
            eprintln!("{}", error);
        }
        std::process::exit(1);
    }

    if fs::metadata(&opts.outpath).is_ok() {
        if let Err(err) = fs::remove_dir_all(&opts.outpath) {
            eprintln!("Error old index directory: {}", err);
            std::process::exit(1);
        } else {
            println!("Old index directory removed successfully!");
        }

        if let Err(err) = fs::create_dir_all(&opts.outpath) {
            eprintln!("Error creating directory: {}", err);
            std::process::exit(1);
        } else {
            println!("New index directory created successfully!");
        }
    }else if let Err(err) = fs::create_dir_all(&opts.outpath) {
        eprintln!("Error creating directory: {}", err);
        std::process::exit(1);
    } else {
        println!("New index directory created successfully!");
    }

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
    let ofile = Ofiles::new( 1, "NOT_USED", "A.gz", "B.gz",  opts.outpath.as_str() );

    let mut report = MappingInfo::new( Some(log_file), 32.0 , 0, Some(ofile) );
    report.start_counter();

    //// created the report object /////////////////////////////////////

    fs::create_dir_all(&opts.outpath).expect("AlreadyExists");

    // read the fasta data
    let mut expr_file = parse_fastx_file(&opts.file).expect("valid path/file");
    let mut i= 0 ;
    let delimiter: &[u8] = b"  ";  // The byte sequence you want to split by
    println!("Indexing...");
    let mut index = GenesMapper::new( 0 ); // only needs the offset

    let _num_threads = match &opts.num_threads{
        Some(_n) => {
            eprintln!("At the moement we only support 1 thread");
            1
        },
        None => {
            //num_cpus::get(),
            1
        },   
    };

    while let Some(e_record) = expr_file.next() {
        i +=1;
        let seqrec = e_record.expect("invalid record");
        let id = std::str::from_utf8( 
                seqrec.id().split(|&x| x == delimiter[0]).collect::<Vec<&[u8]>>()[0]
            ).unwrap().to_string();

        //id = seqrec.id().split(|&x| x == delimiter[0]).collect()[0];
        let keys = index.add( &seqrec.seq(), &id, &id, &id, 0 );
        eprintln!("\t {keys} keys on chromosome '{}':", id );
    }
    eprintln!();

    // replicate Subreads way to get only unique mappers.
    index.purge( 24 );
    index.make_index_te_ready(); 

    report.stop_single_processor_time();

    let (h,m,s,_ms) = MappingInfo::split_duration( report.absolute_start.elapsed().unwrap() );

    eprintln!("For {} sequence regions we needed {} h {} min and {} sec to process.",i, h, m, s );
    
    //eprintln!("{h} h {m} min {s} sec and {ms} millisec since start");
    eprintln!("{}", report.program_states_string() );

    eprintln!("We created this fast_mapper object:");
    index.print();

    index.write_index( &opts.outpath ).unwrap();

    if opts.text{
        index.write_index_txt( &opts.outpath ).unwrap();
        eprintln!("A text version of the index was written to {} - you can simply remove that after an optional inspection.",opts.outpath );
    }

    //index.write_index_txt( opts.outpath.to_string() ).unwrap();
    //eprintln!("THIS IS STILL IN TEST MODE => TEXT INDEX WRITTEN!!! {}",opts.outpath.to_string() );
    eprintln!("{}", report.program_states_string() );

    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            eprintln!("finished in {milli} h {min} min {sec} sec {mil} milli sec");
        },
        Err(e) => {println!("Error: {e:?}");}
    }

}
