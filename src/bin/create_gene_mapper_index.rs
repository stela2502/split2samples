use clap::Parser;

use regex::Regex;

use rustody::genes_mapper::GenesMapper;
//use rustody::fast_mapper::FastMapper;
use rustody::gene::Gene;

use needletail::parse_fastx_file;

//use this::sampleids::SampleIds;
//use this::analysis::

use std::path::PathBuf;
use std::fs;
//use std::path::Path;

use std::time::SystemTime;

use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::collections::HashSet;

use rustody::mapping_info::MappingInfo;
use rustody::ofiles::Ofiles;

use std::thread;
use rayon::slice::ParallelSlice;
use rayon::iter::ParallelIterator;

use indicatif::{ProgressStyle, ProgressBar, MultiProgress};

//static EMPTY_VEC: Vec<String> = Vec::new();

// use std::convert::TryInto;

/// Create a binary indes for the quantify_rhapsody program.
/// Make sure to not cat any gzipped files for Rust.
/// At least here we will only read from the first stream in a catted gz file!

#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the gff/gtf gene information table (ONE gzipped stream - do NOT cat these files!)
    #[clap(short, long)]
    gtf: String,
    /// the fasta genome data
    #[clap(short, long)]
    file: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// the mapping kmer length
    #[clap(default_value_t=32, long)]
    gene_kmers: usize,
    /// create text outfile instead of binary
    #[clap(default_value_t=false, long)]
    text: bool,
    /// how many threads to use to analyze this (default 1)
    #[clap(short, long)]
    num_threads: Option<usize>,
    /// the string to check for gene level names (default gene_name)
    #[clap(default_value="gene_name", long)]
    genename: String,
    /// the string to check for gene id (default gene_id)
    #[clap(default_value="gene_id", long)]
    geneid: String,
    /// the string to check for transcript levels names (transcript_id)
    #[clap(default_value="gene_name", long)]
    transcript: String,
    /// report the expression for a set of genes ("g1 g2 g3")?
    #[clap(long)]
    report4genes: Option<String>,
}

/*
    /// UMI min count - use every umi (per gene; 1) or only reoccuring ones (>1)
    #[clap(default_value_t=1,short, long)]
    umi_count: u8,
*/




fn process_lines ( gtf: &str, re_gene_name: &Regex, 
    re_gene_id: &Regex, re_transcript_id: &Regex) -> HashMap::<String, Gene> {


    if ! gtf.ends_with(".gz") {
        panic!("Please gzip your gtf file - thank you! {}", gtf );
    }

    let f1 = match File::open( gtf ){
        Ok(file) => file,
        Err(err) => panic!("The file {} does not exists: {err}", gtf ),
    };
    let file1 = GzDecoder::new(f1);
    let reader = BufReader::new( file1 );


    let mut genes = HashMap::<String, Gene>::new();
    let mut gene_name:String;
    let mut transcript_id:String;
    

    for rec in reader.lines() {

        let line = rec.expect("Error reading record.");

        let parts: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if parts.len() < 8{
            continue;
        }

        if parts[2] == "transcript"{
            // capture the parts I need using my regexp modules
            if let Some(captures) = re_gene_name.captures( &parts[8].to_string() ){
                gene_name = captures.get(1).unwrap().as_str().to_string();
            }else if let Some(_captures) = re_gene_id.captures( &parts[8].to_string() ){
                continue; // this likely clutters up the data anyhow.
            }
            else {
                panic!("I could not identify a gene_name in the attributes {:?}", &parts[8] );
            }
            
            // if let Some(captures) = re_gene_id.captures( &parts[8].to_string() ){
            //     gene_id = captures.get(1).unwrap().as_str().to_string();
            // }else {
            //     panic!("I could not identify a gene_id in the attributes {:?}", &parts[8].to_string() );
            // }
            if let Some(captures) = re_transcript_id.captures( &parts[8].to_string() ){
                transcript_id = captures.get(1).unwrap().as_str().to_string();
            }else {
                panic!("I could not identify a transcript_id in the attributes {:?}", &parts[8] );
            }
            
            // and add a gene
            // pub fn new(chrom:String, start_s:String, end_s:String, sense_strand_s:String, name:String, id:String )
            let gene = Gene::new( 
                &parts[0],  
                &parts[3], 
                &parts[4], 
                &parts[6],
                &transcript_id, 
                &gene_name, 
                vec![gene_name.to_string(), transcript_id.to_string()] 
            );
            #[cfg(debug_assertions)]
            println!("I am inserting a new gene: {} {}",transcript_id, gene );
            genes.insert( transcript_id, gene );
        }

        if parts[2] == "exon"{
            // capture the parts I need
            //eprintln!("I found an exon!");
            if let Some(captures) = re_transcript_id.captures( &parts[8] ){
                transcript_id = captures.get(1).unwrap().as_str().to_string();
            }else {
                eprintln!("I could not identify a gene_id in the attributes {:?}", &parts[8] );
                continue;
            }
            // and add an exon
            match genes.get_mut( &transcript_id ){
                Some(gene) => gene.add_exon( &parts[3], &parts[4] ),
                None => eprintln!( "ignoring transcript! ({})", transcript_id  )
            }
        }

        // genes are collected from transcripts! to also catch aternative transcription ends.
        // if parts[2] == "gene"{
        //     // here we should start to 'clean out the old ones'!
        //     //eprintln!("We got a gene!");
        //     let start = match parts[3].parse::<usize>(){
        //         Ok(v) => v,
        //         Err(e) => panic!("I could not parse the start of the transcript as usize: {e:?}"),
        //     };
        // }
    }

    genes
}


fn process_genes_multi ( genes: &[Gene], index: &mut GenesMapper ,
    seq_records: &HashMap< String, Vec<u8>>, chr:&Regex, genes2print: &HashSet<String> ){
    const COVERED_AREA:usize = 500; // cover 400 bp of the transcript
    for gene in genes {
        let gene_name: String = gene.name.to_string();
        // Do something with the gene, e.g. remove it
        match seq_records.get( &gene.chrom.to_string() ){
            Some(seq) => {
                gene.add_to_genes_mapper( seq, index, COVERED_AREA, genes2print.contains( &gene_name ) );
                //println!("The genes detected: {:?}", index.names_store );
            },
            None => {
                if chr.is_match ( &gene.chrom.to_string() ){
                    match seq_records.get( &gene.chrom.to_string()[3..] ){
                        Some(seq) => {
                            gene.add_to_genes_mapper( seq, index, COVERED_AREA, genes2print.contains( &gene_name ) );
                            //eprintln!("The genes detected: {:?}", index.names_store );
                        },
                        None => {
                            eprintln!("I do not have the sequence for the chromosome {}", gene.chrom );
                        }
                    }
                }else {
                    match seq_records.get( &format!("chr{}", &gene.chrom ) ){
                        Some(seq) => {
                            gene.add_to_genes_mapper( 
                                seq, index, COVERED_AREA, genes2print.contains( &gene_name ) 
                            );
                            //println!("The genes detected: {:?}", index.names_store );
                        },
                        None => {
                            eprintln!("I do not have the sequence for the chromosome {}", gene.chrom );
                        }
                    }

                }
            }

        }
    }
}


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
    check_file_existence(&opts.gtf, "gtf", &mut errors);
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

    //let kmer_size = opts.gene_kmers.min(32);

    fs::create_dir_all(&opts.outpath).expect("AlreadyExists");

    // let log_file_str = PathBuf::from(&opts.outpath).join(
    //     "index_log.txt"
    // );

    // println!( "the log file: {}", log_file_str.file_name().unwrap().to_str().unwrap() );
    
    // let log_file = match File::create( log_file_str ){
    //     Ok(file) => file,
    //     Err(err) => {
    //         panic!("Error: {err:#?}" );
    //     }
    // };

    let mut genes2print = HashSet::<String>::new();
    if let Some(genes) = opts.report4genes {
        for gene in genes.split_whitespace(){
            genes2print.insert( gene.to_string() );
        }
    }


    // read the fasta data in!
    let mut seq_records = HashMap::< String, Vec<u8>>::new();

    let mut expr_file = parse_fastx_file(&opts.file).expect("valid path/file");
    let mut id: &[u8];
    let delimiter: &[u8] = b"  ";  // The byte sequence you want to split by
    eprint!("sequence names like: ");

    while let Some(e_record) = expr_file.next() {
        let seqrec = e_record.expect("invalid record");
        id = seqrec.id().split(|&x| x == delimiter[0]).collect::<Vec<&[u8]>>()[0];
        //id = seqrec.id().split(|&x| x == delimiter[0]).collect()[0];
        eprint!("'{}', ", std::str::from_utf8(id).unwrap() );
        seq_records.insert( std::str::from_utf8( id ).unwrap().to_string() , seqrec.seq().to_vec());
    }
    eprintln!();



    // and init the Index:

    let mut index = GenesMapper::new( 0 ); // only needs the offset

    // in short I need to get an internal model of a gene to work.
    // I want to know where the gene starts ans ends (likely transcripts)
    // I want to get intron/exon bounds and then add both a spliced version of the transcript 
    // as well as an unspliced version

    /*
    seqname   0
    source    1    
    Feature   2
    start     3
    end       4
    score     5
    strand    6
    frame     7
    attribute 8 
    */
    
    let gtf = Regex::new(r".*gtf.?g?z?$").unwrap();
    let chr = Regex::new(r"^chr").unwrap();

    let re_gene_name: Regex;
    let re_gene_id: Regex;
    let re_transcript_id : Regex;


    match gtf.is_match( &opts.gtf ){
        true => {
            eprintln!("gtf mode");
            let gene_name = format!(r#"{} "([\(\)/\w\d\-\._]*)""#, opts.genename );
            re_gene_name =  Regex::new(&gene_name).unwrap();

            let gene_id = format!(r#"{} "([\(\)/\w\d\-\._]*)""#, opts.geneid );
            re_gene_id =  Regex::new(&gene_id).unwrap();

            let transcript_id = format!(r#"{} "([\(\)/\w\d\-\._]*)""#, opts.transcript );
            re_transcript_id =  Regex::new(&transcript_id).unwrap();
        },
        false => {
            eprintln!("gff mode.");
            let gene_name = format!(r#"{}=([\(\)/\w\d\-\._]*);"#, opts.genename );
            re_gene_name =  Regex::new(&gene_name).unwrap();

            let gene_id = format!(r#"{}=([\(\)/\w\d\-\._]*);"#, opts.geneid );
            re_gene_id =  Regex::new(&gene_id).unwrap();

            let transcript_id = format!(r#"{}=([\(\)/\w\d\-\._]*);"#, opts.transcript );
            re_transcript_id =  Regex::new(&transcript_id).unwrap();
        },
    }
    //let mut gene_id:String;
    

    // now we need to read the gtf info and get all genes out of that
    // we specificly need the last 64bp of the end of the gene as that is what we are going to map
    // if there is a splice entry in that area we need to create a spliced and an unspliced entry

    let num_threads = match &opts.num_threads{
        Some(n) => *n,
        None => num_cpus::get(),
    };
    // this needs to be run like that as the gene info is split over multiple lines. Therefore if we split
    // the gtf data up into random slices we loose the gene transcript exon connections!
    // In the future I could crete the genes and then create the index using multi processor approach.
    // let num_threads = 1;

    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));
    let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");
    pb.set_style(spinner_style);

    let reads_per_chunk = 100_000;
    eprintln!("Starting with data collection");
    
    let max_dim = reads_per_chunk * num_threads;
    //let mut lines = Vec::<String>::with_capacity( max_dim );
    //let mut batch = 0;

    let genes_hash = process_lines( &opts.gtf, &re_gene_name, &re_gene_id, &re_transcript_id );
    let genes: Vec<Gene> = genes_hash.into_values().collect();
    report.stop_file_io_time();
    
    let results:Vec<GenesMapper> = genes.par_chunks(genes.len() / num_threads + 1) // Split the data into chunks for parallel processing
        .map(|data_split| {
            // Get the unique identifier for the current thread
            let _thread_id = thread::current().id();
        
            // Convert the thread ID to a string for use in filenames or identifiers
            //let thread_id_str = format!("{:?}",thread_id );
            //let ofile = Ofiles::new( 1, &("Umapped_with_cellID".to_owned()+&thread_id_str), "R2.fastq.gz", "R1.fastq.gz",  outpath );
            //let log_file_str = PathBuf::from(outpath).join(
            //    format!("Mapping_log_{}.txt",thread_id_str )
            //);
    
            //let log_file = match File::create( log_file_str ){
            //        Ok(file) => file,
            //        Err(err) => {
            //        panic!("thread {thread_id_str} Error: {err:#?}" );
            //    }
            //};
            let mut idx = GenesMapper::new( 0 );

            // Clone or create a new thread-specific report for each task      
            process_genes_multi( data_split, &mut idx, &seq_records, &chr, &genes2print );
            idx

        }) // Analyze each chunk in parallel
    .collect(); // Collect the results into a Vec
    report.stop_multi_processor_time();
    eprintln!("Integrating multicore results");
    for idx in results{
        index.merge(&idx);
        //report.merge( &gex.1 );
    }
    // remove all links that link to more than 24 different locations (rep elements)
    index.purge( 24 );
    index.make_index_te_ready(); 

    report.stop_single_processor_time();

    let (h,m,s,_ms) = MappingInfo::split_duration( report.absolute_start.elapsed().unwrap() );

    eprintln!("For {} sequence regions we needed {} h {} min and {} sec to process.",max_dim, h, m, s );
    
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
