use clap::Parser;

//use this::sampleids::SampleIds;
use rustody::mapping_info::MappingInfo;
use rustody::analysis::Analysis;
//use this::last5::Last5;

use std::path::PathBuf;
use std::fs;

use std::time::SystemTime;

use std::fs::File;

use rustody::ofiles::Ofiles;


// use std::collections::HashSet;
// use std::convert::TryInto;

/// Quantifies a DB Rhapsody experiment and creates sparse matrix outfiles.
/// You need quite long R1 and R2 reads for this! (>70R1 and >70R2 \[v1\] and 52 bp reads for v2.96 and v2.384)

#[derive(Parser)]
#[clap(version = "1.1.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    reads: String,
    /// the input R2 samples file
    #[clap(short, long)]
    file: String,
    /// a pre-defined index folder produced by the cerateIndex scipt
    #[clap(short, long)]
    index: Option<String>,
    /// the specie of the library [mouse, human]
    #[clap(short, long)]
    specie: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// the fasta database containing the genes
    #[clap(short, long)]
    expression: Option<String>,
    /// the fasta database containing the antibody tags
    #[clap(short, long)]
    antibody: Option<String>,
    /// the minimum (UMI) reads per cell (sample + genes + antibody combined)
    #[clap(short, long)]
    min_umi: usize,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(short, long)]
    version: String,
    /// Optional: end the analysis after processing <max_reads> cell fastq entries 
    #[clap(default_value_t=usize::MAX, long)]
    max_reads: usize,
    /// minimal sequencing quality 
    #[clap(default_value_t=25.0, long)]
    min_quality: f32,
    /// minimal sequencing quality 
    #[clap(default_value_t=32, long)]
    gene_kmers: usize,
    /// how many threads to use to analyze this (default max available)
    #[clap(short, long)]
    num_threads: Option<usize>,
    /// this is a BD rhapsody or a 10x expression experiment? 
    #[clap(default_value="bd", long)]
    exp: String,
    /// how many sequences should be analyzed in one chunk
    #[clap(default_value_t=1_000_000, long)]
    chunk_size: usize,
    /// Mapper how many times does a 40bp read combo need to match to any given gene to be reported (default=1)
    #[clap( long)]
    min_matches: Option<usize>,
    /// What is the highest acceptable needleman wush inspired cut off (default 0.5)
    #[clap( long)]
    highest_nw_val: Option<f32>,
    /// What is the highest acceptable humming distance to even run NW (default 0.6)
    #[clap( long)]
    highest_humming_val: Option<f32>,
    /// report the reads matching to a set of genes?
    #[clap(long)]
    report4genes: Option<String>,
}

/*
    /// UMI min count - use every umi (per gene; 1) or only reoccuring ones (>1)
    #[clap(default_value_t=1,short, long)]
    umi_count: u8,
*/

/*// Function to check if a file exists
fn check_file_existence(file_path: &str, option: &str, errors: &mut Vec<String>) {
    if !fs::metadata(file_path).is_ok() {
        errors.push(format!("Option {option} - File not found: {file_path}"));
    }
}*/

// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

    if fs::metadata(&opts.outpath).is_err() {
        if let Err(err) = fs::create_dir_all(&opts.outpath) {
            eprintln!("Error creating directory {}: {}", &opts.outpath, err);
        } else {
            println!("New output directory created successfully!");
        }
    }

    // let max_reads = match &opts.max_reads {
    //     Some(val) => val,
    //     None => &usize::MAX,
    // };

    println!("Analysis will stop after having processed {} fastq entries containing a cell info\n", opts.max_reads);


    println!("init models");    
    let log_file_str = PathBuf::from(&opts.outpath).join(
        "Mapping_log.txt"
    );

    println!( "the log file: {}", log_file_str.file_name().unwrap().to_str().unwrap() );
    
    let log_file = match File::create( log_file_str ){
        Ok(file) => file,
        Err(err) => {
            panic!("Error: {err:#?}" );
        }
    };

    let ofile = Ofiles::new( 1, "Umapped_with_cellID", "R2.fastq.gz", "R1.fastq.gz",  opts.outpath.as_str() );
    
    // needs log_writer:BufWriter<File>, min_quality:f32, max_reads:usize, ofile:Ofiles
    let mut results = MappingInfo::new( Some(log_file), opts.min_quality, opts.max_reads, Some(ofile) );
    

    let pos:&[usize;8];
    let min_sizes:&[usize;2];
    

    if &opts.exp =="10x" {
        min_sizes = &[ 20, 20 ];
        pos = &[0,1, 2,3, 4,5, 6,7 ];
    }
    else if &opts.version == "v1"{
        pos = &[0,9, 21,30, 43,52, 52,60 ];
        min_sizes = &[ 66, 60 ];
    }
    else {
        pos = &[0,9, 13,22, 26,35 , 36,42 ];
        min_sizes = &[ 51, 51 ];
    }


    let num_threads = match opts.num_threads{
        Some(n) => {
            if n < 1{
                1
            } 
            else {
                n
            }
        },
        None => num_cpus::get(),
    };
    
    let save= opts.index.is_none();


    // wants gene_kmers:usize, version:String, expression:String, antibody:String, specie:String
    let mut worker = Analysis::new( opts.gene_kmers, opts.version, opts.expression,
        opts.antibody, opts.specie, opts.index, num_threads, &opts.exp);

    if let Some(genes) = opts.report4genes{
        let slice_str: Vec<&str> = genes.split_whitespace().collect();
        worker.report4gname( &slice_str )
    }
    if let Some(min_matches) = opts.min_matches{
        worker.set_min_matches( min_matches );
        println!("Setting the mapper min_matches to {min_matches}")
    }
    if let Some(highest_nw_val) = opts.highest_nw_val{
       worker.set_highest_nw_val( highest_nw_val );
       println!("Setting the mapper highest_nw_val to {highest_nw_val}")
    }
    if let Some(highest_humming_val) = opts.highest_humming_val{
         worker.set_highest_humming_val( highest_humming_val);
         println!("Setting the mapper highest_humming_val to {highest_humming_val}")
    }


    if save{
        worker.write_index( &opts.outpath );
    }

    let mut split1 = opts.reads.split(',');
    let mut split2 = opts.file.split(',');
    let mut id = 0;
    for f1 in split1.by_ref(){
        if let Some(f2) = split2.next(){
            id += 1;
            println!("\nParsing file pair {id}\n");
            worker.parse_parallel( f1, f2, &mut results, pos, min_sizes, &opts.outpath, opts.max_reads ,opts.chunk_size );
        }
        
    }

    //worker.parse_parallel( &opts.reads, &opts.file, &mut results, pos, min_sizes, &opts.outpath );

    println!("\n\nWriting outfiles ...");

    worker.write_data( opts.outpath, &mut results, opts.min_umi );



    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            println!("quantify_rhapsody finished in {milli}h {min}min {sec} sec {mil}milli sec\n" );},
       Err(e) => {println!("Error: {e:?}");}
    }

}
