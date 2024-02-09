use clap::Parser;

//use this::sampleids::SampleIds;
use rustody::mapping_info::MappingInfo;
use rustody::analysis_te::AnalysisTE;
//use this::last5::Last5;

use std::path::PathBuf;
use std::fs;

use std::time::SystemTime;

use std::fs::File;

use rustody::ofiles::Ofiles;
use num_cpus;


// use std::collections::HashSet;
// use std::convert::TryInto;

/// Quantifies a 10x expression or DB Rhapsody experiment and creates sparse matrix outfiles.
/// You need quite long R1 and R2 reads for this! (>70R1 and >70R2 \[v1\] and 52 bp reads for v2.96 and v2.384)

#[derive(Parser)]
#[clap(version = "1.2.5", author = "Stefan L. <stefan.lang@med.lu.se>")]
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
    /// the index build from the TE transcripts
    #[clap(short, long)]
    tfs: Option<String>,
    /// the index containing the 'normal' transcripts
    #[clap(short, long)]
    expression: Option<String>,
    /// the minimum (UMI) reads per cell (sample + genes + expression combined)
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
    /// this is a BD rhapsody (bd) or a 10x expression experiment( default 10x)? 
    #[clap(default_value="10x", long)]
    exp: String,
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

    if ! fs::metadata(&opts.outpath).is_ok() {
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
    


    // wants gene_kmers:usize, version:String, tfs:String, expression:String, specie:String
    let mut worker = AnalysisTE::new( opts.gene_kmers, opts.version, opts.tfs,
        opts.expression, opts.specie, num_threads, &opts.exp);


    let mut split1 = opts.reads.split(',');
    let mut split2 = opts.file.split(',');
    let mut id = 0;
    for f1 in split1.by_ref(){
        if let Some(f2) = split2.next(){
            id += 1;
            println!("\nParsing file pair {id}\n");
            worker.parse_parallel( f1, f2, &mut results, pos, min_sizes, &opts.outpath, opts.max_reads );
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
