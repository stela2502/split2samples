use clap::Parser;

//use this::sampleids::SampleIds;
use rustody::mapping_info::MappingInfo;
use rustody::analysis_genemapper::AnalysisGeneMapper;
use rustody::genes_mapper::sequence_record::SeqRec;
//use this::last5::Last5;

use std::path::PathBuf;
use std::fs;

use std::time::SystemTime;

use std::fs::File;

use rustody::ofiles::Ofiles;


// use std::collections::HashSet;
// use std::convert::TryInto;

/// Uses the logic from the quantify_ scripts to test exactly one read sequence (just the fasta part without the name line)

#[derive(Parser)]
#[clap(version = "1.1.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the dna you want to match
    #[clap(short, long)]
    dna: String,
    /// a pre-defined index folder produced by the cerateIndex scipt
    #[clap(short, long)]
    index: Option<String>,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// the specie of the library [mouse, human]
    #[clap(short, long)]
    specie: String,
    /// the fasta database containing the genes
    #[clap(short, long)]
    expression: Option<String>,
    /// the fasta database containing the antibody tags
    #[clap(short, long)]
    antibody: Option<String>,
    /// how many threads to use to analyze this (default max available)
    #[clap(short, long)]
    num_threads: Option<usize>,
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
    let mut results = MappingInfo::new( Some(log_file), 20.0, 10, Some(ofile) );
    

    let pos = &[0,9, 21,30, 43,52, 52,60 ];


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
    

    // wants gene_kmers:usize, version:String, expression:String, antibody:String, specie:String
    let mut worker = AnalysisGeneMapper::new( 32, "v1".to_string(), opts.expression,
        opts.antibody, opts.specie, opts.index, num_threads, "bd");

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

    worker.debug( Some(true) );
    // data:&[(Vec<u8>, Vec<u8>)], report:&mut MappingInfo, pos: &[usize;8]
    //                                      FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    let r1 = SeqRec::new( b"SomeRead1", b"AGGAGATTAACTGGCCTGCGAGCCTGTTCAGGTAGCGGTGACGACTACATATGCTGCACATTTTTT", b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" );
    let qual: Vec<u8>  = vec![b'F'; opts.dna.len()];
    let r2 = SeqRec::new( b"SomeRead2", &opts.dna.as_bytes(), qual.as_slice() );
    let data= vec![ (r1, r2 )];
    worker.analyze_paralel( &data, &mut results, pos );

    //worker.parse_parallel( &opts.reads, &opts.file, &mut results, pos, min_sizes, &opts.outpath );

    println!("\n\nWriting outfiles ...");

    //worker.write_data( opts.outpath, &mut results, opts.min_umi );



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
