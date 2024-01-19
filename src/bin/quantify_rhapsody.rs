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
#[clap(version = "1.2.4", author = "Stefan L. <stefan.lang@med.lu.se>")]
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
    /// this is a BD rhapsody or a 10x expression experiment? 
    #[clap(default_value="bd", long)]
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
    let mut results = MappingInfo::new( log_file, opts.min_quality, opts.max_reads, ofile );
    

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
    


    let save= opts.index.is_none();


    // wants gene_kmers:usize, version:String, expression:String, antibody:String, specie:String
    let mut worker: Analysis = Analysis::new( opts.gene_kmers, opts.version, opts.expression,
                opts.antibody, opts.specie, opts.index, 1, &opts.exp );


    // match opts.exp.as_str(){
    //     "bd" => {
    //         let mut idx =  Analysis::<CellIds>::new( opts.gene_kmers, opts.version, opts.expression,
    //             opts.antibody, opts.specie, opts.index, 1 );
    //          Box::new( idx )
    //     },
    //     "10x" => {
    //         let mut idx = Analysis::<CellIds10x>::new( opts.gene_kmers, opts.version, opts.expression,
    //             opts.antibody, opts.specie, opts.index, 1 );
    //          Box::new( idx )
    //     },
    //     _ => panic!( "Only the options 'bp' and '10x' are supported for option exp"),
    // };


    if save{
        worker.write_index( &opts.outpath );
    }
    //worker.parse(&mut self,  f1:String, f2:String,  report:&mut MappingInfo,pos: &[usize;8], min_sizes: &[usize;2]  )
    let mut split1 = opts.reads.split(',');
    let mut split2 = opts.file.split(',');
    for f1 in split1.by_ref(){
        if let Some(f2) = split2.next(){
            worker.parse( f1.to_string(), f2.to_string(), &mut results, pos, min_sizes );
        }
        
    }

    println!("\n\nWriting outfiles ...");

    worker.write_data( opts.outpath, &mut results, opts.min_umi ); // added one thread!


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
