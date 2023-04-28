use clap::Parser;

//use this::sampleids::SampleIds;
use this::analysis::MappingInfo;
use this::analysis::Analysis;
//use this::last5::Last5;

use std::path::PathBuf;
use std::fs;

use std::time::SystemTime;

use std::fs::File;

use this::ofiles::Ofiles;


// use std::collections::HashSet;
// use std::convert::TryInto;

/// Just run a test case for speed reproducability and simplicity

#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the input R1 reads file
    #[clap(default_value= "testData/OneSingleCell.11521.R2.fastq.gz",short, long)]
    reads: String,
    /// the input R2 samples file
    #[clap( default_value=  "testData/OneSingleCell.11521.R1.fastq.gz",short, long)]
    file: String,
    /// the specie of the library [mouse, human]
    #[clap(default_value=  "mouse", short, long)]
    specie: String,
    /// the outpath
    #[clap(default_value=  "testData/mapperTest",short, long)]
    outpath: String,
    /// the fasta database containing the genes
    #[clap(default_value= "testData/genes.fasta",short, long)]
    expression: String,
    /// the fasta database containing the antibody tags
    #[clap(default_value= "testData/MyAbSeqPanel.fasta", short, long)]
    antibody: String,
    /// the minimum reads per cell (sample + genes + antibody combined)
    #[clap(default_value_t= 10, short, long)]
    min_umi: usize,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(default_value= "v2.96",short, long)]
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



    // let max_reads = match &opts.max_reads {
    //     Some(val) => val,
    //     None => &usize::MAX,
    // };

    println!("Analysis will stop after having processed {} fastq entries containing a cell info\n", opts.max_reads);


    println!("init models");    

    fs::create_dir_all(&opts.outpath).expect("AlreadyExists");

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

    let ofile = Ofiles::new( 1, "Unknown", "R2.fastq.gz", "R1.fastq.gz",  opts.outpath.as_str() );
    
    // needs log_writer:BufWriter<File>, min_quality:f32, max_reads:usize, ofile:Ofiles
    let mut results = MappingInfo::new( log_file, opts.min_quality, opts.max_reads, ofile );
    

    let pos:&[usize;8];
    let min_sizes:&[usize;2];

    if &opts.version == "v1"{
        pos = &[0,9, 21,30, 43,52, 52,60 ];
        min_sizes = &[ 66, 60 ];
    }
    else {
        pos = &[0,9, 12,22, 26,35 , 36,42 ];
        min_sizes = &[ 51, 51 ];
    }

    // wants gene_kmers:usize, version:String, expression:String, antibody:String, specie:String
    let mut worker = Analysis::new( opts.gene_kmers, opts.version, opts.expression,
        opts.antibody, opts.specie);

    //worker.parse(&mut self,  f1:String, f2:String,  report:&mut MappingInfo,pos: &[usize;8], min_sizes: &[usize;2]  )
    let mut split1 = opts.reads.split(',');
    let mut split2 = opts.file.split(',');
    for f1 in split1.by_ref(){
        if let Some(f2) = split2.next(){
            worker.parse( f1.to_string(), f2.to_string(), &mut results, pos, min_sizes );
        }
        
    }
    

    println!("\n\nWriting outfiles ...");

    worker.write_data( opts.outpath, &results, opts.min_umi );



    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            println!("quantify_rhapsody finished in {milli}h {min}min {sec} sec {mil}milli sec\n" );
            println!("The original implementation took ~0.68 +-0.01 sec to finish");
            println!("  with added genes2 mapping ~0.69 +-0.01 sec to finish");
        },
       Err(e) => {println!("Error: {e:?}");}
    }

}