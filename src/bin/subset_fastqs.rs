use clap::Parser;

use indicatif::MultiProgress;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use rustody::cellids::CellIds;
use rustody::cellids10x::CellIds10x;
use rustody::traits::CellIndex;
use needletail::parse_fastx_file;

use std::io::BufReader;
use std::io::BufRead;
use std::io;
use std::path::Path;
use std::collections::HashMap;

use std::fs;

use std::time::SystemTime;

use std::fs::File;

use rustody::ofiles::{Ofiles, Fspot};


// use std::collections::HashSet;
// use std::convert::TryInto;

/// Quantifies a DB Rhapsody experiment and creates sparse matrix outfiles.
/// You need quite long R1 and R2 reads for this! (>70R1 and >70R2 \[v1\] and 52 bp reads for v2.96 and v2.384)

#[derive(Parser)]
#[clap(version = "0.0.1", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    reads: String,
    /// the input R2 samples file
    #[clap(short, long)]
    file: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(short, long)]
    version: String,
    /// Optional: end the analysis after writing <max_reads> cell fastq entries 
    #[clap(default_value_t=usize::MAX, long)]
    max_reads: usize,
    /// this is a BD rhapsody or a 10x expression experiment? 
    #[clap(default_value="bd", long)]
    exp: String,
    // A space-separated list of cellID files (one ID per line)
    #[clap(short, long, required = true, use_delimiter = true)]
    id_files: Vec<String>,
}

/*
    /// UMI min count - use every umi (per gene; 1) or only reoccuring ones (>1)
    #[clap(default_value_t=1,short, long)]
    umi_count: u8,
*/


/// Read cell IDs from a file and fill in a HashMap with a u64 cell ID and a usize ID
fn read_cell_ids(cell_id_file: &str, id_map: &mut HashMap<u32, usize>,  id:usize) -> io::Result<()> {
    let file = File::open(cell_id_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        if let Ok(cell_id_str) = line {
            //println!("I have a id {cell_id_str}");
            let cell_id = cell_id_str.trim().parse::<u32>().unwrap();
            //println!("I have a id {cell_id}");
            id_map.insert(cell_id, id);
        }
    }
    Ok(())
}


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
    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));
    let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");
    pb.set_style(spinner_style);


    let mut cell_ids = HashMap::<u32, usize>::new();
    let mut ofiles = Vec::<Ofiles>::with_capacity( opts.id_files.len() );
    let mut id = 0;
    for file in opts.id_files{
        if ! fs::metadata(&file).is_ok() {
            eprintln!("Id file {file} seam to not exists!" );
            std::process::exit(1);
        }
        let base_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
        let _=read_cell_ids( &file, &mut cell_ids, id );
        ofiles.push( Ofiles::new( id, base_name, &opts.reads, &opts.file, &opts.outpath ) );
        id +=1;
    }

    println!("I have a total of {} id(s) to search for.", cell_ids.len());

    let cells:  Box::<dyn CellIndex + Sync> = match opts.exp.as_str() {
        "bd" => Box::new(CellIds::new( &opts.version)),
        "10x" => Box::new(CellIds10x::new( &opts.version )),
        _ => panic!("Only 'bd' or '10x' systems are supported in the exp option")
    };

    let mut r2 = parse_fastx_file(&opts.file).expect("valid path/file"); // R2
    let mut r1 = parse_fastx_file(&opts.reads).expect("valid path/file"); // R1
    let mut lines= 0;
    let mut written = 0;
    let mut filtered = 0;
    'main: while let Some(record1) = r1.next() {
        if let Some(record2) = r2.next() {
            let seqrec1 = match record1{
                Ok(res) => res,
                Err(err) => {
                    eprintln!("could not read from R1:\n{err}");
                    continue 'main;
                }
            };
            let seqrec2 = match record2{
                Ok( res ) => res,
                Err(err) => {
                    eprintln!("could not read from R2:\n{err}");
                    continue 'main;
                }
            };
            lines +=1;

            if lines % 100_000 == 0{
                let log_str = format!( "{:.2} mio entries processed - {} reads written",  lines / 1_000_000, written );
                pb.set_message( log_str.clone() );
                pb.inc(1);
            }
            if written == opts.max_reads{
                eprintln!("Max lines written - finishing");
                break 'main
            }
            if &seqrec2.seq().len() < &52 {
                //println!("read filtered");
                filtered  +=1;
                continue;
            }
            match cells.to_cellid( &seqrec1.seq() ){
                Ok( (cell_id, _umi) ) => {
                    //println!("cell id {cell_id} detected");
                    if let Some(file_id) = cell_ids.get( &cell_id ) {
                        written +=1;
                        let str1 = format!( "@{}\n{}\n+\n{}\n",
                            String::from_utf8_lossy(seqrec1.id()), 
                            String::from_utf8_lossy(&seqrec1.seq()), 
                            String::from_utf8_lossy(seqrec1.qual().unwrap()) );
                        //eprintln!("Is this a fastq entry?! {str1}");
                        
                        let str2 = format!( "@{}\n{}\n+\n{}\n",
                            String::from_utf8_lossy(seqrec2.id()), 
                            String::from_utf8_lossy(&seqrec2.seq()), 
                            String::from_utf8_lossy(seqrec2.qual().unwrap()) );
                        ofiles[*file_id].write_to_oufile( Fspot::Buff2, str2 );
                        ofiles[*file_id].write_to_oufile( Fspot::Buff1, str1 );
                    }else {
                        //eprintln!("I have gotten cellid {cell_id} - not in my list!")
                    }
                }
                Err(_err) => {
                    //eprintln!("I have not gotten a cell id!")
                },
            }
        }
    }

    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            println!("subset_fastqs finished in {milli}h {min}min {sec} sec {mil}milli sec
                written {written}/{lines} fastq reads to {} file pairs\n{filtered} reads filtered\n", ofiles.len() );},
       Err(e) => {println!("Error: {e:?}");}
    }

}
