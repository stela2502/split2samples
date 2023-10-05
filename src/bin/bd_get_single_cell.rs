use clap::Parser;
use needletail::parse_fastx_file;
use this::cellids::CellIds;

use this::ofiles::Ofiles;
use std::fs;

use indicatif::{ProgressBar, ProgressStyle, MultiProgress};


#[derive(Parser)]
#[clap(version = "0.1.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    reads: String,
    /// the input R2 samples file
    #[clap(short, long)]
    file: String,
    /// the specie of the library [mouse, human]
    #[clap(short, long)]
    id: u32,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(short, long)]
    version: String,
}

fn main() {
    
    let opts: Opts = Opts::parse();

    match fs::create_dir_all(&opts.outpath){
        Ok(_) => (),
        Err(e) => panic!("I could not create the outpath: {e}")
    };

    let cells = CellIds::new(&opts.version, 9);

    let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
    let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

    let mut ofile = Ofiles::new( opts.id as usize, "OneSingleCell", "R2.fastq.gz", "R1.fastq.gz",  opts.outpath.as_str() );

    let pos:Vec<usize>;
    let min_sizes:Vec<usize>;

    if &opts.version == "v1"{
        pos = vec![0,9, 21,30, 43,52, 52,60 ];
        min_sizes = vec![ 66, 60 ];
    }
    else {
        pos = vec![0,9, 12,22, 26,35 , 36,42 ];
        min_sizes = vec![ 51, 51 ];
    }


    let mut count = 0;

    println!("writing all reads from the cell {}", opts.id );

    let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");
    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));
    pb.set_style(spinner_style);
    pb.set_prefix(format!("[{}/?]", 100));

    let split = 1000 * 1000;
    //let split = 1000 ;
    let mut total = 0;

    'main: while let Some(record2) = readefile.next() {
        if let Some(record1) = readereads.next() {

            let seqrec = match record2{
                Ok( res ) => res,
                Err(err) => {
                    eprintln!("could not read from R2:\n{err}");
                    continue 'main;
                }
            };
            let seqrec1 = match record1{
                Ok(res) => res,
                Err(err) => {
                    eprintln!("could not read from R1:\n{err}" );
                    continue 'main;
                }
            };

            if seqrec1.seq().len() < min_sizes[0] {
                continue;
            }
            if seqrec.seq().len() < min_sizes[1] {
                continue;
            }
            //let seq = seqrec.seq().into_owned();
            for nuc in &seqrec1.seq()[pos[6]..pos[7]] {  
                if *nuc ==b'N'{
                    continue 'main;
                }
            }

            //let umi = Kmer::from( &seqrec1.seq()[52..60]).into_u64();
            // let umi = Kmer::from( &seqrec1.seq()[pos[6]..pos[7]]).into_u64();

            // first match the cell id - if that does not work the read is unusable
            //match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
            
        
            if let Ok(cell_id) = cells.to_cellid( &seqrec1.seq(), vec![pos[0],pos[1]], vec![pos[2],pos[3]], vec![pos[4],pos[5]]) {
                total += 1;
                if total % split == 0{
                    let log_str = format!("cell read (any/{}) {total}/{count} here I have cell {cell_id}", opts.id );
                    pb.set_message( log_str.clone() );
                    pb.inc(1);
                }
                if cell_id == opts.id{
                    count += 1;
                    //let mut id = seqrec1.id();
                    //id.push( cell_id.to_u8() )
                    match seqrec1.write(&mut ofile.buff1, None){
                        Ok(_) => (),
                        Err(err) => println!("{err}")
                    };
                    match seqrec.write( &mut ofile.buff2, None){
                        Ok(_) => (),
                        Err(err) => println!("{err}")
                    };
                }
            }
        } 
    }

    println!("I found {count} reads for the cell {}", opts.id );
}