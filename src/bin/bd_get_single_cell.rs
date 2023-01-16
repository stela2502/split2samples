use clap::Parser;
use needletail::parse_fastx_file;
use this::cellids::CellIds;

use this::ofiles::Ofiles;
use std::fs;


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
        Err(e) => panic!("I could not create the outpath: {}", e)
    };

    let mut cells = CellIds::new(&opts.version, 9);

    let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
    let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

    let mut ofile = Ofiles::new( opts.id as usize, "OneSingleCell", "R2.fastq.gz", "R1.fastq.gz",  &opts.outpath.as_str() );

    'main: while let Some(record2) = readefile.next() {

        if let Some(record1) = readereads.next() {
            let seqrec = match record2{
                Ok( res ) => res,
                Err(err) => {
                    eprintln!("could not read from R2:\n{}", err);
                    continue 'main;
                }
            };
            let seqrec1 = match record1{
                Ok(res) => res,
                Err(err) => {
                    eprintln!("could not read from R1:\n{}", err);
                    continue 'main;
                }
            };


            if seqrec1.seq().len() < 66 {
                continue;
            }
            if seqrec.seq().len() < 64 {
                continue;
            }

            match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
            	Ok(cell_id) => {
            		if cell_id == opts.id{
            			match seqrec1.write(&mut ofile.buff1, None){
                			Ok(_) => (),
                			Err(err) => println!("{}",err)
            			};
            			match seqrec.write( &mut ofile.buff2, None){
                			Ok(_) => (),
                			Err(err) => println!("{}",err)
            			};
            		}
            	},
            	Err(_) => (),
            };
        }
    }

}