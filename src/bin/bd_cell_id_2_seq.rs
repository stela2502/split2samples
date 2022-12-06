use clap::Parser;
use this::cellids::CellIds;

#[derive(Parser)]
#[clap(version = "0.1.0", author = "Stefan L. <stefan.lang@med.lu.se>, Rob P. <rob@cs.umd.edu>")]
struct Opts {
    /// the cell id you want converted
    #[clap(short, long)]
    id: u32,
}

fn main() {
	let opts: Opts = Opts::parse();

	let cells = CellIds::new();
	let seq: Vec<&[u8; 9]> = cells.to_sequence( opts.id );
	

	println!( "The sequence is:\n{:?}\n{:?}\n{:?}\n", std::str::from_utf8(seq[0]), std::str::from_utf8(seq[1]), std::str::from_utf8(seq[2]) );
}

