use clap::Parser;
use this::cellids::CellIds;

#[derive(Parser)]
#[clap(version = "0.1.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the cell id you want converted
    #[clap(short, long)]
    id: u32,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(short, long)]
    version: String,
}

fn main() {
	let opts: Opts = Opts::parse();

	let cells = CellIds::new(&opts.version);
	let seq: Vec<&[u8; 9]> = cells.to_sequence( opts.id );
	

	println!( "The sequence is:\n{:?}\n{:?}\n{:?}\n", std::str::from_utf8(seq[0]), std::str::from_utf8(seq[1]), std::str::from_utf8(seq[2]) );
}

