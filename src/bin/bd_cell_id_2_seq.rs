use clap::Parser;
use rustody::cellids::CellIds;
use rustody::int_to_str::IntToStr;

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
	let seq: Vec<u64> = cells.to_sequence( opts.id );

	let mut info = format!("The sequence is:\n");

	let mut s = String::from("");
	let tool = IntToStr::new( b"AAAAAAAAA".to_vec(), 9 );
	tool.u64_to_str(9, &seq[0], &mut s);
	info += &(s.clone() +"---");
	//info += &s.clone();
	s.clear();
	tool.u64_to_str(9, &seq[1], &mut s);
	info += &(s.clone() +"---");
	s.clear();
	//info += &s.clone();
	tool.u64_to_str(9, &seq[2], &mut s);
	info += &s.clone();
	//info += &s.clone();

	println!( "{}", info );
}

