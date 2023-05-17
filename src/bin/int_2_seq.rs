use clap::Parser;
use this::int_to_str::IntToStr;

#[derive(Parser)]
#[clap(version = "0.1.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the cell id you want converted
    #[clap(short, long)]
    int: u128,
    /// the sequence lengs - default to 64bp
    #[clap(default_value_t=64, short, long)]
    length: usize,
}


fn main() {
    let opts: Opts = Opts::parse();

    let tool = IntToStr::new( b"".to_vec(), 32 );

    let mut seq = "".to_string();

    tool.u8_array_to_str( opts.length, opts.int.to_le_bytes().to_vec(), &mut seq );

    println!( "The sequence is:\n{}\n", seq );
}
