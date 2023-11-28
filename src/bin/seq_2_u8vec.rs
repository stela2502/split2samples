use clap::Parser;
use rustody::int_to_str::IntToStr;

#[derive(Parser)]
#[clap(version = "0.1.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the cell id you want converted
    #[clap(short, long)]
    seq: String,
    /// output format (u8 u16, 32, u64, u128)
    #[clap(default_value="u8", short, long)]
    out: String,
}


fn main() {
    let opts: Opts = Opts::parse();

    let tool = IntToStr::new( opts.seq.as_bytes().to_vec(), 32 );

    match opts.out.as_str(){
    	"u8" => tool.print(),
    	"u16" => println!( "The sequence is:\n{}\n", tool.into_u16( ) ),
    	"u64" => println!( "The sequence is:\n{}\n", tool.into_u64( ) ),
    	//"u128" => println!( "The sequence is:\n{}\n", tool.into_u128( seq) ),
    	_ => panic!("Sorry the out format {} is not supported at the moment", {opts.out} ),
    };

}