
use std::io::BufReader;
use flate2::write::GzDecoder;

use std::fs::File;
use std::path::PathBuf;

pub struct Ifiles {
    pub count: u32,
    //pub file1: GzEncoder<File>,
    //pub file2: GzEncoder<File>,
    pub buff1: BufReader<GzDecoder<File>>,
    pub buff2: BufReader<GzDecoder<File>>
}

/// Ofiles encapsulates two BufWriter<GzEncoder<File>> objects to make handling of 20 of these more convenient.
impl Ifiles{
    pub fn new(id: usize, mode:&str, reads:&str, file:&str, outpath:&str )->Self {

        let fp1 = PathBuf::from(reads);
        let fp2 = PathBuf::from(file);
        let file1_path = PathBuf::from(outpath).join(format!("{}.{}.{}", mode, id, fp1.file_name().unwrap().to_str().unwrap() ));
        let file2_path = PathBuf::from(outpath).join(format!("{}.{}.{}", mode, id, fp2.file_name().unwrap().to_str().unwrap() ));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let f1 = match File::open(file1_path){
            Ok(file) => file,
            Err(err) => panic!("The file {mode}.{id}.{} does not exists: {err}", fp1.file_name().unwrap().to_str().unwrap() )
        };
        let f2 = match File::open(file2_path){
            Ok(file) => file,
            Err(err) => panic!("The file {mode}.{id}.{} does not exists: {err}", fp2.file_name().unwrap().to_str().unwrap() )
        };
        
        let file1 = GzDecoder::new(f1);
        let file2 = GzDecoder::new(f2);

        let buff1 = BufReader::new( file1 );
        let buff2 = BufReader::new( file2 );

        let count:u32 = 0;
        Self{
            count,
            //file1,
            //file2,
            buff1,
            buff2
        }
    }
    pub fn close( &mut self ){

        // match self.buff1.close(){
        //     Ok(_) => (),
        //     Err(e) => eprintln!("Could not flush R1: {e}"),
        // };
        // match self.buff2.close(){
        //     Ok(_) => (),
        //     Err(e) => eprintln!("Could not flush R2: {e}"),
        // };
    }
}

pub struct Ifilesr {
    pub count: u32,
    //pub file1: GzEncoder<File>,
    //pub file2: GzEncoder<File>,
    pub buff1: BufReader<File>,
    pub buff2: BufReader<File>
}

/// Ofiles encapsulates two BufWriter<GzEncoder<File>> objects to make handling of 20 of these more convenient.
impl Ifilesr{
    pub fn new(id: usize, mode:&str, reads:&str, file:&str, outpath:&str )->Self {

        let fp1 = PathBuf::from(reads);
        let fp2 = PathBuf::from(file);
        let file1_path = PathBuf::from(outpath).join(format!("{}.{}.{}", mode, id, fp1.file_name().unwrap().to_str().unwrap() ));
        let file2_path = PathBuf::from(outpath).join(format!("{}.{}.{}", mode, id, fp2.file_name().unwrap().to_str().unwrap() ));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let f1 = match File::open(file1_path){
            Ok(file) => file,
            Err(err) => panic!("The file {mode}.{id}.{} does not exists: {err}", fp1.file_name().unwrap().to_str().unwrap() )
        };
        let f2 = match File::open(file2_path){
            Ok(file) => file,
            Err(err) => panic!("The file {mode}.{id}.{} does not exists: {err}", fp2.file_name().unwrap().to_str().unwrap() )
        };
        
        //let file1 = GzDecoder::new(f1);
        //let file2 = GzDecoder::new(f2);

        let buff1 = BufReader::new( f1 );
        let buff2 = BufReader::new( f2 );

        let count:u32 = 0;
        Self{
            count,
            //file1,
            //file2,
            buff1,
            buff2
        }
    }
    pub fn close( &mut self ){

        // match self.buff1.close(){
        //     Ok(_) => (),
        //     Err(e) => eprintln!("Could not flush R1: {e}"),
        // };
        // match self.buff2.close(){
        //     Ok(_) => (),
        //     Err(e) => eprintln!("Could not flush R2: {e}"),
        // };
    }
}