
use std::io::BufWriter;
use flate2::Compression;
use flate2::write::GzEncoder;

use std::fs::File;
use std::path::PathBuf;
use std::io::Write;

pub struct Ofiles {
    pub count: u32,
    //pub file1: GzEncoder<File>,
    //pub file2: GzEncoder<File>,
    pub buff1: BufWriter<GzEncoder<File>>,
    pub buff2: BufWriter<GzEncoder<File>>
}

/// Ofiles encapsulates two BufWriter<GzEncoder<File>> objects to make handling of 20 of these more convenient.
impl Ofiles{
    pub fn new(id: usize, mode:&str, reads:&str, file:&str, outpath:&str )->Self {

        let fp1 = PathBuf::from(reads);
        let fp2 = PathBuf::from(file);
        let file1_path = PathBuf::from(outpath).join(format!("{}.{}.{}", mode, id, fp1.file_name().unwrap().to_str().unwrap() ));
        let file2_path = PathBuf::from(outpath).join(format!("{}.{}.{}", mode, id, fp2.file_name().unwrap().to_str().unwrap() ));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let f1 = match File::create(file1_path){
            Ok(file) => file,
            Err(err) => panic!("The file {mode}.{id}.{} cound not be created: {err}", fp1.file_name().unwrap().to_str().unwrap() )
        };
        let f2 = match File::create(file2_path){
            Ok(file) => file,
            Err(err) => panic!("The file {mode}.{id}.{} cound not be created: {err}", fp2.file_name().unwrap().to_str().unwrap() )
        };
        
        let file1 = GzEncoder::new(f1, Compression::default());
        let file2 = GzEncoder::new(f2, Compression::default());

        let buff1 = BufWriter::new( file1 );
        let buff2 = BufWriter::new( file2 );

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

        match self.buff1.flush(){
            Ok(_) => (),
            Err(e) => eprintln!("Could not flush R1: {e}"),
        };
        match self.buff2.flush(){
            Ok(_) => (),
            Err(e) => eprintln!("Could not flush R2: {e}"),
        };
    }
}

pub struct Ofilesr {
    pub count: u32,
    //pub file1: GzEncoder<File>,
    //pub file2: GzEncoder<File>,
    pub buff1: BufWriter<File>,
    pub buff2: BufWriter<File>
}

/// Ofiles encapsulates two BufWriter<GzEncoder<File>> objects to make handling of 20 of these more convenient.
impl Ofilesr{
    pub fn new(id: usize, mode:&str, reads:&str, file:&str, outpath:&str )->Self {

        let fp1 = PathBuf::from(reads);
        let fp2 = PathBuf::from(file);
        let file1_path = PathBuf::from(outpath).join(format!("{}.{}.{}", mode, id, fp1.file_name().unwrap().to_str().unwrap() ));
        let file2_path = PathBuf::from(outpath).join(format!("{}.{}.{}", mode, id, fp2.file_name().unwrap().to_str().unwrap() ));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let f1 = match File::create(file1_path){
            Ok(file) => file,
            Err(err) => panic!("The file {mode}.{id}.{} cound not be created: {err}", fp1.file_name().unwrap().to_str().unwrap() )
        };
        let f2 = match File::create(file2_path){
            Ok(file) => file,
            Err(err) => panic!("The file {mode}.{id}.{} cound not be created: {err}", fp2.file_name().unwrap().to_str().unwrap() )
        };
        
        //let file1 = GzEncoder::new(f1, Compression::default());
        //let file2 = GzEncoder::new(f2, Compression::default());

        let buff1 = BufWriter::new( f1 );
        let buff2 = BufWriter::new( f2 );

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

        match self.buff1.flush(){
            Ok(_) => (),
            Err(e) => eprintln!("Could not flush R1: {e}"),
        };
        match self.buff2.flush(){
            Ok(_) => (),
            Err(e) => eprintln!("Could not flush R2: {e}"),
        };
    }
}