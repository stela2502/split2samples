/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use crate::traits::CellIndex;
use crate::int_to_str::IntToStr;
use crate::cellids10x::CellId10x;
use crate::genes_mapper::SeqRec;

use std::collections::HashSet;

use std::path::PathBuf;
use std::fs::File;
use flate2::read::GzDecoder;
use std::io::BufReader;
use std::io::BufRead;
use std::env;
use crate::errors::CellIdError;

use std::process;

/// CellData here is a storage for the total UMIs. UMIs will be checked per cell
/// But I do not correct the UMIs here - even with sequencing errors 
/// we should get a relatively correct picture as these folow normal distribution.
pub struct CellIds10x{
    pub ver: String,
    possible: Vec<CellId10x>, // all possible cell tags
    search: HashSet<u32>,
    //detected: HashSet<CellId10x>, // those that we already saw in this analysis
    size: usize, // how long are the sequences normally?
}

impl CellIndex for CellIds10x{

    fn to_cellid (&self, r1: &SeqRec )-> Result<( u32, u64, SeqRec, SeqRec ), CellIdError >{

        let cell = match  r1.slice(0, self.size ){
            Some(c) => c,
            None => return Err(CellIdError::TooShort),
        };
        let umi = match r1.slice( self.size, r1.seq().len() - self.size ){
            Some(c) => c,
            None => return Err(CellIdError::TooShort),
        };

        if let Some(checked_cell) = self.process_sequence( cell.to_u32() ){
            return Ok( (checked_cell, umi.to_u64(), cell, umi ) )
        }

        Err( CellIdError::NoMatch )
    }
}




impl CellIds10x{
    pub fn new(ver:&str )-> Self {
        //println!("Yes you initialized a cellIds10x object!");

        /*
        3M-febrary-2018.txt.gz  Single Cell 3' v3, Single Cell 3' v3.1, Single Cell 3' HT v3.1
        737k-august-2016.txt    Single Cell 3' v2, Single Cell 5' v1 and v2, Single Cell 5' HT v2
        737k-april-2014_rc.txt  Single Cell 3' v1 
        737k-arc-v1.txt.gz  Single Cell Multiome (ATAC+GEX) v1
        737-cratac-v1.txt.gz    Single Cell ATAC 
        (located within cellranger-atac installation directory: cellranger-atac-2.1.0/lib/python/atac/barcodes/)
        9K-LT-march-2021.txt.gz     Single Cell 3' LT
        737k-fixed-rna-profiling.txt.gz     Fixed RNA Profiling (Present starting from Cell Ranger v7.0)
        */
        let filename = match ver{
            "Single Cell 3' v3" => "3M-febrary-2018.txt.gz",
            "Single Cell 3' v3.1" => "3M-febrary-2018.txt.gz",
            "Single Cell 3' HT v3.1" => "3M-febrary-2018.txt.gz",
            "Single Cell 3' v2" => "737k-august-2016.txt.gz",
            "Single Cell 5' v1 and v2" => "737k-august-2016.txt.gz",
            "Single Cell 5' v1" => "737k-august-2016.txt.gz",
            "Single Cell 5' v2" => "737k-august-2016.txt.gz",
            "Single Cell 5' HT v2" => "737k-august-2016.txt.gz",
            "Single Cell 3' v1" => "737k-april-2014_rc.txt.gz",
            "Single Cell Multiome (ATAC+GEX) v1" => "737k-arc-v1.txt.gz",
            "Single Cell ATAC" => "737k-arc-v1.txt.gz",
            //"9K-LT-march-2021.txt.gz" => "Single Cell 3' LT",
            //"Fixed RNA Profiling" => "737k-fixed-rna-profiling.txt.gz",
            _ => {

                let mapping = r#"
                version         =>  file to load
            "Single Cell 3' v3" => "3M-febrary-2018.txt.gz",
            "Single Cell 3' v3.1" => "3M-febrary-2018.txt.gz",
            "Single Cell 3' HT v3.1" => "3M-febrary-2018.txt.gz",
            "Single Cell 3' v2" => "737k-august-2016.txt.gz",
            "Single Cell 5' v1 and v2" => "737k-august-2016.txt.gz",
            "Single Cell 5' v1" => "737k-august-2016.txt.gz",
            "Single Cell 5' v2" => "737k-august-2016.txt.gz",
            "Single Cell 5' HT v2" => "737k-august-2016.txt.gz",
            "Single Cell 3' v1" => "737k-april-2014_rc.txt.gz",
            "Single Cell Multiome (ATAC+GEX) v1" => "737k-arc-v1.txt.gz",
            "Single Cell ATAC" => "737k-arc-v1.txt.gz",
"#;           
                eprintln!("CellId10x does not supprt the 10x version {ver}!");
                eprintln!("All/Only versions supported by CellRanger 6.0 are also available here:\n{mapping}");
                process::exit(0);   
            }
        };

        let mut tool = IntToStr::new(b"AGCT".to_vec(), 16);

        let mut filepath = PathBuf::new();
        if let Ok(path) = env::var("RustodyFiles") {
            filepath.push(path);
            filepath.push(filename);
        } else {
            // Handle case where RustodyFiles environment variable is not set
            eprintln!("RustodyFiles environment variable not set!");
            process::exit(0);   
        }

        let mut possible: Vec<CellId10x> = Vec::new();
        //let mut size: usize;
        #[allow(unused_assignments)]
        let mut size = 0;

        if let Ok(file) = File::open(&filepath) {
            // Check if it's a gz file
            if filepath.extension() == Some("gz".as_ref()) {
                println!("I for a file and I know where it would be");
                let gz = BufReader::new(GzDecoder::new(file));
                size = Self::read_sequences(gz, &mut possible, &mut tool);
            } else {
                size = Self::read_sequences(BufReader::new(file), &mut possible, &mut tool);
            }
        } else {
            // Handle case where file does not exist
            eprintln!("The expected CellRanger CellIDs file {} does not exist!", filepath.display());
            process::exit(0);   
        }

        let search: HashSet<u32> = possible.clone().into_iter().map(|c| c.0).collect();

        Self {
            ver: ver.to_string(),
            possible,
            search,
            size,
        }
    }

    fn process_sequence(&self, sequence: u32) -> Option<u32> {
        // Check if the sequence is already detected
        /*if self.detected.contains(&CellId10x(sequence)) {
            return Some(sequence);
        }*/

        // Check if the sequence matches any possible cell tag

        if self.search.contains( &sequence ){
            //println!("Fast find!");
            return Some( sequence )
        }
        self.possible.iter().find(|&&cell| cell.0 == sequence).map(|cell_id| cell_id.0)
    }

    fn read_sequences<R: BufRead>(reader: R, sequences: &mut Vec<CellId10x>, tool: &mut IntToStr) -> usize {
        let mut size = 0;
        for line in reader.lines() {
            if let Ok(sequence) = line {
                size = sequence.len();
                tool.from_vec_u8( sequence.as_bytes().to_vec() );
                sequences.push( CellId10x( tool.into_u32() ) ) ;
            } else {
                panic!("I could not read a line in the file!")
            }
        }
        size
    }

}



// #[cfg(test)]
// mod tests {
//     #[test]
//     fn getsamples() {
//         let cells = crate::cell_ids::new();

//         let mut primer = b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCT";
//         let mut id:u32 = 1;
//         let mut exp= ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
//         match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
//             Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
//             Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
//         };
        
//         let mut exp2 = vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"];
//         assert_eq!( cells.to_sequence( exp ), exp2 );
//         // 3, 3, 3
//         primer = b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTAT";
//         id = 3;
//         exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
//         exp2 = vec![b"CTTCACATA", b"TGTGAAGAA", b"CACAAGTAT"];
//         match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
//             Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
//             Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
//         };
//         //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
//         assert_eq!( cells.to_sequence( exp ), exp2 );

//         // and the last one
//         primer = b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCA";
//         id = 96;
//         exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
//         exp2 = vec![b"TGCGATCTA", b"CAACAACGG", b"CATAGGTCA"];
//         assert_eq!( 884735+1 , exp);
//         match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
//             Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
//             Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
//         };
//         //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
//         assert_eq!( cells.to_sequence( exp ), exp2 );        
//     }
// }


