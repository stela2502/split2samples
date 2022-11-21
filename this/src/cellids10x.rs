/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;
//use kmers::naive_impl::Kmer;
use std::collections::HashSet;

//use std::thread;
use crate::geneids::GeneIds;
//use crate::geneids::Info;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;
use std::io::BufWriter;
use std::fs::File;
use std::io::Write;

use std::path::PathBuf;


/// CellData here is a storage for the total UMIs. UMIs will be checked per cell
/// But I do not correct the UMIs here - even with sequencing errors 
/// we should get a relatively correct picture as these folow normal distribution.
pub struct CellData{
    pub kmer_size: usize,
    pub name: std::string::String,
    pub genes: BTreeMap<usize, HashSet<u64>>
}

impl CellData{
    pub fn new( kmer_size:usize, name: std::string::String ) -> Self{
        let genes =  BTreeMap::new(); // to collect the sample counts
        let loc_name = name.clone();
        Self{
            kmer_size,
            name: loc_name,
            genes
        }
    }

    pub fn add(&mut self, geneid: usize, umi:u64 ){
        //println!("adding gene id {}", geneid );
        match self.genes.get_mut( &geneid ) {
            Some( gene ) => {
                gene.insert( umi ); // the gene has already been added - check if umi matters
                }, 
            None => {
                let mut gc:HashSet<u64> = HashSet::new(); //to store the umis
                gc.insert( umi );
                self.genes.insert( geneid, gc );
            }
        }
    }
    
    pub fn to_str<'live>(&mut self, gene_info:&GeneIds ) -> std::string::String {

        let mut data = Vec::<std::string::String>::with_capacity( gene_info.names.len()+3 );
        data.push(self.name.clone());

        // here our internal data already should be stored with the same ids as the gene names.
        let mut total = 0;
        let mut max = 0;
        let mut max_name:std::string::String = "na".to_string();

        for (name, id) in &gene_info.names {
            //println!("I collect expression for id {}", id);
            match self.genes.get( id  ){
                Some(hash) => {
                    let n = hash.len();
                    total += n;
                    if  n > max{
                        max = n;
                        max_name = name.clone();
                    }
                    data.push( n.to_string() )
                },
                None => {
                    data.push( 0.to_string() )
                }
            }
        }
        data.push( max_name.clone() ); // max expressing gene (or sample id in an HTO analysis)
        data.push( (max as f32 / total as f32 ).to_string()); // fraction of reads for the max gene

        let ret = data.join( "\t" );
        format!( "{}",ret)
    }
}



// This CellIds10x needs to copy some of the logics from split2samples - no it actually is totally dufferent
// Here we look for new sample ids and each sample id needs to be a total match to the previousely identified sample id
// Of cause I could also implement something with a whitelist. But that is for the future.
pub struct CellIds10x{    
    kmer_size: usize,
    //kmers: BTreeMap<u64, u32>,
    cells: BTreeMap<u64, CellData>
}


// here the functions
impl <'a> CellIds10x{

    pub fn new(kmer_size:usize )-> Self {

        let cells = BTreeMap::new();
        // let kmer_len = 5;
        // TACAGAACA

        Self {
            kmer_size,
            cells
        }
    }

    /// here the get checks for a complete match of the cell ID
    /// and if that fails we need to add
    pub fn get(&mut self, cell_id: u64, name: std::string::String ) -> Result< &mut CellData, &str>{
        
        //println!("CellIDs::get cell_id: {}", cell_id );
        if ! self.cells.contains_key( &cell_id ){
            let data = CellData::new(self.kmer_size, name );
            self.cells.insert( cell_id, data );
        }

        let ret = match self.cells.get_mut(&cell_id){
            Some(c1) => c1, 
            None => return Err::< &mut CellData, &str>("BTreeMap Upstream error")
        };
        Ok( ret )
    }

    pub fn write (&mut self, file_path: PathBuf, genes: &GeneIds) -> Result< (), &str>{
        
        let file = match File::create( file_path ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error: {:#?}", err);
            }
        };
        let mut writer = BufWriter::new(&file);

        match writeln!( writer, "{}", genes.to_header() ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {}", err);
                return Err::<(), &str>("Header could not be written")
            }
        };

        for ( _id,  cell_obj ) in &mut self.cells {

            //println!( "get something here?: {}", cell_obj.to_str( &gene_ids ) );

            match writeln!( writer, "{}", cell_obj.to_str( genes )){
            // the compiler thought this might be more correct...
            //match writeln!( writer, "{}", cell_obj.to_str( <Vec<u64> as Borrow<Borrowed>>::borrow(&gene_ids).clone() ) ){
             Ok(_) => (),
             Err(err) => {
                 eprintln!("write error: {}", err);
                 return Err::<(), &str>("cell data could not be written")   
             }
            }
        }
        Ok( () )
    }
}


// macro_rules! to_cellid {
//    ($r1: expr, $r2: expr) => {
//       //1-9,22-30,44-52
//       $r1.to_cellid($r1, $r2, [0,8], [21,29], [43,51]  )
//    };
// }

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


