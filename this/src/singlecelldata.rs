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
use std::fs;

use std::path::PathBuf;
use std::path::Path;


/// CellData here is a storage for the total UMIs. UMIs will be checked per cell
/// But I do not correct the UMIs here - even with sequencing errors 
/// we should get a relatively correct picture as these folow normal distribution.
pub struct CellData{
    pub kmer_size: usize,
    pub name: std::string::String,
    pub genes: BTreeMap<usize, HashSet<u64>>,
    pub passing: bool // check if this cell is worth exporting. Late game
}

impl CellData{
    pub fn new( kmer_size:usize, name: std::string::String ) -> Self{
        let genes =  BTreeMap::new(); // to collect the sample counts
        let loc_name = name.clone();
        let passing = false;
        Self{
            kmer_size,
            name: loc_name,
            genes,
            passing
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

    pub fn n_umi( &self ) -> usize {
        let mut n = 0;
        for (_id, hash ) in &self.genes {
            n += hash.len();
        }
        return n; 
    }

    
    pub fn to_str<'live>(&mut self, gene_info:&GeneIds, names: &Vec<String>, min_count:usize ) -> Result< String, &str> {

        let mut data = Vec::<std::string::String>::with_capacity( gene_info.names.len()+3 );
        data.push(self.name.clone());

        // here our internal data already should be stored with the same ids as the gene names.
        let mut total = 0;
        let mut max = 0;
        let mut max_name:std::string::String = "na".to_string();
        let mut id: &usize;

        for name in names {
            //println!("I collect expression for gene {}", name);
            id = match gene_info.names.get( &name.to_string() ){
                Some(g_id) => g_id,
                None => return Err::<String, &str>("gene could not be resolved to id")
            };
            match self.genes.get( id  ){
                Some(hash) => {
                    let n = hash.len();
                    total += n;
                    if  n > max{
                        max = n;
                        max_name = name.to_string();
                    }
                    data.push( n.to_string() )
                },
                None => {
                    data.push( 0.to_string() )
                }
            }
        }
        if total < min_count{
            return Err::<String, &str>("not enough data");
        }
        data.push( max_name.clone() ); // max expressing gene (or sample id in an HTO analysis)
        data.push( (max as f32 / total as f32 ).to_string()); // fraction of reads for the max gene

        let ret = data.join( "\t" );
        Ok( format!( "{}",ret) )
    }

}



// This SingleCellData needs to copy some of the logics from split2samples - no it actually is totally dufferent
// Here we look for new sample ids and each sample id needs to be a total match to the previousely identified sample id
// Of cause I could also implement something with a whitelist. But that is for the future.
pub struct SingleCellData{    
    kmer_size: usize,
    //kmers: BTreeMap<u64, u32>,
    cells: BTreeMap<u64, CellData>,
    checked: bool
}


// here the functions
impl <'a> SingleCellData{

    pub fn new(kmer_size:usize )-> Self {

        let cells = BTreeMap::new();
        let checked:bool = false;

        Self {
            kmer_size,
            cells,
            checked
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

    pub fn write (&mut self, file_path: PathBuf, genes: &mut GeneIds, min_count:usize) -> Result< (), &str>{

        let mut names: Vec<String> = Vec::with_capacity(genes.names.len());
        for ( name, _id ) in &genes.names {
            names.push( name.to_string() );
        }
        return self.write_sub( file_path, genes, &names, min_count);
    }

    pub fn write_sub (&mut self, file_path: PathBuf, genes: &mut GeneIds, names: &Vec<String>, min_count:usize) -> Result< (), &str>{

        let mut rs:bool=true;
    
        rs = Path::new( &file_path.clone() ).exists();
        if rs{
            fs::remove_file( &file_path );
        }
        
        let file = match File::create( file_path ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error: {:#?}", err);
            }
        };
        let mut writer = BufWriter::new(&file);

        match writeln!( writer, "{}", genes.to_header_n( names ) ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {}", err);
                return Err::<(), &str>("Header could not be written")
            }
        };

        let mut passed = 0;
        let mut failed = 0;

        if ! self.checked{
            self.mtx_counts( genes, names, min_count );
        }

        for ( _id,  cell_obj ) in &mut self.cells {
            if ! cell_obj.passing {
                //println!("failed cell {}", cell_obj.name );
                failed +=1;
                continue;
            }
            //println!( "get something here?: {}", cell_obj.to_str( &gene_ids ) );
            match cell_obj.to_str( genes, names,  0){
                Ok(text) => match writeln!( writer, "{}",text ){
                    Ok(_) => passed +=1,
                    Err(err) => {
                        eprintln!("write error: {}", err);
                        return Err::<(), &str>("cell data could not be written")   
                    }
                },
                Err(_) => failed +=1 ,
            }
        }
        println!( "dense matrix: {} cell written - {} cells too view umis", passed, failed );
        Ok( () )
    }


    /// this will create a path and populate that with 10x kind of files.
    pub fn write_sparse (&mut self, file_path: PathBuf, genes: &mut GeneIds, min_count:usize) -> Result< (), &str>{
        let mut names: Vec<String> = Vec::with_capacity(genes.names.len());
        for ( name, _id ) in &genes.names {
            names.push( name.to_string() );
        }
        return self.write_sparse_sub( file_path, genes, &names, min_count);
    }

    pub fn write_sparse_sub (&mut self, file_path: PathBuf, genes: &mut GeneIds, names: &Vec<String>, min_count:usize) -> Result< (), &str>{
        
        self.checked = false;

        let mut rs:bool=true;
    
        rs = Path::new( &file_path.clone() ).exists();
        if ! rs {
            match fs::create_dir ( file_path.clone() ){
                Ok(_file) => (),
                Err(err) => {
                     eprintln!("Error?: {:#?}", err);
                 }
            };
        }

        let file = match File::create( file_path.clone().join("matrix.mtx") ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error creating the path?: {:#?}", err);
            }
        };
        let mut writer = BufWriter::new(&file);

        let file_b = match File::create( file_path.clone().join("barcodes.tsv") ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error creating the path?: {:#?}", err);
            }
        };
        let mut writer_b = BufWriter::new(&file_b);

        match writeln!( writer, "{}\n{}", 
            "%%MatrixMarket matrix coordinate integer general",
             self.mtx_counts( genes, names, min_count ) ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {}", err);
                return Err::<(), &str>("Header could not be written")
            }
        };

        let file_f = match File::create( file_path.clone().join("features.tsv") ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error creating the path?: {:#?}", err);
            }
        };
        let mut writer_f = BufWriter::new(&file_f);

        for (name, _id) in &genes.names4sparse {
            match writeln!( writer_f, "{}\t{}\t{}",name, name , "Gene Expression" ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {}", err);
                    return Err::<(), &str>("feature could not be written")   
                }
            }
        }

        let mut cell_id = 0;
        let mut failed = 0;
        let mut gene_id;
        let mut passed = 0;
        let mut entry = 0;
        let mut g_id;
        gene_id = 0;
        for ( _id,  cell_obj ) in &self.cells {
            if ! cell_obj.passing {
                //println!("failed cell {}", cell_obj.name );
                failed +=1;
                continue;
            }
            match writeln!( writer_b, "{}",cell_obj.name ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {}", err);
                    return Err::<(), &str>("cell barcode could not be written")   
                }
            };
            //println!("got the cell {}", cell_obj.name );
            cell_id += 1;
            passed +=1;
            gene_id = 0;
            for (name, _id) in &genes.names4sparse {
            //for name in names{
                g_id = match genes.names.get( &name.to_string() ){
                    Some(g_id) => {
                        //println!("I got the id {} for the gene {} - it will get id {} in the sparse M", g_id, &name, gene_id );
                        g_id
                    },
                    None => return Err::<(), &str>("gene could not be resolved to id")
                };
                gene_id +=1;

                match cell_obj.genes.get(  g_id ){
                    Some(hash) => {
                        //println!("   got the gene {} and id {}", name, gene_id );
                        match writeln!( writer, "{} {} {}", gene_id, cell_id, hash.len() ){
                            Ok(_) => entry +=1,
                            Err(err) => {
                                eprintln!("write error: {}", err);
                                return Err::<(), &str>("cell data could not be written")   
                            }   
                        }
                    },
                    None => ()
                }
            }
        }
        println!( "sparse Matrix: {} cell and {} genes written ({} cells too view umis) to path {:?}; n={}", passed, gene_id, failed,  file_path.into_os_string().into_string(), entry);

        return Ok( () );
    }

    pub fn mtx_counts(&mut self, genes: &mut GeneIds, names: &Vec<String>, min_count:usize ) -> String{
        let mut ncell =0 ;
        let mut nentry=0;
        let mut id:&usize;

        let mut gene_id = 0;

        genes.max_id = 0; // reset to collect the passing genes
        genes.names4sparse.clear();

        if self.checked{
            //eprintln!("The cells have already been checked!");
            for ( _id,  cell_obj ) in &mut self.cells {
                if cell_obj.passing {
                    ncell += 1;
                    gene_id = 0;
                    for name in names {
                        id = match genes.names.get( &name.to_string() ){
                            Some(g_id) => g_id,
                            None => panic!("I could not resolve the gene {}", name ),
                        };
                        match cell_obj.genes.get( id  ){
                            Some(_hash) => {
                                //let n = hash.len();
                                if ! genes.names4sparse.contains_key ( name ){
                                    genes.max_id +=1;
                                    genes.names4sparse.insert( name.to_string() , genes.max_id );
                                }
                                nentry +=1;
                            },
                            None => ()
                        }
                    }
                }
            }
        }
        else  {
            //eprintln!("Checking cell for min umi count!");
            'main: for ( _id,  cell_obj ) in &mut self.cells {
                if cell_obj.n_umi() > min_count{
                    cell_obj.passing = true;
                    ncell += 1;

                    for name in names {
                        gene_id += 1;
                        id = match genes.names.get( &name.to_string() ){
                            Some(g_id) => g_id,
                            None => panic!("I could not resolve the gene {}", name ),
                        };
                        match cell_obj.genes.get( id  ){
                            Some(_hash) => {
                                //let n = hash.len();
                                if ! genes.names4sparse.contains_key ( name ){
                                    genes.max_id +=1;
                                    genes.names4sparse.insert( name.to_string() , genes.max_id );
                                }
                                nentry +=1;
                                
                            },
                            None => ()
                        }
                    }
                }else {
                    continue 'main;
                }
                
            }
            self.checked = true;

            let mut new_names:Vec<String> = Vec::with_capacity( genes.names4sparse.len() );
            for (name, _id ) in &genes.names4sparse{
                new_names.push( name.to_string() );
            }
            //let ret = format!("{} {} {}", genes.names4sparse.len(), ncell, nentry );
            //println!("mtx_counts: {}", ret );
            //println!("restart mtx_counts using only {} genes instead of {} and the genes max id = {}", new_names.len(), names.len(), genes.max_id);
            return self.mtx_counts( genes, &new_names, min_count)
        }
        let ret = format!("{} {} {}", genes.names4sparse.len(), ncell, nentry );
        //println!("mtx_counts: {}", ret );
        return ret;
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


