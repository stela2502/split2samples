/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;

use crate::geneids::GeneIds;

use std::io::BufWriter;
use std::fs::File;
use std::io::Write;
use std::fs;

use flate2::Compression;
use flate2::write::GzEncoder;


use std::path::PathBuf;
use std::path::Path;


/// CellData here is a storage for the total UMIs. UMIs will be checked per cell
/// But I do not correct the UMIs here - even with sequencing errors 
/// we should get a relatively correct picture as these folow normal distribution.
pub struct CellData{
    pub kmer_size: usize,
    pub name: std::string::String,
    pub genes: BTreeMap<usize, BTreeMap<u64, u8 >>, // I want to know how many times I got the same UMI
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

    pub fn add(&mut self, geneid: usize, umi:u64 ) -> bool{
        //println!("adding gene id {}", geneid );
        return match self.genes.get_mut( &geneid ) {
            Some( gene ) => {
                match gene.get_mut( &umi ) {
                    Some( count ) => {
                        *count += 1;
                        false
                    },
                    None =>{
                        gene.insert( umi, 1 ); 
                        true
                    }
                }
            }, 
            None => {
                let mut gc:BTreeMap<u64, u8> = BTreeMap::new(); //to store the umis
                gc.insert( umi, 1 );
                self.genes.insert( geneid, gc );
                true
            }
        }
    }

    pub fn n_umi( &self, gene_info:&GeneIds, gnames: &Vec<String>,  min_umi_count:u8 ) -> usize {
        let mut n = 0;

        for name in gnames{
            n += self.n_umi_4_gene( gene_info, name, min_umi_count );
        }
        //println!("I got {} umis for cell {}", n, self.name );
        return n; 
    }

    pub fn n_umi_4_gene( &self, gene_info:&GeneIds, gname:&String, min_umi_count:u8 ) -> usize {
        let mut n = 0;
        let id = match gene_info.names.get( gname ){
            Some(g_id) => g_id,
            None => panic!("I could not resolve the gene name {}", gname ),
        };
        n += match self.genes.get( id  ){
            Some( map ) => {
                let mut h = 0;
                for (_key, value) in map.iter() {
                    if value >= &min_umi_count{
                        h += 1;                        
                    }
                }
                h
            }
            None => 0
        };
        // if n > 0 { println!("I got {} umis for gene {}", n, gname ); }
        return n;
    }

    
    pub fn to_str<'live>(&mut self, gene_info:&GeneIds, names: &Vec<String>, min_umi_count:u8 ) -> String {

        let mut data = Vec::<std::string::String>::with_capacity( gene_info.names.len()+3 ); 
        data.push(self.name.clone());

        // here our internal data already should be stored with the same ids as the gene names.
        let mut total = 0;
        let mut max = 0;
        let mut max_name:std::string::String = "na".to_string();

        for name in names {
            //println!("I collect expression for gene {}", name);
            let n = self.n_umi_4_gene(gene_info, name, min_umi_count );
            if max < n {
                max_name = name.to_string();
                max = n;
            }
            data.push( n.to_string() );
            total += n;
        }

        data.push( max_name.clone() ); // max expressing gene (or sample id in an HTO analysis)
        data.push( (max as f32 / total as f32 ).to_string()); // fraction of reads for the max gene

        let ret = data.join( "\t" );
        format!( "{}",ret)
    }
}



// This SingleCellData needs to copy some of the logics from split2samples - no it actually is totally dufferent
// Here we look for new sample ids and each sample id needs to be a total match to the previousely identified sample id
// Of cause I could also implement something with a whitelist. But that is for the future.
pub struct SingleCellData{    
    kmer_size: usize,
    //kmers: BTreeMap<u64, u32>,
    cells: BTreeMap<u64, CellData>,
    checked: bool,
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

    pub fn write (&mut self, file_path: PathBuf, genes: &mut GeneIds, min_count:usize, min_umi_count:u8) -> Result< (), &str>{

        let mut names: Vec<String> = Vec::with_capacity(genes.names.len());
        for ( name, _id ) in &genes.names {
            names.push( name.to_string() );
        }
        return self.write_sub( file_path, genes, &names, min_count, min_umi_count);
    }

    pub fn write_sub (&mut self, file_path: PathBuf, genes: &mut GeneIds, names: &Vec<String>, min_count:usize, min_umi_count:u8) -> Result< (), &str>{
    
        let rs:bool = Path::new( &file_path.clone() ).exists();
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
            println!("This is questionable - please run mtx_counts before write_sub!");
            self.mtx_counts( genes, names, min_count, min_umi_count );
        }

        for ( _id,  cell_obj ) in &mut self.cells {
            if ! cell_obj.passing {
                failed +=1;
                continue;
            }
            let text = cell_obj.to_str( genes, names,  min_umi_count );
            match writeln!( writer, "{}",text ){
                Ok(_) => passed +=1,
                Err(err) => {
                    eprintln!("write error: {}", err);
                    return Err::<(), &str>("cell data could not be written")   
                }
            };
        }
        println!( "dense matrix: {} cell written - {} cells too view umis", passed, failed );
        Ok( () )
    }


    /// this will create a path and populate that with 10x kind of files.
    pub fn write_sparse (&mut self, file_path: PathBuf, genes: &mut GeneIds, min_count:usize, min_umi_count:u8) -> Result< (), &str>{
        let mut names: Vec<String> = Vec::with_capacity(genes.names.len());
        for ( name, _id ) in &genes.names {
            names.push( name.to_string() );
        }
        return self.write_sparse_sub( file_path, genes, &names, min_count, min_umi_count);
    }

    pub fn write_sparse_sub (&mut self, file_path: PathBuf, genes: &mut GeneIds, names: &Vec<String>, min_count:usize, min_umi_count:u8) -> Result< (), &str>{
            
        let rs = Path::new( &file_path.clone() ).exists();

        let mut passed = 0;
        let mut failed = 0;
        if ! rs {
            match fs::create_dir ( file_path.clone() ){
                Ok(_file) => (),
                Err(err) => {
                     eprintln!("Error?: {:#?}", err);
                 }
            };
        }

        let file = match File::create( file_path.clone().join("matrix.mtx.gz") ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error creating the path?: {:#?}", err);
            }
        };

        let file1 = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::new(file1);

        let file_b = match File::create( file_path.clone().join("barcodes.tsv.gz") ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error creating the path?: {:#?}", err);
            }
        };
        let file2 = GzEncoder::new(file_b, Compression::default());
        let mut writer_b = BufWriter::new(file2);

        match writeln!( writer, "{}\n{}", 
            "%%MatrixMarket matrix coordinate integer general",
             self.mtx_counts( genes, names, min_count, min_umi_count ) ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {}", err);
                return Err::<(), &str>("Header could not be written")
            }
        };

        let file_f = match File::create( file_path.clone().join("features.tsv.gz") ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error creating the path?: {:#?}", err);
            }
        };
        let file3 = GzEncoder::new(file_f, Compression::default());
        let mut writer_f = BufWriter::new(file3);

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
        let mut entries = 0;
        for ( _id,  cell_obj ) in &self.cells {
            if ! cell_obj.passing {
                //println!("failed cell {}", cell_obj.name );
                failed +=1;
                continue;
            }
            passed += 1;
            cell_id += 1;
            match writeln!( writer_b, "{}",cell_obj.name ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {}", err);
                    return Err::<(), &str>("cell barcode could not be written")   
                }
            };

            for (name, gene_id) in &genes.names4sparse {
                //if cell_id == 1{ println!("writing gene info -> Gene {} included in output", name ); }
                let n = cell_obj.n_umi_4_gene( genes, name, min_umi_count );
                if n > 0{
                    match writeln!( writer, "{} {} {}", gene_id, cell_id, n ){
                        Ok(_) => {entries += 1;},
                        Err(err) => {
                            eprintln!("write error: {}", err);
                            return Err::<(), &str>("cell data could not be written")   
                        }   
                    }

                }
            }
        }
        println!( "min UMI count in export function: {}", min_umi_count);
        println!( "sparse Matrix: {} cell(s) and {} gene(s) and {} entries written ({} cells too view umis) to path {:?}; ", passed, genes.names4sparse.len(), entries, failed, file_path.into_os_string().into_string());
        return Ok( () );
    }
    /// Update the gene names for export to sparse
    pub fn update_names_4_sparse( &mut self, genes: &mut GeneIds, names:&Vec<String>, min_umi_count:u8 ) -> [usize; 2] {
        
        let mut entries = 0;
        let mut ncell = 0;
        if ! self.checked{
            panic!("Please always run mtx_counts before update_names_4_sparse");
        }
        genes.names4sparse.clear();
        genes.max_id = 0; // reset to collect the passing genes
        let mut n:usize;
        for ( _id,  cell_obj ) in &mut self.cells {
            if ! cell_obj.passing {
                continue;
            }
            ncell += 1;
            for name in names {
                //if ! genes.names4sparse.contains_key ( name ){
                n = cell_obj.n_umi_4_gene( genes, name, min_umi_count );
                if n > 0{
                    if ! genes.names4sparse.contains_key ( name ){
                        genes.max_id +=1;
                        genes.names4sparse.insert( name.to_string() , genes.max_id );
                        //println!("Gene {} included in output", name );
                    } 
                    entries +=1;
                }
                //}
            }
        }
        if genes.max_id  ==0{
            eprintln!( "None of the genes have data:\n{}", names.join( ", " ) );
        }else {
            println!("{} genes requested and {} with data found", names.len(), genes.max_id);
        }
        if names.len() != genes.max_id{
            // better to run this once more - this does somehow not create the same count if more genes are checked for
            let mut used:Vec<String> = Vec::with_capacity( genes.max_id );
            for (name, gene_id) in &genes.names4sparse {
                used.push(name.to_string());
            }
            return self.update_names_4_sparse(genes, &used, min_umi_count );
        }
        return [  ncell, entries ] ;
    }


    pub fn mtx_counts(&mut self, genes: &mut GeneIds, names: &Vec<String>, min_count:usize, min_umi_count:u8 ) -> String{
        

        if ! self.checked{

            //println!("Checking cell for min umi count!");

            for ( _id,  cell_obj ) in &mut self.cells {
                // total umi check
                let n = cell_obj.n_umi( genes, names, min_umi_count );
                if  n > min_count{
                    cell_obj.passing = true;
                }
            }
            self.checked = true;
            //println!("{} cells have passed the cutoff of {} counts per cell and {} occurances per umi",ncell, min_count, min_umi_count ); 
        }
        
        let ncell_and_entries = self.update_names_4_sparse( genes, names, min_umi_count );

        let ret = format!("{} {} {}", genes.names4sparse.len(), ncell_and_entries[0], ncell_and_entries[1] );
        println!("mtx_counts -> final return: mtx_counts: {}", ret );
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


