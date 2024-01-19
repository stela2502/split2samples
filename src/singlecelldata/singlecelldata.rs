/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;
use std::collections::HashSet;


//use crate::geneids::GeneIds;
use crate::fast_mapper::FastMapper;
use crate::mapping_info::MappingInfo;
use crate::singlecelldata::CellData;

use std::io::BufWriter;
use std::fs::File;
use std::io::Write;
use std::fs;

use flate2::Compression;
use flate2::write::GzEncoder;

use std::path::PathBuf;
use std::path::Path;

use rayon::prelude::*;


// This SingleCellData needs to copy some of the logics from split2samples - no it actually is totally different
// Here we look for new sample ids and each sample id needs to be a total match to the previousely identified sample id
// Of cause I could also implement something with a whitelist. But that is for the future.
pub struct SingleCellData{    
    //kmer_size: usize,
    //kmers: BTreeMap<u64, u32>,
    cells: BTreeMap< u64, CellData>,
    checked: bool,
    passing: usize,
    pub genes_with_data: HashSet<usize>,
    pub num_threads:usize,
}

impl Default for SingleCellData {
    fn default() -> Self {
        Self::new(1)
    }
}

// here the functions
impl SingleCellData{

    //pub fn new(kmer_size:usize )-> Self {
    pub fn new(num_threads:usize )-> Self {

        let cells = BTreeMap::new();
        let checked:bool = false;
        let passing = 0;
        let genes_with_data = HashSet::new();
        Self {
            //kmer_size,
            cells,
            checked,
            passing,
            genes_with_data,
            num_threads,
        }
    }

    /// merge two SingleCellData objects - keep track of the umis!
    pub fn merge( &mut self, other:&SingleCellData ){
        // reset all internal measurements
        self.checked= false;
        self.passing = 0;
        self.genes_with_data = HashSet::new();

        for other_cell in other.cells.values() {
            match self.cells.get_mut( &other_cell.name ){
                Some(cell) => { cell.merge( other_cell ) },
                None => {
                    self.cells.insert( other_cell.name, other_cell.deep_clone() );
                }
            }
        }
    }

    /// the old get funtion. Does not work in the new Analysis package. So only to make the tests work again.
    pub fn get(&mut self, name: &u64 ) -> Result< &mut CellData, &str>{
        
        //println!("CellIDs::get cell_id: {}", cell_id );
        self.checked= false;

        self.cells.entry(*name).or_insert_with( || { CellData::new( *name ) });

        let ret = match self.cells.get_mut( name ){
            Some(c1) => c1, 
            None => return Err::< &mut CellData, &str>("BTreeMap Upstream error")
        };
        Ok( ret )
    }


    /// here the get checks for a complete match of the cell ID
    /// and if that fails we need to add
    pub fn try_insert(&mut self, name: &u64, 
        gene_id:&usize,umi:&u64, report: &mut MappingInfo ) ->bool{
        
        //println!("CellIDs::get cell_id: {}", cell_id );
        self.checked= false;
        self.genes_with_data.insert( *gene_id );
        self.cells.entry(*name).or_insert_with( || { CellData::new( *name ) });
        match self.cells.get_mut(&name){
            Some(cell_info) =>  {
                //println!("I got gene {} for cell {}", gene_id , cell_info.name);
                report.ok_reads += 1;
                if ! &cell_info.add( *gene_id, *umi ) {
                    //println!("  -> pcr duplicate");
                    report.pcr_duplicates += 1;
                    report.local_dup += 1;
                    false
                }else {
                    true
                }
                //println!("I got gene {} for cell {} and managed to add this data? {}", gene_id , cell_info.name, ret);
            },
            None => panic!("Could not add a gene expression: gene_id = {gene_id}, umi = {umi}" ),
        }
    }


    /// here the get checks for a complete match of the cell ID
    /// and if that fails we need to add
    pub fn try_insert_multimapper(&mut self, name: &u64, 
        gene_id:&Vec<usize>,umi:&u64, report: &mut MappingInfo ) ->bool{
        
        //println!("CellIDs::get cell_id: {}", cell_id );
        self.checked= false;
        //self.genes_with_data.insert( *gene_id );
        self.cells.entry(*name).or_insert_with( || { CellData::new( *name ) });
        match self.cells.get_mut(&name){
            Some(cell_info) =>  {
                //println!("I got gene {} for cell {}", gene_id , cell_info.name);
                report.ok_reads += 1;
                if ! &cell_info.add_multimapper( gene_id.clone(), *umi ) {
                    //println!("  -> pcr duplicate");
                    report.pcr_duplicates += 1;
                    report.local_dup += 1;
                    false
                }else {
                    true
                }
                //println!("I got gene {} for cell {} and managed to add this data? {}", gene_id , cell_info.name, ret);
            },
            None => panic!("Could not add a gene expression: gene_id = {gene_id:?}, umi = {umi}" ),
        }
    }

    pub fn write (&mut self, file_path: PathBuf, genes: &mut FastMapper, min_count:usize) -> Result< (), &str>{

        let mut names: Vec<String> = Vec::with_capacity(genes.names.len());
        for name in genes.names.keys() {
            names.push( name.to_string() );
        }
        return self.write_sub( file_path, genes, &names, min_count);
    }

    pub fn write_sub (&mut self, file_path: PathBuf, genes: &mut FastMapper, names: &Vec<String>, min_count:usize ) -> Result< (), &str>{
    
        let rs:bool = Path::new( &file_path ).exists();
        if rs && fs::remove_file(  &file_path ).is_ok(){};
        
        let file = match File::create( file_path ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error: {err:#?}");
            }
        };
        let mut writer = BufWriter::new(&file);

        match writeln!( writer, "{}", genes.to_header_n( names ) ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {err}" );
                return Err::<(), &str>("Header could not be written")
            }
        };

        let mut passed = 0;

        if ! self.checked{
            println!("Both mRNA and antibody did not generate data - useless, but exporting all cells!");
            self.mtx_counts( genes, names, min_count, self.num_threads );
        }

        for  cell_obj in self.cells.values() {

            let text = cell_obj.to_str( genes, names );
            match writeln!( writer, "{text}" ){
                Ok(_) => passed +=1,
                Err(err) => {
                    eprintln!("write error: {err}");
                    return Err::<(), &str>("cell data could not be written")   
                }
            };
        }

        println!( "dense matrix: {passed} cell written");
        Ok( () )
    }


    /// this will create a path and populate that with 10x kind of files.
    pub fn write_sparse (&mut self, file_path: PathBuf, genes: &mut FastMapper, min_count:usize) -> Result< (), &str>{
        let names: Vec<String> = genes.names.keys().map(|k| k.to_string()).collect();

        return self.write_sparse_sub( file_path, genes, &names, min_count);
    }

    pub fn write_sparse_sub (&mut self, file_path: PathBuf, genes: &mut FastMapper, names: &Vec<String>, min_count:usize) -> Result< (), &str>{
            
        let rs = Path::new( &file_path ).exists();

        if names.is_empty(){
            eprintln!("No genes to report on - no data written to path {:?}", file_path.to_str());
            return Ok(())
        }

        let mut passed = 0;
        if ! rs {
            match fs::create_dir ( file_path.clone() ){
                Ok(_file) => (),
                Err(err) => {
                     eprintln!("Error?: {err:#?}");
                 }
            };
        }

        if fs::remove_file(file_path.join("matrix.mtx.gz") ).is_ok(){};

        let file = match File::create( file_path.join("matrix.mtx.gz") ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error creating the path?: {err:#?}");
            }
        };
        let file1 = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::with_capacity(4096,file1);


        //rs = Path::new( &file_path.clone().join("barcodes.tsv.gz") );
        //if rs {
        //    fs::remove_file( rs );
        //}
        if  fs::remove_file(file_path.join("barcodes.tsv.gz") ).is_ok(){};

        let file_b = match File::create( file_path.join("barcodes.tsv.gz") ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error creating the path?: {err:#?}");
            }
        };
        let file2 = GzEncoder::new(file_b, Compression::default());
        let mut writer_b = BufWriter::with_capacity(4096,file2);
        match writeln!( writer, "%%MatrixMarket matrix coordinate integer general\n{}", 
             self.mtx_counts( genes, names, min_count, self.num_threads ) ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {err}");
                return Err::<(), &str>("Header could not be written")
            }
        };

        if fs::remove_file(file_path.join("features.tsv.gz") ).is_ok(){};
 
        let file_f = match File::create( file_path.join("features.tsv.gz") ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error creating the path?: {err:#?}");
            }
        };
        let file3 = GzEncoder::new(file_f, Compression::default());
        let mut writer_f = BufWriter::with_capacity(4096,file3);

        for name in genes.names4sparse.keys() {
            match writeln!( writer_f, "{name}\t{name}\tGene Expression"  ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {err}" );
                    return Err::<(), &str>("feature could not be written")   
                }
            }
        }

        let mut entries = 0;
        let mut passing_cells: Vec<&CellData> = Vec::with_capacity( self.passing );
        //let mut bad = 0; // complier warning is incorrect!

        for cell_obj in self.cells.values() {
            if ! cell_obj.passing {
                //println!("bad cell {}", cell_obj.name );
                //bad +=1;
                continue;
            }
            passed += 1;

            match writeln!( writer_b, "Cell{}",cell_obj.name ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {err}");
                    return Err::<(), &str>("cell barcode could not be written")   
                }
            };

            passing_cells.push( cell_obj );
        }

        let mut gene_id = 0;
        for name in genes.names4sparse.keys() {
            gene_id += 1; // the one in the object is crap!
            //if cell_id == 1{ println!("writing gene info -> Gene {} included in output", name ); }
            for (id, cell ) in passing_cells.iter().enumerate().take(self.passing) {
                let n = cell.n_umi_4_gene( genes, name );
                if n > 0{
                    match writeln!( writer, "{gene_id} {} {n}",  id+1){
                        Ok(_) => {entries += 1;},
                        Err(err) => {
                            eprintln!("write error: {err}");
                            return Err::<(), &str>("cell data could not be written")   
                        }   
                    }
                }
            }
        }

        println!( "sparse Matrix: {} cell(s), {} gene(s) and {} entries written to path {:?}; ", passed, genes.names4sparse.len(), entries, file_path.into_os_string().into_string());
        Ok( () )
    }
    /// Update the gene names for export to sparse
    pub fn update_names_4_sparse( &mut self, genes: &mut FastMapper, names:&Vec<String>) -> [usize; 2] {
        
        let mut entries = 0;
        let _ncell = 0;
        if ! self.checked{
            panic!("Please always run mtx_counts before update_names_4_sparse");
        }
        genes.names4sparse.clear();
        genes.max_id = 0; // reset to collect the passing genes

        let cell_keys:Vec<u64> = self.cells.keys().cloned().collect();
        let chunk_size = cell_keys.len() / self.num_threads +1;

        let gene_data:Vec<(BTreeMap<std::string::String, usize>, usize)> = cell_keys
        .par_chunks(chunk_size)
        .map(|chunk| {
            // Your parallel processing logic here...
            let mut names4sparse:  BTreeMap::<String, usize> = BTreeMap::new();
            let mut n:usize;
            let mut entry = 0;
            for key in chunk {
            let cell_obj = self.cells.get(key).unwrap();
                for name in names {
                    n = cell_obj.n_umi_4_gene( genes, name );
                    if n > 0{
                        entry +=1;
                        if ! names4sparse.contains_key ( name ){
                            names4sparse.insert( name.to_string() , names4sparse.len() + 1 );
                        } 
                    }
                }
            }
            (names4sparse, entry)
        }).collect();

        // Merge gene data from different chunks
        for (genelist, n) in &gene_data{
            entries += n;
            for name in genelist.keys() {
                if ! genes.names4sparse.contains_key ( name ){
                    // Insert gene names into the FastMapper
                    genes.max_id +=1;
                    genes.names4sparse.insert( name.to_string() , genes.max_id );
                }
            }
        }

        if genes.max_id  ==0 && ! names.is_empty() {
            let mut to = 10;
            if names.len() < 10{
                to = names.len() -1;
            }
            eprintln!( "None of the genes have data:\n{} ...", names[0..to].join( ", " ) );
        }
        //else { println!("{} genes requested and {} with data found", names.len(), genes.max_id); }
        if names.len() != genes.max_id{
            // better to run this once more - this does somehow not create the same count if more genes are checked for
            let mut used:Vec<String> = Vec::with_capacity( genes.max_id );
            for name in genes.names4sparse.keys() {
                used.push(name.to_string());
            }
            return self.update_names_4_sparse(genes, &used );
        }
        [ self.cells.len(), entries ]
    }


    pub fn mtx_counts(&mut self, genes: &mut FastMapper, names: &Vec<String>, min_count:usize, num_threads: usize ) -> String{
        

        if ! self.checked{

            self.passing= 0;

            //println!("Checking cell for min umi count!");


            // Split keys into chunks and process them in parallel
            let keys: Vec<u64> = self.cells.keys().cloned().collect();
            let chunk_size = keys.len() / num_threads +1; // You need to implement num_threads() based on your requirement


            let results: Vec<(u64, bool)>  = keys
                .par_chunks(chunk_size) 
                .flat_map(|chunk| {
                    let genes = &genes;
                    let min_count = min_count;
  
                    let mut ret= Vec::<(u64, bool)>::with_capacity(chunk_size);
                    for key in chunk {
                        if let Some(cell_obj) = &self.cells.get(key) {
                            let n = cell_obj.n_umi( genes, names );
                            ret.push( (key.clone(), n >= min_count) );
                        }
                    }
                    return ret
                })
                .collect();

            // let mut results = Vec::new();
            // for handle in handles {
            //     let result = handle.join().expect("Thread panicked");
            //     results.extend(result);
            // }

            let mut bad_cells = 0;
            for ( key, passing ) in results{
                if passing{
                    if let Some(cell_obj) = self.cells.get_mut(&key) {
                        cell_obj.passing = passing;
                    }
                }else {
                    bad_cells +=1;
                }
            }

            println!("Dropping cell with too little counts (n={bad_cells})");
            self.cells.retain( |&_key, cell_data| cell_data.passing );

            for cell_obj in self.cells.values_mut(){
                cell_obj.fix_multimappers();
            }
            // for cell_obj in self.cells.values_mut() {
            //     // total umi check
            //     let n = cell_obj.n_umi( genes, names );
            //     if  n >= min_count{
            //         cell_obj.passing = true;
            //     }
            // }
            self.checked = true;
            println!("{} cells have passed the cutoff of {} umi counts per cell",self.cells.len(), min_count ); 
        }
        
        let ncell_and_entries = self.update_names_4_sparse( genes, names );
        self.passing = ncell_and_entries[0];

        let ret = format!("{} {} {}", genes.names4sparse.len(), ncell_and_entries[0], ncell_and_entries[1] );
        //println!("mtx_counts -> final return: mtx_counts: {}", ret );
        ret
    }

    pub fn n_reads( &mut self, genes: &mut FastMapper, names: &Vec<String> ) -> usize {
        let mut count = 0;

        for cell_obj in self.cells.values() {
            count += cell_obj.n_reads( genes, names )
        }
        count
    }
}



// #[cfg(test)]
// mod tests {
//     use crate::singlecelldata::CellData;
//     use crate::singlecelldata::GeneIds;

//      #[test]
//     fn cell_to_str() {
//         let mut cell1 = CellData::new( "Cell1".to_string() );
//         let mut genes = GeneIds::new( 7 );

//         let mut geneid = 0;
        
//         genes.add( b"AGCTGCTAGCCGATATT", "Gene1".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         genes.add( b"CTGTGTAGATACTATAGATAA", "Gene2".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         // the next two should not be in the output
//         genes.add( b"CGCGATCGGATAGCTAGATAGG", "Gene3".to_string() );
//         genes.add( b"CATACAACTACGATCGAATCG", "Gene4".to_string() );

//         geneid = genes.get_id( "Gene1".to_string() );
//         println!("Gene id == {}", geneid);
//         for umi in 0..20 {
//             println!("I add Gene1 ({}) umi {}", geneid, umi );
//             cell1.add( geneid, umi as u64);
//         }
//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..20 {
//             cell1.add( geneid, umi as u64);
//         }



//         //to_str<'live>(&mut self, gene_info:&GeneIds, names: &Vec<String> ) 
//         let names= vec!("Gene1".to_string(), "Gene2".to_string() );
//         let exp2:String = "Cell1\t20\t0\tGene1\t1".to_string();
//         let val = cell1.to_str( &genes, &names).to_string();
//         println!( "{}", val );
//         assert_eq!( val,  exp2 ); 
//     }

//     use crate::singlecelldata::SingleCellData;
//     #[test]
//     fn singlecelldata_to_sparse() {
//         let mut celldata = SingleCellData::new(  );

//         let cell1 = match celldata.get( 1 , "Cell1".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         let mut genes = GeneIds::new( 7 );

//         let mut geneid = 0;
        
//         genes.add( b"AGCTGCTAGCCGATATT", "Gene1".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         genes.add( b"CTGTGTAGATACTATAGATAA", "Gene2".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         // the next two should not be in the output
//         genes.add( b"CGCGATCGGATAGCTAGATAGG", "Gene3".to_string() );
//         genes.add( b"CATACAACTACGATCGAATCG", "Gene4".to_string() );

//         geneid = genes.get_id( "Gene1".to_string() );
//         println!("Gene id == {}", geneid);
//         for umi in 0..20 {
//             println!("I add Gene1 ({}) umi {}", geneid, umi );
//             cell1.add( geneid, umi as u64);
//         }
//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..20 {
//             cell1.add( geneid, umi as u64);
//         }


//         //to_str<'live>(&mut self, gene_info:&GeneIds, names: &Vec<String> ) 
//         let  names= vec!("Gene1".to_string(), "Gene3".to_string() );
//         let  exp2:String = "2 1 2".to_string();
//         let  val = celldata.mtx_counts( &mut genes, &names, 1 as usize );
//         //to_str( &genes, &names, 1 as u8 ).to_string();

//         assert_eq!( val,  exp2 ); 

//     }

//     #[test]
//     fn singlecelldata_to_sparse2cells() {
//         let mut celldata = SingleCellData::new( );

//         let mut cell1 = match celldata.get( 1 , "Cell1".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         let mut genes = GeneIds::new( 7 );

//         let mut geneid = 0;
        
//         genes.add( b"AGCTGCTAGCCGATATT", "Gene1".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         genes.add( b"CTGTGTAGATACTATAGATAA", "Gene2".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         // the next two should not be in the output
//         genes.add( b"CGCGATCGGATAGCTAGATAGG", "Gene3".to_string() );
//         genes.add( b"CATACAACTACGATCGAATCG", "Gene4".to_string() );

//         geneid = genes.get_id( "Gene1".to_string() );
//         println!("Gene id == {}", geneid);
//         for umi in 0..20 {
//             println!("I add Gene1 ({}) umi {}", geneid, umi );
//             cell1.add( geneid, umi as u64);
//         }
//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..20 {
//             cell1.add( geneid, umi as u64);
//         }

//         cell1 = match celldata.get( 2 , "Cell2".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..20 {
//             cell1.add( geneid, umi as u64);
//         }

//         let names = vec!("Gene3".to_string());
//         let val = celldata.mtx_counts( &mut genes, &names, 1 as usize );
//         let exp2 = "1 2 2".to_string();
        
//         assert_eq!( val,  exp2 ); 

//     }
//     use std::fs;
//     use std::fs::File;
//     use std::path::PathBuf;
//     use flate2::read::GzDecoder;
//     use std::io::Read;

//     #[test]
//     fn singlecelldata_to_sparse2cells_outfiles() {
//         let mut celldata = SingleCellData::new( );

//         let mut cell1 = match celldata.get( 1 , "Cell1".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         let mut genes = GeneIds::new( 12 );

//         let mut geneid = 0;
        
//         genes.add( b"AGCTGCTAGCCGATATT", "Gene1".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         genes.add( b"CTGTGTAGATACTATAGATAA", "Gene2".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         // the next two should not be in the output
//         genes.add( b"CGCGATCGGATAGCTAGATAGG", "Gene3".to_string() );
//         genes.add( b"CATACAACTACGATCGAATCG", "Gene4".to_string() );

//         geneid = genes.get_id( "Gene1".to_string() );
//         println!("Gene id == {}", geneid);

//         for umi in 0..20 {
//             //println!("I add Gene1 ({}) umi {}", geneid, umi );
//             assert_eq!( cell1.add( geneid, umi as u64), true);
//         }
//         assert_eq!( cell1.add( geneid, 1 as u64), false);

//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..20 {
//             cell1.add( geneid, umi as u64);
//         }

//         cell1 = match celldata.get( 2 , "Cell2".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..10 {
//             cell1.add( geneid, umi as u64);
//         }

//         let names = vec!("Gene3".to_string());

//         let file_path_sp = PathBuf::from( "../testData/output_sparse");

//         match celldata.write_sparse_sub ( file_path_sp, &mut genes, &names, 1) {
//             Ok(_) => (),
//             Err(err) => panic!("Error in the data write: {}", err)
//         };

//         match fs::create_dir("../testData/output_sparse/"){
//             Ok(_) => (),
//             Err(err) => println!("{}", err),
//         };

//         let file = match File::open("../testData/output_sparse/matrix.mtx.gz"){
//             Ok(f) => f,
//             Err(_err) => panic!("expected outfile is missing"),
//         };
//         let mut gz = GzDecoder::new(file);

//         let val = "%%MatrixMarket matrix coordinate integer general\n1 2 2\n1 1 20\n1 2 10\n".to_string();

//         let mut buffer = String::new();

//         match gz.read_to_string(&mut buffer) {
//             Ok(_string) => {},
//             Err(_err) => {},
//         };

//         assert_eq!( val,  buffer ); 

//     }


// #[test]
//     fn singlecelldata_only_one_call_to_cellid_add() {
//         let mut celldata = SingleCellData::new( );

//         let mut cell1 = match celldata.get( 1 , "Cell1".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         let mut genes = GeneIds::new( 12 );

//         let mut geneid = 0;
        
//         genes.add( b"AGCTGCTAGCCGATATT", "Gene1".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         genes.add( b"CTGTGTAGATACTATAGATAA", "Gene2".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         // the next two should not be in the output
//         genes.add( b"CGCGATCGGATAGCTAGATAGG", "Gene3".to_string() );
//         genes.add( b"CATACAACTACGATCGAATCG", "Gene4".to_string() );

//         geneid = genes.get_id( "Gene1".to_string() );
//         assert_eq!( cell1.add( geneid, 1 as u64), true);

//         geneid = genes.get_id( "Gene3".to_string() );
//         assert_eq!( cell1.add( geneid, 1 as u64), true);

//         cell1 = match celldata.get( 2 , "Cell2".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..10 {
//             cell1.add( geneid, umi as u64);
//         }

//         let names = vec!("Gene1".to_string(), "Gene3".to_string());

//         let file_path_sp = PathBuf::from( "../testData/output_sparse2");

//         match celldata.write_sparse_sub ( file_path_sp, &mut genes, &names, 1) {
//             Ok(_) => (),
//             Err(err) => panic!("Error in the data write: {}", err)
//         };

//         match fs::create_dir("../testData/output_sparse2/"){
//             Ok(_) => (),
//             Err(err) => println!("{}", err),
//         };

//         let file = match File::open("../testData/output_sparse2/matrix.mtx.gz"){
//             Ok(f) => f,
//             Err(_err) => panic!("expected outfile is missing"),
//         };
//         let mut gz = GzDecoder::new(file);

//         let val = "%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 1\n2 1 1\n2 2 10\n".to_string();

//         let mut buffer = String::new();

//         match gz.read_to_string(&mut buffer) {
//             Ok(_string) => {},
//             Err(_err) => {},
//         };

//         assert_eq!( val,  buffer ); 

//     }
// }


