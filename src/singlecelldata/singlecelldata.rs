/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;
use std::collections::HashSet;


//use crate::geneids::GeneIds;
use crate::fast_mapper::FastMapper;
use crate::mapping_info::MappingInfo;
use crate::singlecelldata::cell_data::GeneUmiHash;
//use crate::singlecelldata::ambient_rna_detect::AmbientRnaDetect;
use crate::singlecelldata::CellData;
//use crate::cellids::CellIds;
use crate::singlecelldata::IndexedGenes;

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
    data: [BTreeMap< u64, CellData>; u8::MAX as usize],
    genes_to_print: Vec::<String>,
    //ambient_cell_content: BTreeMap< u64, CellData>,
    checked: bool,
    passing: usize,
    pub genes_with_data: HashSet<usize>,
    pub num_threads:usize,
    //ambient_store:AmbientRnaDetect,
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

        const EMPTY_MAP: BTreeMap<u64, CellData> = BTreeMap::new();
        let data = [EMPTY_MAP ;u8::MAX as usize];
        let checked:bool = false;
        let passing = 0;
        let genes_with_data = HashSet::new();
        Self {
            //kmer_size,
            data,
            genes_to_print: Vec::<String>::with_capacity(4),
            //ambient_cell_content: BTreeMap::new(),
            checked,
            passing,
            genes_with_data,
            num_threads,
            //ambient_store:AmbientRnaDetect::new(),
        }
    }

    fn to_key(&self, name: &u64 ) -> usize{
        // 48 = 64 -16
        (*name >> 54) as usize
    }

    fn values(&self) -> impl Iterator<Item = &CellData> {
        self.data.iter().flat_map(|map| map.values())
    }

    fn keys(&self) -> Vec<u64> {
        let mut all_keys = Vec::new();
        for map in &self.data {
            all_keys.extend(map.keys().copied());
        }
        all_keys
    }

    pub fn is_empty(&self) -> bool{
        self.data.is_empty()
    }

    pub fn len(&self) -> usize {
        let mut size = 0;
        for map in &self.data {
            size += map.len();
        }
        size
    } 

    fn get_mut(&mut self, key: &u64) -> Option<&mut CellData> {
        let index = self.to_key(key); // Extracting the first u8 of the u64 key
        self.data[index].get_mut(key)
    }

    pub fn get(&self, key: &u64) -> Option<&CellData> {
        let index = self.to_key(key); // Extracting the first u8 of the u64 key
        self.data[index].get(key)
    }

    pub fn keep_only_passing_cells(&mut self) {
        for map in &mut self.data {
            map.retain(|_, cell_data| cell_data.passing);
        }
    }

    /// merge two SingleCellData objects - keep track of the umis!
    pub fn merge( &mut self, other:&SingleCellData ){
        // reset all internal measurements
        self.checked= false;
        self.passing = 0;
        self.genes_with_data = HashSet::new();

        for other_cell in other.values() {
            let index = self.to_key(&other_cell.name); // Extracting the first u8 of the u64 key
            match self.data[index].entry(other_cell.name) {
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    // If cell exists, merge with existing cell
                    let cell = entry.get_mut();
                    cell.merge(other_cell);
                }
                std::collections::btree_map::Entry::Vacant(entry) => {
                    // If cell doesn't exist, insert new cell
                    entry.insert(other_cell.deep_clone());
                }
            }
        }
    }

    /// checks if there are cells that seam to contain data from another cell.
    /// Meaning there are the same gene/umi combinations in both cells
    pub fn find_cell_cell_obverlaps( &mut self ) {

    }


    /// try_insert now is more slick and insterts the data in 256 different sub-areas.
    pub fn try_insert(&mut self, name: &u64, data: GeneUmiHash, report: &mut MappingInfo) -> bool {
        let index = self.to_key(name); // Extracting the first u8 of the u64 key
        self.genes_with_data.insert(data.0);
        self.checked = false;

        let cell_info = self.data[index]
            .entry(*name)
            .or_insert_with(|| CellData::new( *name )); // Insert new cell if not exists

        report.ok_reads += 1;
        if !cell_info.add(data) {
            report.pcr_duplicates += 1;
            report.local_dup += 1;
            false
        } else {
            true
        }
    }

    pub fn write (&mut self, file_path: PathBuf, genes:&IndexedGenes, min_count:usize) -> Result< (), &str>{

        let names = genes.get_all_gene_names();
        return self.write_sub( file_path, genes, &names, min_count);
    }

    pub fn write_sub (&mut self, file_path: PathBuf, genes:&IndexedGenes, names: &Vec<String>, min_count:usize) -> Result< (), &str>{
    
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
            self.update_genes_to_print( genes, names);
            self.mtx_counts( genes, min_count, self.num_threads );
        }

        //println!("We are here exporting a samples table and want these samples to be included: {:?}", names );
        //println!("And we have these ids for them: {:?} using the offset", genes.ids_for_gene_names( names) );

        for cell_obj in self.values() {
            if ! cell_obj.passing{
                continue;
            }
            let text = cell_obj.to_str( genes, names );
            //println!("this should contain some info {}", text);
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
    pub fn write_sparse (&mut self, file_path: PathBuf, genes: &IndexedGenes, min_count:usize) -> Result< (), &str>{
        let names= genes.get_all_gene_names();
        return self.write_sparse_sub( file_path, genes, &names, min_count);
    }

    pub fn write_sparse_sub (&mut self, file_path: PathBuf, genes:&IndexedGenes, names: &Vec<String>, min_count:usize) -> Result< (), &str>{
            
        let rs = Path::new( &file_path ).exists();

        self.update_genes_to_print( genes, names);

        if self.genes_to_print.is_empty(){
            eprintln!("No genes to report on - no data written to path {:?}", file_path.to_str());
            return Ok(())
        }

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
            self.mtx_counts( genes, min_count, self.num_threads ) ){
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

        for name in &self.genes_to_print {
            match writeln!( writer_f, "{name}\t{name}\tGene Expression"  ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {err}" );
                    return Err::<(), &str>("feature could not be written")   
                }
            }
        }

        self.keep_only_passing_cells();

        let mut entries = 0;
        let mut cell_id = 0;

        let gene_ids = genes.ids_for_gene_names( &self.genes_to_print );

        for cell_obj in self.values() {
            if ! cell_obj.passing {
                continue;
            }
            cell_id += 1;

            match writeln!( writer_b, "Cell{}",cell_obj.name ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {err}");
                    return Err::<(), &str>("cell barcode could not be written")   
                }
            };

            for (gene_id, id) in gene_ids.iter().enumerate() {
                let n = cell_obj.n_umi_4_gene_id( id);
                if n > 0 {
                    match writeln!(writer, "{} {cell_id} {n}", gene_id+1) {
                        Ok(_) => { entries += 1; },
                        Err(err) => {
                            eprintln!("write error: {err}");
                            return Err::<(), &str>("cell data could not be written")
                        }
                    }
                }
            }
        } 

        println!( "sparse Matrix: {} cell(s), {} gene(s) and {} entries written to path {:?}; ", cell_id, gene_ids.len(), entries, file_path.into_os_string().into_string());
        Ok( () )
    }
    /// Update the gene names for export to sparse
    /// returns the count of cells and the count of total gene values
    pub fn update_genes_to_print( &mut self, genes:&IndexedGenes, names:&Vec<String>) -> [usize; 2] {
        
        let mut entries = 0;

        self.genes_to_print.clear();

        let cell_keys:Vec<u64> = self.keys();
        let chunk_size = cell_keys.len() / self.num_threads +1;

        let gene_data:Vec<(BTreeMap<std::string::String, usize>, usize)> = cell_keys
        .par_chunks(chunk_size)
        .map(|chunk| {
            // Your parallel processing logic here...
            let mut names4sparse:  BTreeMap::<String, usize> = BTreeMap::new();
            let mut n:usize;
            let mut entry = 0;
            let gene_ids = genes.ids_for_gene_names( names );
            for key in chunk {
                if let Some(cell_obj) = self.get(key){
                    if self.checked & ! cell_obj.passing{
                        println!("ignoring cell");
                        continue;
                    }
                    for (id, int_id) in gene_ids.iter().enumerate() {
                    //for name in names {
                        n = *cell_obj.total_reads.get( int_id  ).unwrap_or(&0);
                        if n > 0{
                            entry +=1;
                            if ! names4sparse.contains_key ( &names[id] ){
                                names4sparse.insert( names[id].to_string() , names4sparse.len() + 1 );
                            } 
                        }
                    }
                }
            }
            (names4sparse, entry)
        }).collect();

        // Merge gene data from different chunks
        let mut names4sparse = HashSet::<String>::new();

        for (genelist, n) in &gene_data{
            entries += n;
            for name in genelist.keys() {
                names4sparse.insert( name.to_string() );
            }
        }

        self.genes_to_print = names4sparse.iter().cloned().collect();


        /*if genes.max_id  ==0 && ! names.is_empty() {
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
        }*/
        [ self.len(), entries ]
    }


    pub fn mtx_counts(&mut self, genes:&IndexedGenes,  min_count:usize, num_threads: usize ) -> String{
        

        if ! self.checked{

            self.passing= 0;

            //println!("Checking cell for min umi count!");

            // here we should firt check if there is some strange similarity over the cells.
            // do we have overlapping gene_id umi connections?
            // should be VERY rare - let's check that!

            // Split keys into chunks and process them in parallel
            let keys = self.keys();
            let chunk_size = keys.len() / num_threads +1; // You need to implement num_threads() based on your requirement

            let results: Vec<(u64, bool)>  = keys
            .par_chunks(chunk_size) 
            .flat_map(|chunk| {
                let min_count = min_count;

                let mut ret= Vec::<(u64, bool)>::with_capacity(chunk_size);
                for key in chunk {
                    if let Some(cell_obj) = &self.get(key) {
                        let n = cell_obj.n_umi( );
                        ret.push( (*key, n >= min_count) );
                    }
                }
                ret
            })
            .collect();

            let mut bad_cells = 0;
            for ( key, passing ) in results{
                if passing{
                    if let Some(cell_obj) = self.get_mut(&key) {
                        cell_obj.passing = passing;
                    }
                }else {
                    bad_cells +=1;
                }
            }

            println!("Dropping cell with too little counts (n={bad_cells})");
            self.keep_only_passing_cells();

            println!("{} cells have passed the cutoff of {} umi counts per cell.\n\n",self.len(), min_count ); 
            self.checked = true;

        }
        
        let ncell_and_entries = self.update_genes_to_print( genes, &self.genes_to_print.clone() );
        self.passing = ncell_and_entries[0];

        let ret = format!("{} {} {}", &self.genes_to_print.len(), ncell_and_entries[0], ncell_and_entries[1] );
        //println!("mtx_counts -> final return: mtx_counts: {}", ret );
        ret
    }

    pub fn n_reads( &mut self, genes:&IndexedGenes, names: &Vec<String> ) -> usize {
        let mut count = 0;

        for cell_obj in self.values() {
            count += cell_obj.n_reads( genes, names )
        }
        count
    }
}






