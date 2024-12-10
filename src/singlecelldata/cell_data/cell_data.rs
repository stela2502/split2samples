use std::collections::BTreeMap;
use std::collections::HashSet;

use std::collections::hash_map::DefaultHasher;
use crate::singlecelldata::cell_data::GeneUmiHash;
use crate::singlecelldata::ambient_rna_detect::AmbientRnaDetect;

use std::hash::Hash;
use std::hash::Hasher;

//use crate::geneids::GeneIds;
use crate::singlecelldata::IndexedGenes;
use core::fmt;


/// CellData here is a storage for the total UMIs. UMIs will be checked per cell
/// But I do not correct the UMIs here - even with sequencing errors 
/// we should get a relatively correct picture as these folow normal distribution.
#[derive(Default)]
pub struct CellData{
    //pub kmer_size: usize,
    pub name: u64,
    pub genes: BTreeMap<GeneUmiHash, usize>, // I want to know how many times I got the same UMI
    pub genes_with_data: HashSet<String>, // when exporting genes it is helpfule to know which of the possible genomic genes actually have an expression reported...
    pub total_reads: BTreeMap<usize, usize>, // instead of the per umi counting
    //pub cell_counts: BTreeMap<usize, usize>,
    pub passing: bool, // check if this cell is worth exporting. Late game
    pub total_umis:usize, // brute force adding all different gene types together to speed cell subsetting up
    // Some sequences do not uniquely map to only one gene
    // I want to make sure I find a good way to deal with this problem.
    // I hope that some reads would also map to one of the genes uniquely.
    // This way we could easily add this data there
    pub multimapper: BTreeMap<u64, (HashSet<u64>, Vec<usize>) >,
}



// Implementing Display trait for MapperResult
impl fmt::Display for CellData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "CellData with {} genes and {} total_umis", 
            self.genes.len(), self.total_umis )
    }
}

impl CellData{
    pub fn new(  name: u64 ) -> Self{
        let genes =  BTreeMap::new(); // to collect the sample counts
        let total_reads = BTreeMap::new();
        let genes_with_data = HashSet::new();
        let passing = false;
        let total_umis = 0;
        let multimapper = BTreeMap::new();
        Self{
            name,
            genes,
            genes_with_data,
            total_reads,
            passing,
            total_umis,
            multimapper,
        }
    }

    pub fn split_ambient( &self, ambient:&AmbientRnaDetect ) -> (Self, Self){
        let mut non_ambient_data = Self::new(self.name);
        let mut ambient_data = Self::new(self.name);

        for (gene_hash, count) in &self.genes {
            if ambient.is_ambient(gene_hash) {
                ambient_data.add_count(*gene_hash, *count);
            } else {
                non_ambient_data.add_count(*gene_hash, *count);
            }
        }
        non_ambient_data.passing = true;
        ambient_data.passing = true;
        (non_ambient_data, ambient_data)
    }

    fn _hash_classes(vals: &Vec<usize> ) -> u64 {
        let mut hasher = DefaultHasher::new();
        for val in vals {
            val.hash(&mut hasher);
        }
        // Finalize and return the hash value
        hasher.finish()
    }

    pub fn deep_clone(&self) -> CellData {
        let cloned_genes: BTreeMap<GeneUmiHash , usize> = self
            .genes
            .iter()
            .map(|(k, v)| (*k, *v)) // Clone each HashSet within the BTreeMap
            .collect();

        let cloned_genes_with_data: HashSet<String> = self.genes_with_data.clone();

        let cloned_total_reads: BTreeMap<usize, usize> = self.total_reads.clone();

        let cloned_multimapper:  BTreeMap<u64, (HashSet<u64>, Vec<usize>) > = self.multimapper.clone();

        CellData {
            name: self.name,
            genes: cloned_genes,
            genes_with_data: cloned_genes_with_data,
            total_reads: cloned_total_reads,
            passing: self.passing,
            total_umis: self.total_umis,
            multimapper: cloned_multimapper,
        }
    }

    /// adds the other values into this object
    pub fn merge(&mut self, other: &mut CellData) {
        self.total_umis += other.total_umis;
        let mut too_much = BTreeMap::<usize, usize>::new();

        let other_genes = std::mem::take(&mut other.genes);
        for (gene_umi_combo, counts) in other_genes {
            match self.genes.entry(gene_umi_combo) {
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    *entry.get_mut() += counts;
                    let counter = too_much.entry(gene_umi_combo.0).or_insert(0);
                    *counter += 1;
                }
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(counts);
                }
            }
        }

        let other_total_reads = std::mem::take(&mut other.total_reads);
        for (gene_id, count) in other_total_reads {
            let double_counts = too_much.get(&gene_id).unwrap_or(&0);
            let mine = self.total_reads.entry(gene_id).or_insert(0);
            *mine += count - double_counts;
        }
    }

    /// adds the other values into this object
    pub fn merge_re_id_genes(&mut self, other: &CellData, other_genes: &Vec::<usize> ) {
        self.total_umis += other.total_umis;
        let mut too_much = BTreeMap::<usize, usize>::new();

        for (gene_umi_combo, counts) in &other.genes {
            // re-id the UMI count touple
            let re_ided_umi_combo= GeneUmiHash( 
                other_genes[gene_umi_combo.0],
                gene_umi_combo.1
            );
            match self.genes.entry(re_ided_umi_combo) {
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    *entry.get_mut() += counts;
                    let counter = too_much.entry(re_ided_umi_combo.0).or_insert(0);
                    *counter += 1;
                }
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(*counts);
                }
            }
        }

        //let other_total_reads = std::mem::take(&mut other.total_reads);
        for (gene_id, count) in &other.total_reads {
            let double_counts = too_much.get(&gene_id).copied().unwrap_or(0);
            self.total_reads.entry(*gene_id).and_modify(|e| *e += count - double_counts);
        }
    }

    /* my old function
    pub fn merge(&mut self, other:&CellData ){

        self.total_umis += other.total_umis;
        let mut too_much = BTreeMap::<usize, usize>::new();

        for (gene_umi_combo, counts) in &other.genes {
            match self.genes.get_mut( gene_umi_combo ) {
                Some(count) => {
                    *count += counts;
                    let counter = too_much.entry(gene_umi_combo.0).or_insert(0);
                    *counter += 1;
                },
                None => {
                    self.genes.insert( *gene_umi_combo, 1);
                }
            }
        }
        for (gene_id, count) in &other.total_reads {
            let double_counts = *too_much.entry(*gene_id).or_insert_with(|| 0);
            let mine = self.total_reads.entry(*gene_id).or_insert(0);
            *mine += count - double_counts;
        }
        /*for (hash, (umis, ids) ) in &other.multimapper {
            let double_counts = 0;
            for umi in umis {
                self.add_multimapper( ids.clone(), *umi );
                double_counts += 1;
            }
            match self.total_reads.get_mut( &(*hash as usize) ){
                Some( val ) => {
                    let add = other.total_reads.get( &(*hash as usize) ).unwrap_or(&double_counts) - double_counts;
                    *val += add
                },
                None => panic!("merge could not copy the total gene values for multimapper hash {hash}",  ),
            };
        }*/
    }*/

    // returns false if the gene/umi combo has already been recorded!
    pub fn add(&mut self, gh: GeneUmiHash ) -> bool{
        //println!("adding gene id {}", geneid );
        return match self.genes.get_mut( &gh ) {
            Some( gene ) => {
                *gene +=1;
                false
            }, 
            None => {
                self.genes.insert(gh, 1);
                let counts = self.total_reads.entry( gh.0 ).or_insert_with(|| 0);
                *counts +=1;
                self.total_umis += 1;
                true
            }
        }
    }

    // returns false if the gene/umi combo has already been recorded!
    pub fn add_count(&mut self, gh: GeneUmiHash, count:usize ) -> bool{
        //println!("adding gene id {}", geneid );
        return match self.genes.get_mut( &gh ) {
            Some( gene ) => {
                *gene += count;
                false
            }, 
            None => {
                self.genes.insert(gh, count);
                let counts = self.total_reads.entry( gh.0 ).or_insert_with(|| 0);
                *counts +=1;
                self.total_umis += 1;
                true
            }
        }
    }

    /*
    // returns false if the gene/umi combo has already been recorded!
    pub fn add_multimapper(&mut self, geneids: Vec::<usize>, umi:u64 ) -> bool{
        //println!("adding gene id {}", geneid );
        //    pub multimapper: BTreeMap<u64, (HashSet<u64>, Vec<usize>) >,

        let hash = Self::hash_classes( &geneids );

        return match self.multimapper.get_mut( &hash ) {
            Some( (gene, _ids) ) => {
                if gene.insert( umi ){
                    match self.total_reads.get_mut( &(hash as usize) ){
                        Some( count ) => *count += 1,
                        None => panic!("This must not happen - libraray error"),
                    }
                    //eprintln!("Good data");
                    self.total_umis += 1;
                    return true
                }
                //eprintln!("UMI already known {umi}");
                false
            }, 
            None => {
                let mut gc:HashSet<u64> = HashSet::new(); //to store the umis
                gc.insert( umi );
                self.total_reads.insert( hash as usize , 1 );
                self.multimapper.insert( hash, (gc, geneids) );
                //eprintln!("Good data2");
                self.total_umis += 1;
                true
            }
        }
    }*/

    /// n_umi gets all umis stored for one cell
    pub fn n_umi( &self ) -> usize {
        self.total_umis

        // let mut n = 0;

        // for name in gnames{
        //     n += self.n_umi_4_gene( gene_info, name );
        // }
        // //println!("I got {} umis for cell {}", n, self.name );
        // n
    }

    /// This is used to calculate the subset specific umi counts after the
    /// crap cells have been discarded - definitely necessary to get it gene specific!
    pub fn n_reads( &self, gene_info:&IndexedGenes, gnames: &Vec<String> ) -> usize {
        //println!("I got {} umis for cell {}", n, self.name );
        //return self.total_umis;
        let mut n = 0;
        let gene_ids = gene_info.ids_for_gene_names( gnames );
        for id in gene_ids{
            n += match self.total_reads.get( &id  ){
                Some( count ) => {
                    *count
                },
                None => 0
            };
        }
        n
    }
    

    pub fn n_umi_4_gene_id( &self, gene_id:&usize ) -> usize{
        *self.total_reads.get( gene_id ).unwrap_or(&0)
    }
    
    pub fn to_str(&self, gene_info:&IndexedGenes, names: &Vec<String> ) -> String {

        let mut data = Vec::<std::string::String>::with_capacity( names.len()+4 ); 
        data.push(format!( "Cell{}", self.name ) );

        // here our internal data already should be stored with the same ids as the gene names.
        let mut total = 0;
        let mut max = 0;
        let mut max_name:std::string::String = "na".to_string();

        let gene_ids = gene_info.ids_for_gene_names( names );
        let mut dist2max = usize::MAX;
        for (i, id) in gene_ids.iter().enumerate(){
            let n = *self.total_reads.get( id  ).unwrap_or(&0);
            //println!("I collected expression for gene {}: n={}", name, n);
            if max < n {
                max_name = names[i].to_string();
                max = n;
            }
            data.push( n.to_string() );
            total += n;
        }
        for id in gene_ids{
            let n = *self.total_reads.get( &id  ).unwrap_or(&0);
            if (dist2max > max - n) & (max - n != 0) {
                dist2max = max - n;
                //println!("setting dist2max to new max: {dist2max} ( {max} - {n}");
            }
        }

        data.push( max_name ); // max expressing gene (or sample id in an HTO analysis)
        data.push( (max as f32 / total as f32 ).to_string()); // fraction of reads for the max gene
        data.push( ( total ).to_string());
        data.push( (dist2max as f32 / max as f32 ).to_string()); // percent max reads distance to nr.2
        //println!("This should print one LINE here: {data:?}",);
        data.join( "\t" )
    }
}
