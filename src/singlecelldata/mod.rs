// /singlecelldata/mod.rs


pub mod indexed_genes;
pub mod singlecelldata;
pub mod ambient_rna_detect;
pub mod cell_data;


pub use singlecelldata::SingleCellData as SingleCellData;
pub use crate::singlecelldata::cell_data::CellData  as CellData;
pub use crate::singlecelldata::ambient_rna_detect::AmbientRnaDetect  as AmbientRnaDetect;
pub use crate::singlecelldata::indexed_genes::IndexedGenes as IndexedGenes;