#![allow(
    // clippy is broken and shows wrong warnings
    // clippy on stable does not know yet about the lint name
    unknown_lints,
    // https://github.com/rust-lang/rust-clippy/issues/8560
    clippy::only_used_in_recursion,
    // https://github.com/rust-lang/rust-clippy/issues/8867
    clippy::derive_partial_eq_without_eq,
    // https://github.com/rust-lang/rust-clippy/issues/9101
    clippy::explicit_auto_deref
)]


pub mod analysis;
pub mod errors;
pub mod opts_te;
pub mod cellids10x;
pub mod cellids;
pub mod geneids;
pub mod sampleids;
pub mod singlecelldata;
pub mod ofiles;
pub mod ifiles;
pub mod last5;

pub mod fast_mapper;
pub mod genes_mapper;
//pub mod fast_mapper::mapper_entries;
pub mod gene;

pub mod int_to_str;
pub mod mapping_info;
pub mod gene_family;


pub mod traits;