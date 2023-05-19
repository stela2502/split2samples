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

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}


pub mod cellids;
pub mod geneids;
pub mod sampleids;
pub mod singlecelldata;
pub mod ofiles;
pub mod ifiles;
pub mod analysis;
pub mod last5;
pub mod fast_mapper;
pub mod mapper_entries;
pub mod gene;
pub mod int_to_str;
pub mod mapping_info;