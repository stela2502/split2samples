// fast_mapper/mapper_entries/mod.rs

pub mod name_entries;
pub mod mapper_entries;
pub mod second_seq;

pub use name_entries::NameEntry as NameEntry;
pub use mapper_entries::MapperEntry as MapperEntry;
pub use second_seq::SecondSeq as SecondSeq;

// Import the submodules
// mod name_entries;
// mod mapper_entries;