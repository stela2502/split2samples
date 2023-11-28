// fast_mapper/mapper_entries/mod.rs

pub mod name_entries;
pub mod mapper_entries;

pub use name_entries::NameEntry as NameEntry;
pub use mapper_entries::MapperEntry as MapperEntry;

// Import the submodules
// mod name_entries;
// mod mapper_entries;