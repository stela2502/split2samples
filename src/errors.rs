
#[derive(Debug)]
pub enum FilterError {
    Length,
    //PolyA,
    Ns,
    Quality,
}



#[derive(Debug)]
pub enum MappingError {
    NoMatch,
    MultiMatch,
}
