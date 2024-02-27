
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

#[derive(Debug)]
pub enum CellIdError {
    NoMatch,
    Ns,
}

#[derive(Debug)]
pub enum SeqError {
    Ns,
    End,
    LowComplexity,
}


