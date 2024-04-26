#[derive(Debug)]
pub enum GeneSelectionError {
    NotSame,
    TooView,
}

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
    TooShort,
    NoMatch,
    Ns,
}

#[derive(Debug)]
pub enum SeqError {
    Ns,
    End,
    LowComplexity,
}



