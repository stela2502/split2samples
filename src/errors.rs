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

#[derive(Debug,Clone, Copy, PartialEq)]
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
