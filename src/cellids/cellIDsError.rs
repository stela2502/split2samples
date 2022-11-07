use std::fmt;
use std::error::Error;
use core::str::Utf8Error;


type Result<'a,T> = std::result::Result<T, NnuclError<'a>>;


#[derive(Debug, Clone)]
pub struct NnuclError<'a> {
    seq: &'a str
}

impl NnuclError<'_>{
    pub fn new(seq: &str )-> NnuclError<'_>{
        Self{
            seq
        }
    }
    // pub fn From( err: Utf8Error){
    //     return NnuclError::new( "some updtream UTF8 errors - what kind of sequence have you given me?!");
    // }
}

impl fmt::Display for NnuclError<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.seq )
    }
}

impl Error for NnuclError<'_> {
    fn description(&self) -> &str {
        return  self.seq
    }
}


impl From<Utf8Error> for NnuclError<'_> {
    fn from(T: Utf8Error) -> Self{
        return NnuclError::new( "some upstream UTF8 errors - what kind of sequence have you given me?!");
    }
}