use gromos87;

use std::io;

#[derive(Debug, Fail)]
pub enum ReadError {
    #[fail(display = "Error reading GROMOS87 file ({})", _0)]
    Gromos87(gromos87::IoError),
    #[fail(display = "Error: Could not open file for reading ({})", _0)]
    IoError(io::Error),
}

impl From<io::Error> for ReadError {
    fn from(err: io::Error) -> ReadError {
        ReadError::IoError(err)
    }
}
