use gromos87;

use std::io;

#[derive(Debug, Fail)]
pub enum WriteError {
    #[fail(display = "Could not write GROMOS87 file ({})", _0)]
    Gromos87(gromos87::WriteError),
    #[fail(display = "Could not open file for writing ({})", _0)]
    IoError(io::Error),
}

impl From<io::Error> for WriteError {
    fn from(err: io::Error) -> WriteError {
        WriteError::IoError(err)
    }
}

#[derive(Debug, Fail)]
pub enum ReadError {
    #[fail(display = "Could not read GROMOS87 file ({})", _0)]
    Gromos87(gromos87::ReadError),
    #[fail(display = "Could not open file for reading ({})", _0)]
    IoError(io::Error),
}

impl From<io::Error> for ReadError {
    fn from(err: io::Error) -> ReadError {
        ReadError::IoError(err)
    }
}
