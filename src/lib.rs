#![feature(nll)]

extern crate failure;
#[macro_use] extern crate failure_derive;

mod conf;
mod error;
mod gromos87;
mod rvec;

pub use conf::Conf;
