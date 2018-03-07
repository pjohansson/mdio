#![feature(nll)]

extern crate failure;
#[macro_use] extern crate failure_derive;

mod conf;
mod error;
mod gromos87;
mod rvec;

pub use conf::{Conf, Atom, Residue, ResidueIter, get_or_insert_atom_and_residue};
pub use rvec::RVec;
