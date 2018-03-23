#![feature(nll)]

extern crate failure;
#[macro_use]
extern crate failure_derive;

mod conf;
mod error;
mod gromos87;
mod rvec;

pub use conf::{get_or_insert_atom_and_residue, Atom, Conf, Residue, ResidueIter};
pub use rvec::RVec;
