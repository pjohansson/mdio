use error::{ReadError, WriteError};
use gromos87;
use rvec::RVec;

use std::cell::RefCell;
use std::fs::File;
use std::io::{BufReader, BufWriter};
// use std::ops::Deref;
use std::path::Path;
use std::rc::Rc;

/// A system configuration.
#[derive(Clone, Debug)]
pub struct Conf {
    /// Configuration title.
    pub title: String,
    /// Origin of configuration.
    pub origin: RVec,
    /// Size of configuration.
    pub size: RVec,
    /// A list of residues which exist in the configuration.
    ///
    /// These are shared, mutable references to the objects, since we might want
    /// to modify the residues after their creation.
    pub residues: Vec<Rc<RefCell<Residue>>>,
    /// A list of the atoms of the configuration.
    pub atoms: Vec<Atom>,
}

impl Conf {
    /// Read a configuration from a `Gromos87` formatted file.
    pub fn from_gromos87(path: &Path) -> Result<Conf, ReadError> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        gromos87::read_gromos87_conf(&mut reader).map_err(|err| ReadError::Gromos87(err))
    }

    /// Group atoms as their residues and iterate over them.
    pub fn iter_residues(&self) -> ResidueIter {
        ResidueIter {
            index: 0,
            atoms: &self.atoms,
        }
    }

    /// Extend the configuration along each direction by copying and translating the atoms.
    pub fn pbc_multiply(&self, nx: usize, ny: usize, nz: usize) -> Conf {
        let mut conf = Conf {
            title: self.title.clone(),
            origin: self.origin.clone(),
            size: self.size.pbc_multiply(nx, ny, nz),
            residues: self.residues.clone(),
            atoms: Vec::new(),
        };

        for ix in 1..(nx + 1) {
            for iy in 1..(ny + 1) {
                for iz in 1..(nz + 1) {
                    let dr = self.size.pbc_multiply(ix - 1, iy - 1, iz - 1);

                    self.atoms.iter().for_each(|atom| {
                        conf.atoms.push(Atom {
                            name: Rc::clone(&atom.name),
                            residue: Rc::clone(&atom.residue),
                            position: atom.position + dr,
                            velocity: atom.velocity.clone(),
                        });
                    });
                }
            }
        }

        conf
    }

    /// Write the configuration to a GROMOS87 formatted file.
    pub fn write_gromos87(&self, path: &Path) -> Result<(), WriteError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        gromos87::write_gromos87_conf(self, &mut writer).map_err(|err| WriteError::Gromos87(err))
    }
}

/// Error from iterating over residues.
#[derive(Debug, Fail)]
#[fail(display = "Bad residue starting at index {}", index)]
pub struct ResidueError {
    index: usize,
}

/// An iterator over residues of a collection of `Atom`s.
#[derive(Debug)]
pub struct ResidueIter<'a> {
    index: usize,
    atoms: &'a [Atom],
}

impl<'a> ResidueIter<'a> {
    fn get_iter_error(&mut self, i: usize) -> ResidueError {
        self.index += i;
        ResidueError { index: self.index - i }
    }
}

impl<'a> Iterator for ResidueIter<'a> {
    type Item = Result<Vec<Atom>, ResidueError>;

    fn next(&mut self) -> Option<Self::Item> {
        let atom1 = self.atoms.get(self.index)?.clone();

        let residue = atom1.residue.clone();
        let residue_len = residue.borrow().atoms.len();

        // If the first atom is wrong, return an error and skip it
        if !Rc::ptr_eq(&atom1.name, &residue.borrow().atoms[0]) {
            return Some(Err(self.get_iter_error(1)));
        }

        let mut atom_list = Vec::new();
        atom_list.push(atom1);

        for i in 1..residue_len {
            match self.atoms.get(i + self.index) {
                Some(atom) => {
                    if !Rc::ptr_eq(&atom.name, &residue.borrow().atoms[i]) {
                        return Some(Err(self.get_iter_error(i)));
                    }

                    atom_list.push(atom.clone());
                },
                None => {
                    return Some(Err(self.get_iter_error(i)));
                },
            }
        }

        self.index += residue_len;

        Some(Ok(atom_list))
    }
}

/// Information about a residue.
#[derive(Clone, Debug)]
pub struct Residue {
    /// The residue name.
    pub name: Rc<RefCell<String>>,
    /// Atoms which belong to the residue.
    pub atoms: Vec<Rc<RefCell<String>>>,
}

impl Residue {
    fn get_or_insert_atom(&mut self, atom_name: &str) -> Rc<RefCell<String>> {
        self.atoms
            .iter()
            .find(|name| &*name.borrow() == &atom_name)
            .cloned()
            .unwrap_or_else(|| {
                let atom = Rc::new(RefCell::new(String::from(atom_name)));
                self.atoms.push(atom.clone());

                atom
            })
    }
}

/// A single atom belonging to a residue in the configuration.
#[derive(Clone, Debug)]
pub struct Atom {
    /// A reference to the atom name. Should point to an atom in the `residue`.
    pub name: Rc<RefCell<String>>,
    /// A reference to the residue which owns the atom. Will typicall point to a residue
    /// in the `Conf` in which this atom exists.
    pub residue: Rc<RefCell<Residue>>,
    /// The atom position in configuration-relative coordinates.
    pub position: RVec,
    /// The atom velocity, if it has one.
    pub velocity: Option<RVec>,
}

fn get_or_insert_residue(name: &str, residues: &mut Vec<Rc<RefCell<Residue>>>)
        -> Rc<RefCell<Residue>> {
    residues.iter()
            .find(|res| *res.borrow().name.borrow() == name)
            .cloned()
            .unwrap_or_else(|| {
                let res = Rc::new(RefCell::new(Residue {
                    name: Rc::new(RefCell::new(String::from(name))),
                    atoms: Vec::new(),
                }));

                residues.push(res.clone());
                res
            })
}

pub fn get_or_insert_atom_and_residue(residue_name: &str,
                                      atom_name: &str,
                                      residues: &mut Vec<Rc<RefCell<Residue>>>)
        -> Result<(Rc<RefCell<Residue>>, Rc<RefCell<String>>), String> {
    let residue = get_or_insert_residue(residue_name, residues);
    let atom = residue.borrow_mut().get_or_insert_atom(atom_name);

    Ok((residue, atom))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env::temp_dir;

    #[test]
    fn get_or_insert_residue_from_list() {
        let mut residues = Vec::new();

        let res1_name = "RES1";
        let res1 = get_or_insert_residue(res1_name, &mut residues);

        assert_eq!(*res1.borrow().name.borrow(), res1_name);
        assert!(&res1.borrow().atoms.is_empty());

        assert_eq!(residues.len(), 1);
        assert!(Rc::ptr_eq(&res1, &residues[0]));

        let res1_again = get_or_insert_residue(res1_name, &mut residues);
        assert!(Rc::ptr_eq(&res1, &res1_again));

        let res2_name = "RES2";
        let res2 = get_or_insert_residue(res2_name, &mut residues);

        assert_eq!(*res2.borrow().name.borrow(), res2_name);
        assert!(&res2.borrow().atoms.is_empty());
        assert!(!Rc::ptr_eq(&res1, &res2));

        assert_eq!(residues.len(), 2);
        assert!(Rc::ptr_eq(&res2, &residues[1]));
    }

    #[test]
    fn get_or_insert_atom_from_residue() {
        let mut residue = Residue {
            name: Rc::new(RefCell::new(String::from("RES"))),
            atoms: Vec::new(),
        };

        let atom1_name = "ATOM1";
        let atom1 = residue.get_or_insert_atom(atom1_name);

        assert_eq!(&*atom1.borrow(), atom1_name);
        assert!(Rc::ptr_eq(&atom1, &residue.atoms[0]));

        let atom1_again = residue.get_or_insert_atom(atom1_name);
        assert!(Rc::ptr_eq(&atom1_again, &atom1));

        let atom2_name = "ATOM2";
        let atom2 = residue.get_or_insert_atom(atom2_name);

        assert_eq!(&*atom2.borrow(), atom2_name);
        assert!(Rc::ptr_eq(&atom2, &residue.atoms[1]));
        assert!(!Rc::ptr_eq(&atom1, &atom2));
    }

    #[test]
    fn get_atom_and_residue_from_list() {
        let mut residues = Vec::new();

        let res1_name = "RES1";
        let atom1_name = "AT1";

        let (res1, atom1) = get_or_insert_atom_and_residue(
            res1_name, atom1_name, &mut residues).unwrap();

        assert_eq!(*res1.borrow().name.borrow(), res1_name);
        assert_eq!(&*atom1.borrow(), &atom1_name);
        assert!(Rc::ptr_eq(&atom1, &res1.borrow().atoms[0]));

        let atom2_name = "AT2";
        let (res1_again, atom2) = get_or_insert_atom_and_residue(
            res1_name, atom2_name, &mut residues).unwrap();

        assert!(Rc::ptr_eq(&res1, &res1_again));
        assert_eq!(&*atom2.borrow(), &atom2_name);

        let res2_name = "RES2";
        let atom3_name = "AT3";

        let (res2, atom3) = get_or_insert_atom_and_residue(
            res2_name, atom3_name, &mut residues).unwrap();

        assert!(!Rc::ptr_eq(&res1, &res2));
        assert_eq!(*res2.borrow().name.borrow(), res2_name);
        assert_eq!(&*atom3.borrow(), &atom3_name);

        // An atom with a name of another residue can be added, they will not be the same
        let (res2_again, atom1_not_res1) = get_or_insert_atom_and_residue(
            res2_name, atom1_name, &mut residues).unwrap();

        assert!(Rc::ptr_eq(&res2, &res2_again));
        assert!(!Rc::ptr_eq(&atom1, &atom1_not_res1));
    }

    #[test]
    fn read_bad_filename_gives_error() {
        let mut filename = temp_dir();
        filename.push("_file_should_not_exist_mdio_test_");

        assert!(Conf::from_gromos87(&filename).is_err());
    }

    #[test]
    fn residue_iter_on_empty_conf_returns_none() {
        let conf = Conf {
            title: "A title".to_string(),
            origin: RVec { x: 0.0, y: 0.0, z: 0.0 },
            size: RVec { x: 0.0, y: 0.0, z: 0.0 },
            residues: Vec::new(),
            atoms: Vec::new(),
        };

        let mut iter = conf.iter_residues();

        assert!(iter.next().is_none());
    }

    #[test]
    fn residue_iter_over_two_atoms_of_different_residues() {
        let residues = vec![
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES1".to_string())),
                atoms: vec![Rc::new(RefCell::new("AT1".to_string()))],
            })),
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES2".to_string())),
                atoms: vec![Rc::new(RefCell::new("AT2".to_string()))],
            })),
        ];

        let conf = Conf {
            title: "A title".to_string(),
            origin: RVec { x: 0.0, y: 0.0, z: 0.0 },
            size: RVec { x: 0.0, y: 0.0, z: 0.0 },
            residues: residues.clone(),
            atoms: vec![
                // Residue 2
                Atom {
                    name: Rc::clone(&residues[1].borrow().atoms[0]),
                    residue: Rc::clone(&residues[1]),
                    position: RVec { x: 0.0, y: 1.0, z: 2.0 },
                    velocity: Some(RVec { x: 0.0, y: 0.1, z: 0.2 }),
                },
                // Residue 1
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 3.0, y: 4.0, z: 5.0 },
                    velocity: Some(RVec { x: 0.3, y: 0.4, z: 0.5 }),
                },
            ]
        };

        let mut iter = conf.iter_residues();

        let res = iter.next().unwrap().unwrap();
        assert_eq!(res.len(), 1);
        assert!(Rc::ptr_eq(&res[0].residue, &residues[1]));
        assert!(Rc::ptr_eq(&res[0].name, &residues[1].borrow().atoms[0]));
        assert_eq!(res[0].position, RVec { x: 0.0, y: 1.0, z: 2.0 });
        assert_eq!(res[0].velocity.unwrap(), RVec { x: 0.0, y: 0.1, z: 0.2 });

        let res = iter.next().unwrap().unwrap();
        assert_eq!(res.len(), 1);
        assert!(Rc::ptr_eq(&res[0].residue, &residues[0]));
        assert!(Rc::ptr_eq(&res[0].name, &residues[0].borrow().atoms[0]));
        assert_eq!(res[0].position, RVec { x: 3.0, y: 4.0, z: 5.0 });
        assert_eq!(res[0].velocity.unwrap(), RVec { x: 0.3, y: 0.4, z: 0.5 });

        assert!(iter.next().is_none());
    }

    #[test]
    fn iterate_over_a_residue_with_several_atoms() {
        let residues = vec![
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES1".to_string())),
                atoms: vec![
                    Rc::new(RefCell::new("AT1".to_string())),
                    Rc::new(RefCell::new("AT2".to_string()))
                ],
            })),
        ];

        let conf = Conf {
            title: "A title".to_string(),
            origin: RVec { x: 0.0, y: 0.0, z: 0.0 },
            size: RVec { x: 0.0, y: 0.0, z: 0.0 },
            residues: residues.clone(),
            atoms: vec![
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 0.0, y: 1.0, z: 2.0 },
                    velocity: None,
                },
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[1]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 3.0, y: 4.0, z: 5.0 },
                    velocity: None,
                },
            ]
        };

        let mut iter = conf.iter_residues();

        let res = iter.next().unwrap().unwrap();
        assert_eq!(res.len(), 2);

        assert!(Rc::ptr_eq(&res[0].residue, &residues[0]));
        assert!(Rc::ptr_eq(&res[0].name, &residues[0].borrow().atoms[0]));
        assert_eq!(res[0].position, RVec { x: 0.0, y: 1.0, z: 2.0 });

        assert!(Rc::ptr_eq(&res[1].residue, &residues[0]));
        assert!(Rc::ptr_eq(&res[1].name, &residues[0].borrow().atoms[1]));
        assert_eq!(res[1].position, RVec { x: 3.0, y: 4.0, z: 5.0 });

        assert!(iter.next().is_none());
    }

    #[test]
    fn iterating_over_residues_ensures_that_all_are_consistent() {
        let residues = vec![
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES1".to_string())),
                atoms: vec![
                    Rc::new(RefCell::new("AT1".to_string())),
                    Rc::new(RefCell::new("AT2".to_string()))
                ],
            })),
        ];

        let conf = Conf {
            title: "A title".to_string(),
            origin: RVec { x: 0.0, y: 0.0, z: 0.0 },
            size: RVec { x: 0.0, y: 0.0, z: 0.0 },
            residues: residues.clone(),
            atoms: vec![
                // Complete residue
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 0.0, y: 1.0, z: 2.0 },
                    velocity: None,
                },
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[1]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 3.0, y: 4.0, z: 5.0 },
                    velocity: None,
                },
                // Incomplete residue: misses second atom
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 0.0, y: 1.0, z: 2.0 },
                    velocity: None,
                },
                // A final complete residue
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 6.0, y: 7.0, z: 8.0 },
                    velocity: None,
                },
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[1]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 9.0, y: 10.0, z: 11.0 },
                    velocity: None,
                },
            ]
        };

        let mut iter = conf.iter_residues();

        assert!(iter.next().unwrap().is_ok());

        // Second gives error
        assert!(iter.next().unwrap().is_err());

        // Third recovers (TODO: Decide whether this should be the case)
        let res = iter.next().unwrap().unwrap();
        assert_eq!(res.len(), 2);

        assert!(Rc::ptr_eq(&res[0].residue, &residues[0]));
        assert!(Rc::ptr_eq(&res[0].name, &residues[0].borrow().atoms[0]));
        assert_eq!(res[0].position, RVec { x: 6.0, y: 7.0, z: 8.0 });

        assert!(Rc::ptr_eq(&res[1].residue, &residues[0]));
        assert!(Rc::ptr_eq(&res[1].name, &residues[0].borrow().atoms[1]));
        assert_eq!(res[1].position, RVec { x: 9.0, y: 10.0, z: 11.0 });

        assert!(iter.next().is_none());
    }

    #[test]
    fn iterating_over_residues_ensures_that_they_are_ordered() {
        let residues = vec![
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES1".to_string())),
                atoms: vec![
                    Rc::new(RefCell::new("AT1".to_string())),
                    Rc::new(RefCell::new("AT2".to_string()))
                ],
            })),
        ];

        let conf = Conf {
            title: "A title".to_string(),
            origin: RVec { x: 0.0, y: 0.0, z: 0.0 },
            size: RVec { x: 0.0, y: 0.0, z: 0.0 },
            residues: residues.clone(),
            atoms: vec![
                // Residue begins with wrong atom, and skipped
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[1]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 0.0, y: 1.0, z: 2.0 },
                    velocity: None,
                },
                // This residue (which along with the previous atom is a good residue)
                // is found as incomplete and skipped
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 0.0, y: 1.0, z: 2.0 },
                    velocity: None,
                },
                // The next residue is good
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 6.0, y: 7.0, z: 8.0 },
                    velocity: None,
                },
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[1]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 9.0, y: 10.0, z: 11.0 },
                    velocity: None,
                },
            ]
        };

        let mut iter = conf.iter_residues();

        // First and second residues will be bad (both are incomplete)
        assert!(iter.next().unwrap().is_err());
        assert!(iter.next().unwrap().is_err());

        // This is good
        let res = iter.next().unwrap().unwrap();
        assert_eq!(res.len(), 2);

        assert!(Rc::ptr_eq(&res[0].residue, &residues[0]));
        assert!(Rc::ptr_eq(&res[0].name, &residues[0].borrow().atoms[0]));
        assert_eq!(res[0].position, RVec { x: 6.0, y: 7.0, z: 8.0 });

        assert!(Rc::ptr_eq(&res[1].residue, &residues[0]));
        assert!(Rc::ptr_eq(&res[1].name, &residues[0].borrow().atoms[1]));
        assert_eq!(res[1].position, RVec { x: 9.0, y: 10.0, z: 11.0 });

        assert!(iter.next().is_none());
    }

    #[test]
    fn iterate_over_several_different_residues() {
        let residues = vec![
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES1".to_string())),
                atoms: vec![
                    Rc::new(RefCell::new("AT1".to_string())),
                    Rc::new(RefCell::new("At2".to_string())),
                ],
            })),
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES2".to_string())),
                atoms: vec![
                    Rc::new(RefCell::new("AT3".to_string())),
                ],
            })),
        ];

        // This configuration contains 2 of the first residue, then 2 of the second,
        // and finally 1 of the first
        let atoms = vec![
            Atom {
                name: residues[0].borrow().atoms[0].clone(),
                residue: residues[0].clone(),
                position: RVec { x: 0.0, y: 1.0, z: 2.0 },
                velocity: None,
            },
            Atom {
                name: residues[0].borrow().atoms[1].clone(),
                residue: residues[0].clone(),
                position: RVec { x: 3.0, y: 4.0, z: 5.0 },
                velocity: None,
            },
            Atom {
                name: residues[0].borrow().atoms[0].clone(),
                residue: residues[0].clone(),
                position: RVec { x: 6.0, y: 7.0, z: 8.0 },
                velocity: None,
            },
            Atom {
                name: residues[0].borrow().atoms[1].clone(),
                residue: residues[0].clone(),
                position: RVec { x: 9.0, y: 10.0, z: 11.0 },
                velocity: None,
            },
            Atom {
                name: residues[1].borrow().atoms[0].clone(),
                residue: residues[1].clone(),
                position: RVec { x: 12.0, y: 13.0, z: 14.0 },
                velocity: None,
            },
            Atom {
                name: residues[1].borrow().atoms[0].clone(),
                residue: residues[1].clone(),
                position: RVec { x: 15.0, y: 16.0, z: 17.0 },
                velocity: None,
            },
            Atom {
                name: residues[0].borrow().atoms[0].clone(),
                residue: residues[0].clone(),
                position: RVec { x: 18.0, y: 19.0, z: 20.0 },
                velocity: None,
            },
            Atom {
                name: residues[0].borrow().atoms[1].clone(),
                residue: residues[0].clone(),
                position: RVec { x: 21.0, y: 22.0, z: 23.0 },
                velocity: None,
            },
        ];

        let conf = Conf {
            title: "System".to_string(),
            origin: RVec { x: 0.0, y: 0.0, z: 0.0 },
            size: RVec { x: 1.0, y: 2.0, z: 3.0 },
            residues: residues.clone(),
            atoms,
        };

        let mut iter = conf.iter_residues();

        // Check the fourth and fifth (final) residues
        assert!(iter.next().unwrap().is_ok());
        assert!(iter.next().unwrap().is_ok());
        assert!(iter.next().unwrap().is_ok());

        let res4 = iter.next().unwrap().unwrap();
        assert_eq!(res4.len(), 1);
        assert!(Rc::ptr_eq(&res4[0].residue, &residues[1]));
        assert!(Rc::ptr_eq(&res4[0].name, &residues[1].borrow().atoms[0]));
        assert_eq!(res4[0].position, RVec { x: 15.0, y: 16.0, z: 17.0 });
        assert_eq!(res4[0].velocity, None);

        let res5 = iter.next().unwrap().unwrap();
        assert_eq!(res5.len(), 2);

        assert!(Rc::ptr_eq(&res5[0].residue, &residues[0]));
        assert!(Rc::ptr_eq(&res5[0].name, &residues[0].borrow().atoms[0]));
        assert_eq!(res5[0].position, RVec { x: 18.0, y: 19.0, z: 20.0 });
        assert_eq!(res5[0].velocity, None);

        assert!(Rc::ptr_eq(&res5[1].residue, &residues[0]));
        assert!(Rc::ptr_eq(&res5[1].name, &residues[0].borrow().atoms[1]));
        assert_eq!(res5[1].position, RVec { x: 21.0, y: 22.0, z: 23.0 });
        assert_eq!(res5[1].velocity, None);

        assert!(iter.next().is_none());
    }

    #[test]
    fn multiply_conf_to_extend_it() {
        let size = RVec { x: 10.0, y: 20.0, z: 30.0 };

        let residues = vec![
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES1".to_string())),
                atoms: vec![Rc::new(RefCell::new("AT1".to_string()))],
            })),
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES2".to_string())),
                atoms: vec![Rc::new(RefCell::new("AT2".to_string()))],
            })),
        ];

        let conf = Conf {
            title: "A title".to_string(),
            origin: RVec { x: 0.0, y: 0.0, z: 0.0 },
            size,
            residues: residues.clone(),
            atoms: vec![
                Atom {
                    name: Rc::clone(&residues[1].borrow().atoms[0]),
                    residue: Rc::clone(&residues[1]),
                    position: RVec { x: 0.0, y: 1.0, z: 2.0 },
                    velocity: Some(RVec { x: 0.0, y: 0.1, z: 0.2 }),
                },
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec { x: 3.0, y: 4.0, z: 5.0 },
                    velocity: Some(RVec { x: 0.3, y: 0.4, z: 0.5 }),
                },
            ]
        };

        let (nx, ny, nz) = (2, 3, 4);
        let multiplied_conf = conf.pbc_multiply(nx, ny, nz);

        assert_eq!(multiplied_conf.size,
            RVec { x: 10.0 * (nx as f64), y: 20.0 * (ny as f64), z: 30.0 * (nz as f64) }
        );
        assert_eq!(multiplied_conf.atoms.len(), conf.atoms.len() * nx * ny * nz);

        // The final atom should be from the maximum (nx, ny, nz) image
        assert!(Rc::ptr_eq(
            &multiplied_conf.atoms.last().unwrap().name, &conf.atoms.last().unwrap().name
        ));
        assert!(Rc::ptr_eq(
            &multiplied_conf.atoms.last().unwrap().residue, &conf.atoms.last().unwrap().residue
        ));
        assert_eq!(
            multiplied_conf.atoms.last().unwrap().position,
            conf.atoms.last().unwrap().position + conf.size.pbc_multiply(nx - 1, ny - 1, nz - 1)
        );
        assert_eq!(
            multiplied_conf.atoms.last().unwrap().velocity,
            conf.atoms.last().unwrap().velocity
        );
    }
}
