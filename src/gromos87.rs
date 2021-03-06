use conf::{get_or_insert_atom_and_residue, Atom, Conf};
use rvec::{ParseRVecError, RVec};

use std::io;
use std::io::{BufRead, BufReader, Read, Write};

pub fn write_gromos87_conf<W: Write>(conf: &Conf, mut writer: &mut W) -> Result<(), WriteError> {
    write!(&mut writer, "{}\n{}\n", conf.title, conf.atoms.len())?;

    let mut atom_num = 0;

    for (res_num, residue) in conf.iter_residues().enumerate() {
        // GROMOS-87 wraps indices at 5 digits, ie. at 100_000
        let res_num_wrapped = (res_num + 1) % 100_000;

        for atom in residue
            .map_err(|_| WriteError::BadResidue(res_num + 1))?
            .iter()
        {
            atom_num += 1;
            let atom_num_wrapped = atom_num % 100_000;

            write!(
                &mut writer,
                "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                res_num_wrapped,
                atom.residue.borrow().name.borrow(),
                *atom.name.borrow(),
                atom_num_wrapped,
                atom.position.x,
                atom.position.y,
                atom.position.z
            )?;

            if let Some(velocity) = atom.velocity {
                write!(
                    &mut writer,
                    "{:>8.4}{:>8.4}{:>8.4}",
                    velocity.x, velocity.y, velocity.z
                )?;
            }

            write!(&mut writer, "\n")?;
        }
    }

    write!(
        &mut writer,
        " {:12.5} {:12.5} {:12.5}\n",
        conf.size.x, conf.size.y, conf.size.z
    )?;

    Ok(())
}

struct Line<'a> {
    // residue_number: usize,
    residue_name: &'a str,
    atom_name: &'a str,
    // atom_number: usize,
    position: RVec,
    velocity: Option<RVec>,
}

#[derive(Debug, Fail)]
pub enum WriteError {
    #[fail(display = "Error writing configuration ({})", _0)]
    IoError(io::Error),
    #[fail(display = "Error writing residue {}, which was incomplete", _0)]
    BadResidue(usize),
}

impl From<io::Error> for WriteError {
    fn from(err: io::Error) -> WriteError {
        WriteError::IoError(err)
    }
}

#[derive(Debug, Fail)]
pub enum ReadError {
    #[fail(display = "Could not read line {}: invalid UTF-8", _0)]
    Utf8Error(usize),
    #[fail(display = "Expected a configuration title at line 1")]
    MissingTitle,
    #[fail(display = "Expected a number of atoms entry at line 2")]
    MissingNumAtoms,
    #[fail(display = "Could not parse number of atoms entry at line 2")]
    NumAtomsError,
    #[fail(display = "Expected an atom entry at line {}", _0)]
    MissingAtomLine(usize),
    #[fail(display = "Could not parse atom entry at line {}", _0)]
    LineError(usize),
    #[fail(display = "Expected a configuration box size entry at line {}", _0)]
    NoBoxSize(usize),
    #[fail(display = "Could not parse box size entry at line {}", _0)]
    BoxSizeError(usize),
}

pub fn read_gromos87_conf<R: Read>(reader: R) -> Result<Conf, ReadError> {
    let mut buf_reader = BufReader::new(reader);
    let mut buf = String::new();

    buf_reader
        .read_line(&mut buf)
        .map_err(|_| ReadError::Utf8Error(1))?;
    let title = buf.trim().to_string();
    buf.clear();

    buf_reader
        .read_line(&mut buf)
        .map_err(|_| ReadError::Utf8Error(1))?;
    let num_atoms = buf.trim()
        .parse::<usize>()
        .map_err(|_| ReadError::NumAtomsError)?;
    buf.clear();

    let mut residues = Vec::new();
    let mut atoms = Vec::new();

    for i in 0..num_atoms {
        buf_reader
            .read_line(&mut buf)
            .map_err(|_| ReadError::Utf8Error(2 + i))?;

        let atom_line = parse_atom_line(&buf).map_err(|_| ReadError::LineError(2 + i))?;
        let (residue, atom) = get_or_insert_atom_and_residue(
            atom_line.residue_name,
            atom_line.atom_name,
            &mut residues,
        ).map_err(|_| ReadError::LineError(2 + i))?;

        atoms.push(Atom {
            name: atom,
            residue,
            position: atom_line.position,
            velocity: atom_line.velocity,
        });

        buf.clear();
    }

    buf_reader
        .read_line(&mut buf)
        .map_err(|_| ReadError::Utf8Error(3 + num_atoms))?;
    let size = RVec::from_whitespace(&buf).expect("could not read box size");

    Ok(Conf {
        title,
        origin: RVec {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        size,
        residues,
        atoms,
    })
}

#[derive(Debug, Fail)]
#[fail(display = "Could not parse a line")]
struct ParseLineError;

fn parse_atom_line(line: &str) -> Result<Line, ParseLineError> {
    const GRO_MINLINELEN: usize = 44;
    if line.len() < GRO_MINLINELEN {
        return Err(ParseLineError);
    }

    // let residue_number = line[0..5].trim().parse::<usize>().map_err(|_| ParseLineError)?;
    let residue_name = line[5..10].trim();
    let atom_name = line[10..15].trim();
    // let atom_number = line[15..20].trim().parse::<usize>().map_err(|_| ParseLineError)?;

    let position = RVec::from_fixed(&line[20..], 8).map_err(|_| ParseLineError)?;
    let velocity = match RVec::from_fixed(&line[44..], 8) {
        Ok(rvec) => Some(rvec),
        Err(ParseRVecError::MissingValues) => None,
        _ => return Err(ParseLineError),
    };

    Ok(Line {
        // residue_number,
        residue_name,
        atom_name,
        // atom_number,
        position,
        velocity,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use conf::{Atom, Conf, Residue};
    use std::cell::RefCell;
    use std::io::Cursor;
    use std::rc::Rc;

    #[test]
    fn parse_atom_line_errors() {
        // Too-short strings
        assert!(parse_atom_line("").is_err());
        assert!(parse_atom_line("    1RES   ATOM1    ").is_err());
        assert!(parse_atom_line("    1RES   ATOM1    1000.0012000.002").is_err());

        // Baseline correct line
        assert!(parse_atom_line("    1RES   ATOM1    1000.0012000.0023000.003").is_ok());

        // Bad number values
        // assert!(parse_atom_line("    sRES   ATOM1    1000.0012000.0023000.003").is_err());
        // assert!(parse_atom_line("    1RES   ATOM1 s  1000.0012000.0023000.003").is_err());
        assert!(parse_atom_line("    1RES   ATOM1    100s.0012000.0023000.003").is_err());
        assert!(parse_atom_line("    1RES   ATOM1    1000.00120s0.0023000.003").is_err());
        assert!(parse_atom_line("    1RES   ATOM1    1000.0012000.00230s0.003").is_err());
    }

    #[test]
    fn parse_correct_atom_lines() {
        let s = "    1RES   ATOM1    1000.0012000.0023000.003";
        let line = parse_atom_line(s).unwrap();
        // assert_eq!(line.residue_number, 1);
        // assert_eq!(line.atom_number, 1);
        assert_eq!(line.residue_name, "RES");
        assert_eq!(line.atom_name, "ATOM");
        assert_eq!(
            line.position,
            RVec {
                x: 1000.001,
                y: 2000.002,
                z: 3000.003,
            }
        );
        assert_eq!(line.velocity, None);

        let s = "    12RES12ATO150001 100.01  200.02  300.03  400.04  500.05  600.06 ";
        let line = parse_atom_line(s).unwrap();
        // assert_eq!(line.residue_number, 1);
        // assert_eq!(line.atom_number, 50001);
        assert_eq!(line.residue_name, "2RES1");
        assert_eq!(line.atom_name, "2ATO1");
        assert_eq!(
            line.position,
            RVec {
                x: 100.01,
                y: 200.02,
                z: 300.03,
            }
        );
        assert_eq!(
            line.velocity,
            Some(RVec {
                x: 400.04,
                y: 500.05,
                z: 600.06,
            })
        );
    }

    #[test]
    fn read_correct_file() {
        let title = "A title";

        let res1_name = "RES1";
        let res2_name = "RES2";

        let atom1_name = "AT1";
        let atom1_pos = RVec {
            x: 0.0,
            y: 1.0,
            z: 2.0,
        };
        let atom1_vel = RVec {
            x: 0.0,
            y: 0.1,
            z: 0.3,
        };
        let atom2_name = "AT2";
        let atom2_pos = RVec {
            x: 3.0,
            y: 4.0,
            z: 5.0,
        };
        let atom2_vel = RVec {
            x: 0.3,
            y: 0.4,
            z: 0.5,
        };
        let atom3_name = "AT3";
        let atom3_pos1 = RVec {
            x: 6.0,
            y: 7.0,
            z: 8.0,
        };
        let atom3_vel1 = RVec {
            x: 0.6,
            y: 0.7,
            z: 0.8,
        };
        let atom3_pos2 = RVec {
            x: 9.0,
            y: 10.0,
            z: 11.0,
        };
        let atom3_vel2 = RVec {
            x: 0.9,
            y: 1.0,
            z: 1.1,
        };

        let size = RVec {
            x: 10.0,
            y: 11.0,
            z: 12.0,
        };

        let content = format!(
            "\
{}
4
{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.4}{:>8.4}{:>8.4}
{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.4}{:>8.4}{:>8.4}
{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.4}{:>8.4}{:>8.4}
{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.4}{:>8.4}{:>8.4}
{:12.3} {:12.3} {:12.3}
",
            title,
            1,
            res1_name,
            atom1_name,
            1,
            atom1_pos.x,
            atom1_pos.y,
            atom1_pos.z,
            atom1_vel.x,
            atom1_vel.y,
            atom1_vel.z,
            1,
            res1_name,
            atom2_name,
            2,
            atom2_pos.x,
            atom2_pos.y,
            atom2_pos.z,
            atom2_vel.x,
            atom2_vel.y,
            atom2_vel.z,
            2,
            res2_name,
            atom3_name,
            3,
            atom3_pos1.x,
            atom3_pos1.y,
            atom3_pos1.z,
            atom3_vel1.x,
            atom3_vel1.y,
            atom3_vel1.z,
            3,
            res2_name,
            atom3_name,
            4,
            atom3_pos2.x,
            atom3_pos2.y,
            atom3_pos2.z,
            atom3_vel2.x,
            atom3_vel2.y,
            atom3_vel2.z,
            size.x,
            size.y,
            size.z
        );

        let conf = read_gromos87_conf(content.as_bytes()).unwrap();

        assert_eq!(conf.title, title);
        assert_eq!(conf.atoms.len(), 4);
        assert_eq!(
            conf.origin,
            RVec {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            }
        );
        assert_eq!(conf.size, size);

        // Verify that all residues were correctly constructed
        assert_eq!(conf.residues.len(), 2);
        assert_eq!(*conf.residues[0].borrow().name.borrow(), res1_name);
        assert_eq!(*conf.residues[1].borrow().name.borrow(), res2_name);

        assert_eq!(conf.residues[0].borrow().atoms.len(), 2);
        assert_eq!(&*conf.residues[0].borrow().atoms[0].borrow(), atom1_name);
        assert_eq!(&*conf.residues[0].borrow().atoms[1].borrow(), atom2_name);

        assert_eq!(conf.residues[1].borrow().atoms.len(), 1);
        assert_eq!(&*conf.residues[1].borrow().atoms[0].borrow(), atom3_name);

        // Verify that pointers are correctly set from atoms to their residues
        assert!(Rc::ptr_eq(
            &conf.atoms[0].name,
            &conf.residues[0].borrow().atoms[0]
        ));
        assert!(Rc::ptr_eq(&conf.atoms[1].residue, &conf.residues[0]));
        assert!(Rc::ptr_eq(
            &conf.atoms[2].name,
            &conf.residues[1].borrow().atoms[0]
        ));
        assert!(Rc::ptr_eq(&conf.atoms[2].residue, &conf.residues[1]));

        // Verify one position and velocity
        assert_eq!(conf.atoms[2].position, atom3_pos1);
        assert_eq!(conf.atoms[2].velocity, Some(atom3_vel1));
    }

    #[test]
    fn read_incorrect_file_returns_error() {
        let size_line = "0.1 0.2 0.3";
        let two_atom_lines = "\
    1RES1   AT1    1   0.000   1.000   2.000   0.000   0.100   0.300
    1RES1   AT2    2   3.000   4.000   5.000   0.300   0.400   0.500";

        let content = format!(
            "{}\n{}\n{}\n{}",
            "More atom lines than atoms", 1, two_atom_lines, size_line
        );
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!(
            "{}\n{}\n{}\n{}",
            "Fewer atom lines than atoms", 3, two_atom_lines, size_line
        );
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!("{}\n{}\n{}", "No box size line", 2, two_atom_lines);
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!("{}\n{}\n{}1.0 2.0", "Bad box size line", 2, two_atom_lines);
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!(
            "{}\n{}\n{}1.s 2.0 3.0",
            "Bad box size line", 2, two_atom_lines
        );
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!(
            "{}\n{}\n{}",
            "No number of atoms line", two_atom_lines, size_line
        );
        assert!(read_gromos87_conf(content.as_bytes()).is_err());
    }

    #[test]
    fn write_conf_with_two_different_residues_to_buffer() {
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
            origin: RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            },
            size: RVec {
                x: 10.0,
                y: 20.0,
                z: 30.0,
            },
            residues: residues.clone(),
            atoms: vec![
                // Residue 2
                Atom {
                    name: Rc::clone(&residues[1].borrow().atoms[0]),
                    residue: Rc::clone(&residues[1]),
                    position: RVec {
                        x: 0.0,
                        y: 1.0,
                        z: 2.0,
                    },
                    velocity: Some(RVec {
                        x: 0.0,
                        y: 0.1,
                        z: 0.2,
                    }),
                },
                // Residue 1
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec {
                        x: 3.0,
                        y: 4.0,
                        z: 5.0,
                    },
                    velocity: Some(RVec {
                        x: 0.3,
                        y: 0.4,
                        z: 0.5,
                    }),
                },
            ],
        };

        // Write the configuration to a buffer
        let mut buf = Cursor::new(Vec::<u8>::new());
        assert!(write_gromos87_conf(&conf, &mut buf).is_ok());

        // Verify that we get the same configuration back after reading the output
        buf.set_position(0);
        let read_conf = read_gromos87_conf(buf).unwrap();

        assert_eq!(read_conf.title, conf.title);
        assert_eq!(
            read_conf.origin,
            RVec {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            }
        );
        assert_eq!(read_conf.size, conf.size);

        assert_eq!(read_conf.residues.len(), 2);
        assert_eq!(read_conf.atoms.len(), conf.atoms.len());

        for (read_atom, atom) in read_conf.atoms.iter().zip(conf.atoms.iter()) {
            assert_eq!(*read_atom.name.borrow(), *atom.name.borrow());
            assert_eq!(
                *read_atom.residue.borrow().name,
                *atom.residue.borrow().name
            );
            assert_eq!(read_atom.position, atom.position);
            assert_eq!(read_atom.velocity.unwrap(), atom.velocity.unwrap());
        }
    }

    #[test]
    fn box_size_is_written_in_a_fixed_format_with_leading_space_for_all_dimensions() {
        let conf = Conf {
            title: "A title".to_string(),
            origin: RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            },
            size: RVec {
                x: 10.0,
                y: 20.0,
                z: 30.0,
            },
            residues: Vec::new(),
            atoms: Vec::new(),
        };

        let mut buf = Cursor::new(Vec::<u8>::new());
        assert!(write_gromos87_conf(&conf, &mut buf).is_ok());

        buf.set_position(0);
        let box_size_line = buf.lines().skip(2).next().unwrap().unwrap();

        assert_eq!(
            format!(
                " {:12.5} {:12.5} {:12.5}",
                conf.size.x, conf.size.y, conf.size.z
            ),
            box_size_line
        );
    }

    #[test]
    fn writing_residue_and_atom_numbers_wrap_at_100_000() {
        let residues = vec![
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES1".to_string())),
                atoms: vec![Rc::new(RefCell::new("AT1".to_string()))],
            })),
        ];

        let conf = Conf {
            title: "A title".to_string(),
            origin: RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            },
            size: RVec {
                x: 10.0,
                y: 20.0,
                z: 30.0,
            },
            residues: residues.clone(),

            // Add 100_000 atoms, since indexing begins at 1 the last atom will wrap to 0!
            atoms: vec![
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec {
                        x: 0.0,
                        y: 1.0,
                        z: 2.0,
                    },
                    velocity: None,
                };
                100_000
            ],
        };

        // Write the configuration to a buffer
        let mut buf = Cursor::new(Vec::<u8>::new());
        assert!(write_gromos87_conf(&conf, &mut buf).is_ok());

        buf.set_position(0);

        // Two meta data lines, then 99_999 atom lines are skipped
        let tail = buf.lines().skip(100_001).next().unwrap().unwrap();
        let residue_num = &tail[0..5];
        let atom_num = &tail[15..20];

        assert_eq!(residue_num, "    0");
        assert_eq!(atom_num, "    0");

        // Verify that the names are untouched
        let residue_name = &tail[5..10];
        let atom_name = &tail[10..15];
        assert_eq!(residue_name, "RES1 ");
        assert_eq!(atom_name, "  AT1");
    }

    #[test]
    fn write_conf_with_3_digit_position_precision_and_four_digit_velocity_precision() {
        let residues = vec![
            Rc::new(RefCell::new(Residue {
                name: Rc::new(RefCell::new("RES1".to_string())),
                atoms: vec![Rc::new(RefCell::new("AT1".to_string()))],
            })),
        ];

        let conf = Conf {
            title: "A title".to_string(),
            origin: RVec::default(),
            size: RVec::default(),
            residues: residues.clone(),
            atoms: vec![
                Atom {
                    name: Rc::clone(&residues[0].borrow().atoms[0]),
                    residue: Rc::clone(&residues[0]),
                    position: RVec {
                        x: 0.0,
                        y: 1.0,
                        z: 2.0,
                    },
                    velocity: Some(RVec {
                        x: 0.0,
                        y: 0.1,
                        z: 0.2,
                    }),
                },
            ],
        };

        // Write the configuration to a buffer
        let mut buf = Cursor::new(Vec::<u8>::new());
        assert!(write_gromos87_conf(&conf, &mut buf).is_ok());

        buf.set_position(0);
        let line = buf.lines().skip(2).next().unwrap().unwrap();

        let x = line[20..28].to_string();
        let y = line[28..36].to_string();
        let z = line[36..44].to_string();

        let vx = line[44..52].to_string();
        let vy = line[52..60].to_string();
        let vz = line[60..68].to_string();

        for position in vec![x, y, z] {
            let parts: Vec<_> = position.splitn(2, '.').collect();
            assert_eq!(parts[0].len(), 4);
            assert_eq!(parts[1].len(), 3);
        }

        for velocity in vec![vx, vy, vz] {
            let parts: Vec<_> = velocity.splitn(2, '.').collect();
            assert_eq!(parts[0].len(), 3);
            assert_eq!(parts[1].len(), 4);
        }
    }
}
