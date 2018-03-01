use conf::{Atom, Conf, get_or_insert_atom_and_residue};
use rvec::{RVec, ParseRVecError};

use std::io::{BufRead, BufReader, Read};

struct Line<'a> {
    // residue_number: usize,
    residue_name: &'a str,
    atom_name: &'a str,
    // atom_number: usize,
    position: RVec,
    velocity: Option<RVec>,
}

#[derive(Debug, Fail)]
pub enum IoError {
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

pub fn read_gromos87_conf<R: Read>(reader: R) -> Result<Conf, IoError> {
    let mut buf_reader = BufReader::new(reader);
    let mut buf = String::new();

    buf_reader.read_line(&mut buf).map_err(|_| IoError::Utf8Error(1))?;
    let title = buf.trim().to_string();
    buf.clear();

    buf_reader.read_line(&mut buf).map_err(|_| IoError::Utf8Error(1))?;
    let num_atoms = buf.trim().parse::<usize>().map_err(|_| IoError::NumAtomsError)?;
    buf.clear();

    let mut residues = Vec::new();
    let mut atoms = Vec::new();

    for i in 0..num_atoms {
        buf_reader.read_line(&mut buf).map_err(|_| IoError::Utf8Error(2 + i))?;

        let atom_line = parse_atom_line(&buf).map_err(|_| IoError::LineError(2 + i))?;
        let (residue, atom) = get_or_insert_atom_and_residue(
                atom_line.residue_name, atom_line.atom_name, &mut residues
            ).map_err(|_| IoError::LineError(2 + i))?;;

        atoms.push(Atom {
            name: atom,
            residue,
            position: atom_line.position,
            velocity: atom_line.velocity,
        });

        buf.clear();
    }

    buf_reader.read_line(&mut buf).map_err(|_| IoError::Utf8Error(3 + num_atoms))?;
    let size = RVec::from_whitespace(&buf).expect("could not read box size");

    Ok(Conf {
        title,
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
        assert_eq!(line.position, RVec { x: 1000.001, y: 2000.002, z: 3000.003 });
        assert_eq!(line.velocity, None);

        let s = "    12RES12ATO150001 100.01  200.02  300.03  400.04  500.05  600.06 ";
        let line = parse_atom_line(s).unwrap();
        // assert_eq!(line.residue_number, 1);
        // assert_eq!(line.atom_number, 50001);
        assert_eq!(line.residue_name, "2RES1");
        assert_eq!(line.atom_name, "2ATO1");
        assert_eq!(line.position, RVec { x: 100.01, y: 200.02, z: 300.03 });
        assert_eq!(line.velocity, Some(RVec { x: 400.04, y: 500.05, z: 600.06 }));
    }

    #[test]
    fn read_correct_file() {
        let title = "A title";

        let res1_name = "RES1";
        let res2_name = "RES2";

        let atom1_name = "AT1";
        let atom1_pos = RVec { x: 0.0, y: 1.0, z: 2.0 };
        let atom1_vel = RVec { x: 0.0, y: 0.1, z: 0.3 };
        let atom2_name = "AT2";
        let atom2_pos = RVec { x: 3.0, y: 4.0, z: 5.0 };
        let atom2_vel = RVec { x: 0.3, y: 0.4, z: 0.5 };
        let atom3_name = "AT3";
        let atom3_pos1 = RVec { x: 6.0, y: 7.0, z: 8.0 };
        let atom3_vel1 = RVec { x: 0.6, y: 0.7, z: 0.8 };
        let atom3_pos2 = RVec { x: 9.0, y: 10.0, z: 11.0 };
        let atom3_vel2 = RVec { x: 0.9, y: 1.0, z: 1.1 };

        let size = RVec { x: 10.0, y: 11.0, z: 12.0 };

        let content = format!("\
{}
4
{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.3}{:>8.3}{:>8.3}
{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.3}{:>8.3}{:>8.3}
{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.3}{:>8.3}{:>8.3}
{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.3}{:>8.3}{:>8.3}
{:12.3} {:12.3} {:12.3}
", title,
   1, res1_name, atom1_name, 1, atom1_pos.x, atom1_pos.y, atom1_pos.z,
                                atom1_vel.x, atom1_vel.y, atom1_vel.z,
   1, res1_name, atom2_name, 2, atom2_pos.x, atom2_pos.y, atom2_pos.z,
                                atom2_vel.x, atom2_vel.y, atom2_vel.z,
   2, res2_name, atom3_name, 3, atom3_pos1.x, atom3_pos1.y, atom3_pos1.z,
                                atom3_vel1.x, atom3_vel1.y, atom3_vel1.z,
   3, res2_name, atom3_name, 4, atom3_pos2.x, atom3_pos2.y, atom3_pos2.z,
                                atom3_vel2.x, atom3_vel2.y, atom3_vel2.z,
   size.x, size.y, size.z);

        let conf = read_gromos87_conf(content.as_bytes()).unwrap();

        assert_eq!(conf.title, title);
        assert_eq!(conf.atoms.len(), 4);
        assert_eq!(conf.size, size);

        // Verify that all residues were correctly constructed
        assert_eq!(conf.residues.len(), 2);
        assert_eq!(conf.residues[0].borrow().name, res1_name);
        assert_eq!(conf.residues[1].borrow().name, res2_name);

        assert_eq!(conf.residues[0].borrow().atoms.len(), 2);
        assert_eq!(&*conf.residues[0].borrow().atoms[0].borrow(), atom1_name);
        assert_eq!(&*conf.residues[0].borrow().atoms[1].borrow(), atom2_name);

        assert_eq!(conf.residues[1].borrow().atoms.len(), 1);
        assert_eq!(&*conf.residues[1].borrow().atoms[0].borrow(), atom3_name);

        // Verify that pointers are correctly set from atoms to their residues
        assert!(Rc::ptr_eq(&conf.atoms[0].name, &conf.residues[0].borrow().atoms[0]));
        assert!(Rc::ptr_eq(&conf.atoms[1].residue, &conf.residues[0]));
        assert!(Rc::ptr_eq(&conf.atoms[2].name, &conf.residues[1].borrow().atoms[0]));
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

        let content = format!("{}\n{}\n{}\n{}", "More atom lines than atoms", 1, two_atom_lines, size_line);
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!("{}\n{}\n{}\n{}", "Fewer atom lines than atoms", 3, two_atom_lines, size_line);
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!("{}\n{}\n{}", "No box size line", 2, two_atom_lines);
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!("{}\n{}\n{}1.0 2.0", "Bad box size line", 2, two_atom_lines);
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!("{}\n{}\n{}1.s 2.0 3.0", "Bad box size line", 2, two_atom_lines);
        assert!(read_gromos87_conf(content.as_bytes()).is_err());

        let content = format!("{}\n{}\n{}", "No number of atoms line", two_atom_lines, size_line);
        assert!(read_gromos87_conf(content.as_bytes()).is_err());
    }
}
