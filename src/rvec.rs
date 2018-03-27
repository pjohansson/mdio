use std::default::Default;
use std::f64;
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign};

/// Directions in a carthesian 3-dimensional system.
pub enum Direction {
    X,
    Y,
    Z,
}

#[derive(Clone, Copy, Debug, PartialEq)]
/// Shorthand for a 3-vector.
pub struct RVec {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, PartialEq)]
pub enum ParseRVecError {
    MissingValues,
    ParseFloatError,
}

impl RVec {
    /// Return the absolute distance between two vectors.
    pub fn distance(&self, other: &RVec) -> f64 {
        f64::sqrt(
            (self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2),
        )
    }

    /// Return the cylindrical distance between two vectors and along an input `Direction`
    /// as a (dr, dh) tuple. For the height difference, the second value (the other)
    /// is subtracted from the first (self).
    pub fn distance_cylindrical(&self, other: &RVec, dir: Direction) -> (f64, f64) {
        let dr = |dx: f64, dy: f64| f64::sqrt(dx.powi(2) + dy.powi(2));

        match dir {
            Direction::X => (dr(self.y - other.y, self.z - other.z), self.x - other.x),
            Direction::Y => (dr(self.x - other.x, self.z - other.z), self.y - other.y),
            Direction::Z => (dr(self.x - other.x, self.y - other.y), self.z - other.z),
        }
    }

    pub fn from_fixed(input: &str, length: usize) -> Result<RVec, ParseRVecError> {
        use std::str::from_utf8;

        if input.trim().is_empty() {
            return Err(ParseRVecError::MissingValues);
        }

        let mut iter = input
            .as_bytes()
            .chunks(length)
            .map(|chunk| from_utf8(chunk).map_err(|_| ParseRVecError::ParseFloatError))
            .map(|s| {
                s?.trim()
                    .parse::<f64>()
                    .map_err(|_| ParseRVecError::ParseFloatError)
            });

        Ok(RVec {
            x: iter.next().ok_or(ParseRVecError::MissingValues)??,
            y: iter.next().ok_or(ParseRVecError::MissingValues)??,
            z: iter.next().ok_or(ParseRVecError::MissingValues)??,
        })
    }

    pub fn from_whitespace(input: &str) -> Result<RVec, ParseRVecError> {
        let mut iter = input.split_whitespace().map(|s| {
            s.parse::<f64>()
                .map_err(|_| ParseRVecError::ParseFloatError)
        });

        Ok(RVec {
            x: iter.next().ok_or(ParseRVecError::MissingValues)??,
            y: iter.next().ok_or(ParseRVecError::ParseFloatError)??,
            z: iter.next().ok_or(ParseRVecError::ParseFloatError)??,
        })
    }

    pub fn pbc_multiply(self, nx: usize, ny: usize, nz: usize) -> RVec {
        RVec {
            x: self.x * (nx as f64),
            y: self.y * (ny as f64),
            z: self.z * (nz as f64),
        }
    }

    pub fn to_tuple(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }
}

impl Default for RVec {
    fn default() -> RVec {
        RVec {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

impl Add for RVec {
    type Output = RVec;

    fn add(self, other: RVec) -> Self::Output {
        RVec {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl AddAssign for RVec {
    fn add_assign(&mut self, other: RVec) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl Sub for RVec {
    type Output = RVec;

    fn sub(self, other: RVec) -> Self::Output {
        RVec {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl SubAssign for RVec {
    fn sub_assign(&mut self, other: RVec) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl Neg for RVec {
    type Output = RVec;

    fn neg(self) -> Self::Output {
        RVec {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_rvec_from_fixed_string() {
        assert_eq!(
            RVec::from_fixed("123", 1),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );
        assert_eq!(
            RVec::from_fixed(" 11230", 2),
            Ok(RVec {
                x: 1.0,
                y: 12.0,
                z: 30.0,
            })
        );
        assert_eq!(
            RVec::from_fixed(" 12 3", 2),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );
        assert_eq!(
            RVec::from_fixed("1.02.03.0", 3),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );
        assert_eq!(
            RVec::from_fixed(" 1.0 2.0 3.0", 4),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );
        assert_eq!(
            RVec::from_fixed("123   ", 1),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );

        assert_eq!(RVec::from_fixed("", 1), Err(ParseRVecError::MissingValues));
        assert_eq!(
            RVec::from_fixed("  ", 1),
            Err(ParseRVecError::MissingValues)
        );
        assert_eq!(
            RVec::from_fixed("12 ", 1),
            Err(ParseRVecError::ParseFloatError)
        );
        assert_eq!(
            RVec::from_fixed(" 12", 1),
            Err(ParseRVecError::ParseFloatError)
        );
        assert_eq!(
            RVec::from_fixed(" 123", 1),
            Err(ParseRVecError::ParseFloatError)
        );
        assert_eq!(
            RVec::from_fixed("1s3", 1),
            Err(ParseRVecError::ParseFloatError)
        );
        assert_eq!(
            RVec::from_fixed("1s23", 1),
            Err(ParseRVecError::ParseFloatError)
        );
    }

    #[test]
    fn parse_rvec_from_whitespace_separated_string() {
        assert_eq!(
            RVec::from_whitespace("1 2 3"),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );
        assert_eq!(
            RVec::from_whitespace("    1  2         3"),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );
        assert_eq!(
            RVec::from_whitespace("\t1\t2\t3\t"),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );
        assert_eq!(
            RVec::from_whitespace("1\t2 3"),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );
        assert_eq!(
            RVec::from_whitespace("1 2 3 4 5"),
            Ok(RVec {
                x: 1.0,
                y: 2.0,
                z: 3.0,
            })
        );

        assert_eq!(
            RVec::from_whitespace(""),
            Err(ParseRVecError::MissingValues)
        );
        assert_eq!(
            RVec::from_whitespace("   "),
            Err(ParseRVecError::MissingValues)
        );
        assert_eq!(
            RVec::from_whitespace("   2 3"),
            Err(ParseRVecError::ParseFloatError)
        );
        assert_eq!(
            RVec::from_whitespace("1 s 2 3"),
            Err(ParseRVecError::ParseFloatError)
        );
        assert_eq!(
            RVec::from_whitespace("1,2,3"),
            Err(ParseRVecError::ParseFloatError)
        );
    }

    #[test]
    fn rvec_to_tuple() {
        let (x, y, z) = (1.0, 2.0, 3.0);
        let r = RVec { x, y, z };
        assert_eq!((x, y, z), r.to_tuple());
    }

    #[test]
    fn add_rvec_operator() {
        let r1 = RVec {
            x: 0.0,
            y: 1.0,
            z: 2.0,
        };
        let r2 = RVec {
            x: 3.0,
            y: 4.0,
            z: 5.0,
        };
        assert_eq!(
            r1 + r2,
            RVec {
                x: 3.0,
                y: 5.0,
                z: 7.0,
            }
        );

        let mut r3 = RVec {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        r3 += r1;
        assert_eq!(r3, r1);
    }

    #[test]
    fn sub_rvec_operator() {
        let r1 = RVec {
            x: 0.0,
            y: 1.0,
            z: 2.0,
        };
        let r2 = RVec {
            x: 3.0,
            y: 4.0,
            z: 5.0,
        };
        assert_eq!(
            r2 - r1,
            RVec {
                x: 3.0,
                y: 3.0,
                z: 3.0,
            }
        );

        let mut r3 = r2.clone();
        r3 -= r2;
        assert_eq!(
            r3,
            RVec {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            }
        );
    }

    #[test]
    fn mul_rvec_with_usize() {
        let r = RVec {
            x: 1.0,
            y: 2.0,
            z: 3.0,
        };
        assert_eq!(
            r.pbc_multiply(2, 3, 4),
            RVec {
                x: 2.0,
                y: 6.0,
                z: 12.0,
            }
        );
    }

    #[test]
    fn neg_rvec_operator() {
        let (x, y, z) = (1.0, 2.0, 3.0);
        let r = RVec { x, y, z };
        assert_eq!(
            -r,
            RVec {
                x: -x,
                y: -y,
                z: -z,
            }
        );
    }

    #[test]
    fn distance_between_rvecs() {
        let r1 = RVec {
            x: 1.0,
            y: 2.0,
            z: 3.0,
        };
        let r2 = RVec {
            x: 7.0,
            y: 11.0,
            z: 13.0,
        };
        let dr2 = (r1.x - r2.x).powi(2) + (r1.y - r2.y).powi(2) + (r1.z - r2.z).powi(2);

        assert_eq!(r1.distance(&r2), dr2.sqrt());
    }

    #[test]
    fn cylindrical_distances_between_rvecs() {
        let r1 = RVec {
            x: 1.0,
            y: 2.0,
            z: 3.0,
        };
        let r2 = RVec {
            x: 7.0,
            y: 11.0,
            z: 13.0,
        };

        let dr_z = ((r1.x - r2.x).powi(2) + (r1.y - r2.y).powi(2)).sqrt();
        let dh_z = r1.z - r2.z;
        assert_eq!(r1.distance_cylindrical(&r2, Direction::Z), (dr_z, dh_z));

        let dr_y = ((r1.x - r2.x).powi(2) + (r1.z - r2.z).powi(2)).sqrt();
        let dh_y = r1.y - r2.y;
        assert_eq!(r1.distance_cylindrical(&r2, Direction::Y), (dr_y, dh_y));

        let dr_x = ((r1.z - r2.z).powi(2) + (r1.y - r2.y).powi(2)).sqrt();
        let dh_x = r1.x - r2.x;
        assert_eq!(r1.distance_cylindrical(&r2, Direction::X), (dr_x, dh_x));
    }

    #[test]
    fn rvec_default_is_origo() {
        let origo = RVec {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };

        assert_eq!(origo, RVec::default());
    }
}
