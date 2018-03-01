#[derive(Clone, Copy, Debug, PartialEq)]
/// Shorthand for a 3-vector.
pub struct RVec {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(Debug, PartialEq)]
pub enum ParseRVecError {
    MissingValues,
    ParseFloatError,
}

impl RVec {
    pub fn from_fixed(input: &str, length: usize) -> Result<RVec, ParseRVecError> {
        use std::str::from_utf8;

        if input.trim().is_empty() {
            return Err(ParseRVecError::MissingValues);
        }

        let mut iter = input
            .as_bytes()
            .chunks(length)
            .map(|chunk| from_utf8(chunk).map_err(|_| ParseRVecError::ParseFloatError))
            .map(|s| s?.trim().parse::<f32>().map_err(|_| ParseRVecError::ParseFloatError));

        Ok(RVec {
            x: iter.next().ok_or(ParseRVecError::MissingValues)??,
            y: iter.next().ok_or(ParseRVecError::MissingValues)??,
            z: iter.next().ok_or(ParseRVecError::MissingValues)??,
        })
    }

    pub fn from_whitespace(input: &str) -> Result<RVec, ParseRVecError> {
        let mut iter = input
            .split_whitespace()
            .map(|s| s.parse::<f32>().map_err(|_| ParseRVecError::ParseFloatError));

        Ok(RVec {
            x: iter.next().ok_or(ParseRVecError::MissingValues)??,
            y: iter.next().ok_or(ParseRVecError::ParseFloatError)??,
            z: iter.next().ok_or(ParseRVecError::ParseFloatError)??,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_rvec_from_fixed_string() {
        assert_eq!(RVec::from_fixed("123", 1), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));
        assert_eq!(RVec::from_fixed(" 11230", 2), Ok(RVec { x: 1.0, y: 12.0, z: 30.0 }));
        assert_eq!(RVec::from_fixed(" 12 3", 2), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));
        assert_eq!(RVec::from_fixed("1.02.03.0", 3), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));
        assert_eq!(RVec::from_fixed(" 1.0 2.0 3.0", 4), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));
        assert_eq!(RVec::from_fixed("123   ", 1), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));

        assert_eq!(RVec::from_fixed("", 1), Err(ParseRVecError::MissingValues));
        assert_eq!(RVec::from_fixed("  ", 1), Err(ParseRVecError::MissingValues));
        assert_eq!(RVec::from_fixed("12 ", 1), Err(ParseRVecError::ParseFloatError));
        assert_eq!(RVec::from_fixed(" 12", 1), Err(ParseRVecError::ParseFloatError));
        assert_eq!(RVec::from_fixed(" 123", 1), Err(ParseRVecError::ParseFloatError));
        assert_eq!(RVec::from_fixed("1s3", 1), Err(ParseRVecError::ParseFloatError));
        assert_eq!(RVec::from_fixed("1s23", 1), Err(ParseRVecError::ParseFloatError));
    }

    #[test]
    fn parse_rvec_from_whitespace_separated_string() {
        assert_eq!(RVec::from_whitespace("1 2 3"), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));
        assert_eq!(RVec::from_whitespace("    1  2         3"), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));
        assert_eq!(RVec::from_whitespace("\t1\t2\t3\t"), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));
        assert_eq!(RVec::from_whitespace("1\t2 3"), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));
        assert_eq!(RVec::from_whitespace("1 2 3 4 5"), Ok(RVec { x: 1.0, y: 2.0, z: 3.0 }));

        assert_eq!(RVec::from_whitespace(""), Err(ParseRVecError::MissingValues));
        assert_eq!(RVec::from_whitespace("   "), Err(ParseRVecError::MissingValues));
        assert_eq!(RVec::from_whitespace("   2 3"), Err(ParseRVecError::ParseFloatError));
        assert_eq!(RVec::from_whitespace("1 s 2 3"), Err(ParseRVecError::ParseFloatError));
        assert_eq!(RVec::from_whitespace("1,2,3"), Err(ParseRVecError::ParseFloatError));
    }
}
