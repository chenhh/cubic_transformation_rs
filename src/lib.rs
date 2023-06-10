pub struct Moment{
    mean: f64,
    var: f64,
    skew: f64,
    ex_kurt: f64
}

pub fn cubic_transformation(target_mom: &Moment, n_scenario: u32) -> Vec<f64>{
    let scenarios = Vec::new();
    scenarios
}


pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
