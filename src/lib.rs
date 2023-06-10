use rand::prelude::*;
use rand_distr::StandardNormal;

pub struct Statistics {
    mean: f64,
    var: f64,
    skew: f64,
    ex_kurt: f64,
    n_sample: u32,
}

// impl Statistics {
//     pub fn moments(&self, bias: bool) -> Option<[f64; 4]> {
//         let mut moments = [0f64; 4];
//         if self.n_sample <= 2 {
//             return None;
//         }
//         moments[0] = self.n_sample as f64 * self.mean;
//
//         if bias == true {
//             // biased estimator
//             moments[1] = self.n_sample as f64 * self.var;
//         } else {
//             // unbiased estimator
//             moments[1] = (self.n_sample - 1) as f64 * self.var;
//         }
//
//         Some(moments)
//     }
// }

pub fn cubic_transformation(tgt_stats: &Statistics, n_scenario: usize) -> Vec<f64> {
    let mut scenarios = Vec::with_capacity(n_scenario);

    // generating standard normal distribution samples
    let mut rng = thread_rng();
    let samples: Vec<f64> = StandardNormal
        .sample_iter(&mut rng)
        .take(n_scenario)
        .collect();

    let mut ex = [0f64; 12];
    for idx in 0..ex.len() {
        ex[idx] = samples
            .iter()
            .map(|x| x.powf((idx + 1) as f64))
            .sum::<f64>()
            / tgt_stats.n_sample as f64;
    }

    // to generate samples Y with zero mean, unit variance, the same skew and kurt with tgt.
    let mut y_moments = [0f64, 1f64, tgt_stats.skew, tgt_stats.ex_kurt + 3.];

    // using least square error algorithm to get params
    let mut params = [0f64, 1f64, 0f64, 0f64];

    scenarios
        .iter_mut()
        .map(|x| params[0] + params[1] * x + params[2] * x.powf(2.) + params[3] * x.powf(3.));

    scenarios
}

fn cubic_errors(params: [f64; 4], EX: [f64; 12], EY: [f64; 4]) -> [f64; 4] {
    let (a, b, c, d) = (params[0], params[1], params[2], params[3]);
    let (mut v1, mut v2, mut v3, mut v4) = (0f64, 0f64, 0f64, 0f64);

    v1 = a + b * EX[0] + c * EX[1] + d * EX[2] - EY[0];

    v2 = (d * d * EX[5]
        + 2. * c * d * EX[4]
        + (2. * b * d + c * c) * EX[3]
        + (2. * a * d + 2. * b * c) * EX[2]
        + (2. * a * c + b * b) * EX[1]
        + 2. * a * b * EX[0]
        + a * a
        - EY[1]);

    v3 = ((d * d * d) * EX[8]
        + (3. * c * d * d) * EX[7]
        + (3. * b * d * d + 3. * c * c * d) * EX[6]
        + (3. * a * d * d + 6. * b * c * d + c * c * c) * EX[5]
        + (6. * a * c * d + 3. * b * b * d + 3. * b * c * c) * EX[4]
        + (a * (6. * b * d + 3. * c * c) + 3. * b * b * c) * EX[3]
        + (3. * a * a * d + 6. * a * b * c + b * b * b) * EX[2]
        + (3. * a * a * c + 3. * a * b * b) * EX[1]
        + 3. * a * a * b * EX[0]
        + a * a * a
        - EY[2]);

    v4 = ((d * d * d * d) * EX[11]
        + (4. * c * d * d * d) * EX[10]
        + (4. * b * d * d * d + 6. * c * c * d * d) * EX[9]
        + (4. * a * d * d * d + 12. * b * c * d * d + 4. * c * c * c * d) * EX[8]
        + (12. * a * c * d * d + 6. * b * b * d * d + 12. * b * c * c * d + c * c * c * c) * EX[7]
        + (a * (12. * b * d * d + 12. * c * c * d) + 12. * b * b * c * d + 4. * b * c * c * c)
            * EX[6]
        + (6. * a * a * d * d
            + a * (24. * b * c * d + 4. * c * c * c)
            + 4. * b * b * b * d
            + 6. * b * b * c * c)
            * EX[5]
        + (12. * a * a * c * d + a * (12. * b * b * d + 12. * b * c * c) + 4. * b * b * b * c)
            * EX[4]
        + (a * a * (12. * b * d + 6. * c * c) + 12. * a * b * b * c + b * b * b * b) * EX[3]
        + (4. * a * a * a * d + 12. * a * a * b * c + 4. * a * b * b * b) * EX[2]
        + (4. * a * a * a * c + 6. * a * a * b * b) * EX[1]
        + (4. * a * a * a * b) * EX[0]
        + a * a * a * a
        - EY[3]);

    [v1, v2, v3, v4]
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
