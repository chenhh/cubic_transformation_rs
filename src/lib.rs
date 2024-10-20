use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
// extern crate nalgebra as na;
use nalgebra as na;
use debug_print::debug_println;
use na::{Matrix4, Owned, Vector4, U4};
use rand::prelude::*;
use rand_distr::StandardNormal;

pub struct Statistics {
    mean: f64,
    var: f64,
    skew: f64,
    ex_kurt: f64,
}
pub fn mean(samples: &[f64]) -> f64 {
    let n = samples.len();
    match n {
        0 => 0.0,
        _ => samples.iter().sum::<f64>() / n as f64,
    }
}

pub fn variance(samples: &[f64], bias: bool) -> f64 {
    /* biased estimator, as the same default value of numpy.var */
    let n = samples.len();
    match n {
        0 => 0.0,
        1 if !bias => 0.0,
        _ => {
            let n = n as f64;
            let mu = mean(samples);
            let mut var = samples.iter().map(|x| (*x - mu) * (*x - mu)).sum::<f64>();
            var /= match bias {
                true => n,
                false => n - 1.,
            };
            var
        }
    }
}

pub fn skewness(samples: &[f64], bias: bool) -> f64 {
    /* biased estimator, as the same default value of scipy.stats.kurtosis */
    let n = samples.len();
    match n {
        0 => 0.0,
        1..=2 if !bias => 0.0,
        _ => {
            let n = n as f64;
            let mu = mean(samples);
            let (mut m3, mut s3) = (0f64, 0f64);

            for v in samples {
                let shift = v - mu;
                let shift2 = shift * shift;
                s3 += shift2;
                m3 += shift2 * shift;
            }

            m3 /= n;
            s3 /= n;
            let res = m3 / s3.powf(1.5);
            match bias {
                true => res,
                false => res * ((n - 1.) * n).sqrt() / (n - 2.),
            }
        }
    }
}

pub fn kurtosis(samples: &[f64], bias: bool) -> f64 {
    /* biased estimator, as the same default value of scipy.stats.kurtosis */
    let n = samples.len();
    match n {
        0 => 0.0,
        1..=3 if !bias => 0.0,
        _ => {
            let n = n as f64;
            let mu = mean(samples);
            let (mut m4, mut m2) = (0f64, 0f64);

            for v in samples {
                let shift2 = (v - mu) * (v - mu);
                m4 += shift2 * shift2;
                m2 += shift2;
            }

            m4 /= n;
            m2 /= n;
            match bias {
                true => m4 / m2 / m2 - 3.,
                false => {
                    (n - 1.) / (n - 2.) / (n - 3.) * ((n + 1.) * m4 / m2 / m2 - 3. * (n - 1.))
                }
            }
        }
    }
}

pub fn cubic_transformation_sampling_iter3(
    tgt_stats: &Statistics,
    n_scenario: usize,
    max_stats_mse: f64,
) -> Option<Vec<f64>> {
    cubic_transformation_sampling(tgt_stats, n_scenario, max_stats_mse, 3)
}

pub fn cubic_transformation_sampling_iter5(
    tgt_stats: &Statistics,
    n_scenario: usize,
    max_stats_mse: f64,
) -> Option<Vec<f64>> {
    cubic_transformation_sampling(tgt_stats, n_scenario, max_stats_mse, 5)
}

pub fn cubic_transformation_sampling_iter10(
    tgt_stats: &Statistics,
    n_scenario: usize,
    max_stats_mse: f64,
) -> Option<Vec<f64>> {
    cubic_transformation_sampling(tgt_stats, n_scenario, max_stats_mse, 10)
}

pub fn cubic_transformation_sampling(
    tgt_stats: &Statistics,
    n_scenario: usize,
    max_stats_mse: f64,         // max moment the least square error
    max_cubic_iteration: usize, // max resampling times
) -> Option<Vec<f64>> {
    // sometimes the samples don't converge well because of the bad random samples xs,
    // and it requires resampling the xs until the error converges.

    let mut scenarios: Option<Vec<f64>> = None;

    for _cubic_iter in 0..max_cubic_iteration {
        // generating standard normal distribution samples
        let mut rng = thread_rng();
        let xs: Vec<f64> = StandardNormal
            .sample_iter(&mut rng)
            .take(n_scenario)
            .collect();

        // 1 to 12th moments of the samples
        let ex: [f64; 12] = (1..=12)
            .map(|pdx| xs.iter().map(|&x| x.powi(pdx)).sum::<f64>() / n_scenario as f64)
            .collect::<Vec<_>>()
            .try_into()
            .expect("can't generating 1~12th moments.");

        // to generate samples Z with zero mean, unit variance, the same skew and kurt with tgt.
        let z_moments = [0., 1., tgt_stats.skew, tgt_stats.ex_kurt + 3.];

        // using the least square error algorithm to get params
        // let mut params = [0f64, 1f64, 0f64, 0f64];

        let problem = CubicProblem {
            p: Vector4::new(0., 1., 0., 0.),
            ex,
            ez: z_moments,
        };
        let (result, _) = LevenbergMarquardt::new().minimize(problem);
        let cubic_params = result.p;

        let zs: Vec<f64> = xs
            .iter()
            .map(|x| {
                cubic_params[0]
                    + cubic_params[1] * x
                    + cubic_params[2] * x * x
                    + cubic_params[3] * x * x * x
            })
            .collect();

        scenarios = Some(
            zs.iter()
                .map(|z| tgt_stats.mean + tgt_stats.var.sqrt() * z)
                .collect(),
        );
        let samples = scenarios.clone()?;
        let stats_mse = statistics_mse(tgt_stats, &samples);
        debug_println!(
            "cub_iter[{}] inside scenario statistics:{}, {}, {}, {}, mse:{}",
            _cubic_iter,
            mean(&samples),
            variance(&samples, true),
            skewness(&samples, true),
            kurtosis(&samples, true),
            stats_mse
        );
        if stats_mse < max_stats_mse {
            break;
        }
    }

    scenarios
}

fn statistics_mse(tgt_stats: &Statistics, scenarios: &[f64]) -> f64 {
    let mut errors = [0.; 4];
    errors[0] = tgt_stats.mean - mean(scenarios);
    errors[1] = tgt_stats.var - variance(scenarios, true);
    errors[2] = tgt_stats.skew - skewness(scenarios, true);
    errors[3] = tgt_stats.ex_kurt - kurtosis(scenarios, true);

    let res = errors.iter().map(|x| x * x).sum::<f64>();
    res
}

#[derive(Debug)]
struct CubicProblem {
    // holds current value of the n parameters
    // the 4 dimensions are (x,y,z,w)
    p: Vector4<f64>,
    ex: [f64; 12],
    ez: [f64; 4],
}

impl LeastSquaresProblem<f64, U4, U4> for CubicProblem {
    type ResidualStorage = Owned<f64, U4>;
    type JacobianStorage = Owned<f64, U4, U4>;
    type ParameterStorage = Owned<f64, U4>;
    fn set_params(&mut self, p: &Vector4<f64>) {
        self.p.copy_from(p);
    }

    fn params(&self) -> Vector4<f64> {
        self.p
    }

    fn residuals(&self) -> Option<Vector4<f64>> {
        let [a, b, c, d] = [self.p.x, self.p.y, self.p.z, self.p.w];
        let ex = &self.ex;
        let ez = &self.ez;
        let mut errors = [0f64; 4];
        errors[0] = a + b * ex[0] + c * ex[1] + d * ex[2] - ez[0];
        errors[1] = d * d * ex[5]
            + 2. * c * d * ex[4]
            + (2. * b * d + c * c) * ex[3]
            + (2. * a * d + 2. * b * c) * ex[2]
            + (2. * a * c + b * b) * ex[1]
            + 2. * a * b * ex[0]
            + a * a
            - ez[1];

        errors[2] = d * d * d * ex[8]
            + 3. * c * d * d * ex[7]
            + (3. * b * d * d + 3. * c * c * d) * ex[6]
            + (3. * a * d * d + 6. * b * c * d + c * c * c) * ex[5]
            + (6. * a * c * d + 3. * b * b * d + 3. * b * c * c) * ex[4]
            + (a * (6. * b * d + 3. * c * c) + 3. * b * b * c) * ex[3]
            + (3. * a * a * d + 6. * a * b * c + b * b * b) * ex[2]
            + (3. * a * a * c + 3. * a * b * b) * ex[1]
            + 3. * a * a * b * ex[0]
            + a * a * a
            - ez[2];

        errors[3] = (d * d * d * d) * ex[11]
            + (4. * c * d * d * d) * ex[10]
            + (4. * b * d * d * d + 6. * c * c * d * d) * ex[9]
            + (4. * a * d * d * d + 12. * b * c * d * d + 4. * c * c * c * d) * ex[8]
            + (12. * a * c * d * d + 6. * b * b * d * d + 12. * b * c * c * d + c * c * c * c)
            * ex[7]
            + (a * (12. * b * d * d + 12. * c * c * d) + 12. * b * b * c * d + 4. * b * c * c * c)
            * ex[6]
            + (6. * a * a * d * d
            + a * (24. * b * c * d + 4. * c * c * c)
            + 4. * b * b * b * d
            + 6. * b * b * c * c)
            * ex[5]
            + (12. * a * a * c * d + a * (12. * b * b * d + 12. * b * c * c) + 4. * b * b * b * c)
            * ex[4]
            + (a * a * (12. * b * d + 6. * c * c) + 12. * a * b * b * c + b * b * b * b) * ex[3]
            + (4. * a * a * a * d + 12. * a * a * b * c + 4. * a * b * b * b) * ex[2]
            + (4. * a * a * a * c + 6. * a * a * b * b) * ex[1]
            + (4. * a * a * a * b) * ex[0]
            + a * a * a * a
            - ez[3];

        Some(Vector4::new(errors[0], errors[1], errors[2], errors[3]))
    }
    fn jacobian(&self) -> Option<Matrix4<f64>> {
        let [a, b, c, d] = [self.p.x, self.p.y, self.p.z, self.p.w];
        let ex = &self.ex;

        // first row of Jacobian, derivatives of first residual
        let dxx = 1.;
        let dxy = ex[0];
        let dxz = ex[1];
        let dxw = ex[2];
        let dyx = 2. * a + 2. * b * ex[0] + 2. * c * ex[1] + 2. * d * ex[2];
        let dyy = 2. * a * ex[0] + 2. * b * ex[1] + 2. * c * ex[2] + 2. * d * ex[3];
        let dyz = 2. * a * ex[1] + 2. * b * ex[2] + 2. * c * ex[3] + 2. * d * ex[4];
        let dyw = 2. * a * ex[2] + 2. * b * ex[3] + 2. * c * ex[4] + 2. * d * ex[5];
        let dzx = 3. * a * a
            + 6. * a * b * ex[0]
            + 6. * c * d * ex[4]
            + 3. * d * d * ex[5]
            + ex[1] * (6. * a * c + 3. * b * b)
            + ex[2] * (6. * a * d + 6. * b * c)
            + ex[3] * (6. * b * d + 3. * c * c);
        let dzy = 3. * a * a * ex[0]
            + 6. * a * b * ex[1]
            + 6. * c * d * ex[5]
            + 3. * d * d * ex[6]
            + ex[2] * (6. * a * c + 3. * b * b)
            + ex[3] * (6. * a * d + 6. * b * c)
            + ex[4] * (6. * b * d + 3. * c * c);
        let dzz = 3. * a * a * ex[1]
            + 6. * a * b * ex[2]
            + 6. * c * d * ex[6]
            + 3. * d * d * ex[7]
            + ex[3] * (6. * a * c + 3. * b * b)
            + ex[4] * (6. * a * d + 6. * b * c)
            + ex[5] * (6. * b * d + 3. * c * c);
        let dzw = 3. * a * a * ex[2]
            + 6. * a * b * ex[3]
            + 6. * c * d * ex[7]
            + 3. * d * d * ex[8]
            + ex[4] * (6. * a * c + 3. * b * b)
            + ex[5] * (6. * a * d + 6. * b * c)
            + ex[6] * (6. * b * d + 3. * c * c);
        let dwx = 4. * a * a * a
            + 12. * a * a * b * ex[0]
            + 12. * c * d * d * ex[7]
            + 4. * d * d * d * ex[8]
            + ex[1] * (12. * a * a * c + 12. * a * b * b)
            + ex[2] * (12. * a * a * d + 24. * a * b * c + 4. * b * b * b)
            + ex[3] * (2. * a * (12. * b * d + 6. * c * c) + 12. * b * b * c)
            + ex[4] * (24. * a * c * d + 12. * b * b * d + 12. * b * c * c)
            + ex[5] * (12. * a * d * d + 24. * b * c * d + 4. * c * c * c)
            + ex[6] * (12. * b * d * d + 12. * c * c * d);
        let dwy = 4. * a * a * a * ex[0]
            + 12. * a * a * b * ex[1]
            + 12. * c * d * d * ex[8]
            + 4. * d * d * d * ex[9]
            + ex[2] * (12. * a * a * c + 12. * a * b * b)
            + ex[3] * (12. * a * a * d + 24. * a * b * c + 4. * b * b * b)
            + ex[4] * (a * (24. * b * d + 12. * c * c) + 12. * b * b * c)
            + ex[5] * (24. * a * c * d + 12. * b * b * d + 12. * b * c * c)
            + ex[6] * (12. * a * d * d + 24. * b * c * d + 4. * c * c * c)
            + ex[7] * (12. * b * d * d + 12. * c * c * d);
        let dwz = 4. * a * a * a * ex[1]
            + 12. * a * a * b * ex[2]
            + 12. * c * d * d * ex[9]
            + 4. * d * d * d * ex[10]
            + ex[3] * (12. * a * a * c + 12. * a * b * b)
            + ex[4] * (12. * a * a * d + 24. * a * b * c + 4. * b * b * b)
            + ex[5] * (a * (24. * b * d + 12. * c * c) + 12. * b * b * c)
            + ex[6] * (24. * a * c * d + 12. * b * b * d + 12. * b * c * c)
            + ex[7] * (12. * a * d * d + 24. * b * c * d + 4. * c * c * c)
            + ex[8] * (12. * b * d * d + 12. * c * c * d);
        let dww = 4. * a * a * a * ex[2]
            + 12. * a * a * b * ex[3]
            + 12. * c * d * d * ex[10]
            + 4. * d * d * d * ex[11]
            + ex[4] * (12. * a * a * c + 12. * a * b * b)
            + ex[5] * (12. * a * a * d + 24. * a * b * c + 4. * b * b * b)
            + ex[6] * (a * (24. * b * d + 12. * c * c) + 12. * b * b * c)
            + ex[7] * (24. * a * c * d + 12. * b * b * d + 12. * b * c * c)
            + ex[8] * (12. * a * d * d + 24. * b * c * d + 4. * c * c * c)
            + ex[9] * (12. * b * d * d + 12. * c * c * d);
        Some(Matrix4::new(
            dxx, dxy, dxz, dxw, dyx, dyy, dyz, dyw, dzx, dzy, dzz, dzw, dwx, dwy, dwz, dww,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mean() {
        let zero_sample: [f64; 0] = [];
        assert!((mean(&zero_sample) - 0.) < 1e-3);
    }

    #[test]
    fn test_variance() {
        let zero_sample: [f64; 0] = [];
        let one_sample = [1., ];
        let samples = [
            6.37717487e-01,
            -4.13245759e-01,
            -6.35436489e-01,
            -4.24825535e-03,
            4.14607296e-01,
            9.88113278e-04,
            1.51558792e+00,
            4.81379973e-02,
            -9.51924640e-02,
            1.88047272e+00,
            3.17862022e-01,
            6.92952755e-01,
            1.42847418e+00,
            9.95428104e-01,
            7.56869101e-01,
            3.83725359e-01,
            -2.84257360e-02,
            5.59426799e-01,
            1.46734559e+00,
            7.57219266e-01,
        ];
        assert!((variance(&zero_sample, false) - 0.) < 1e-3);
        assert!((variance(&zero_sample, true) - 0.) < 1e-3);
        assert!((variance(&one_sample, false) - 0.) < 1e-3);
        assert!((variance(&one_sample, true) - 1.) < 1e-3);
        assert!((variance(&samples, true) - 0.4333647219842476) < 1e-3);
        assert!((variance(&samples, false) - 0.45617339156236586) < 1e-3);
    }

    #[test]
    fn test_skewness() {
        let zero_sample: [f64; 0] = [];


        let samples = [
            6.37717487e-01,
            -4.13245759e-01,
            -6.35436489e-01,
            -4.24825535e-03,
            4.14607296e-01,
            9.88113278e-04,
            1.51558792e+00,
            4.81379973e-02,
            -9.51924640e-02,
            1.88047272e+00,
            3.17862022e-01,
            6.92952755e-01,
            1.42847418e+00,
            9.95428104e-01,
            7.56869101e-01,
            3.83725359e-01,
            -2.84257360e-02,
            5.59426799e-01,
            1.46734559e+00,
            7.57219266e-01,
        ];
        assert!((skewness(&zero_sample, false) - 0.) < 1e-3);
        assert!((skewness(&zero_sample, true) - 0.) < 1e-3);
        assert!((skewness(&samples, true) - 0.3027475074365746) < 1e-3);
        assert!((skewness(&samples, false) - 0.327868632598646) < 1e-3)
    }

    #[test]
    fn test_kurtosis() {
        let zero_sample: [f64; 0] = [];
        let samples = [
            6.37717487e-01,
            -4.13245759e-01,
            -6.35436489e-01,
            -4.24825535e-03,
            4.14607296e-01,
            9.88113278e-04,
            1.51558792e+00,
            4.81379973e-02,
            -9.51924640e-02,
            1.88047272e+00,
            3.17862022e-01,
            6.92952755e-01,
            1.42847418e+00,
            9.95428104e-01,
            7.56869101e-01,
            3.83725359e-01,
            -2.84257360e-02,
            5.59426799e-01,
            1.46734559e+00,
            7.57219266e-01,
        ];
        assert!((kurtosis(&zero_sample, false) - 0.) < 1e-3);
        assert!((kurtosis(&zero_sample, true) - 0.) < 1e-3);
        assert!((kurtosis(&samples, true) - -0.6516395832609172) < 1e-3);
        assert!((kurtosis(&samples, false) - -0.47713788797747014) < 1e-3);
    }

    #[test]
    fn test_cubic() {
        // skew的絕對值大於1.5之後的樣本很難收斂
        // 如果是近似於常態分佈的樣本很容易收斂
        let max_stats_mse = 0.01;

        let mut rng = thread_rng();
        for idx in 0..1000 {
            let tgt = Statistics {
                mean: rng.gen(),
                var: rng.gen(),
                skew: rng.gen(),
                ex_kurt: rng.gen(),
            };
            let scenarios =
                cubic_transformation_sampling_iter10(&tgt, 1000, max_stats_mse).unwrap();
            let stats = statistics_mse(&tgt, &scenarios);
            debug_println!("test idx:[{idx}], mse: {stats}");
            assert!(stats < max_stats_mse);
        }
    }
}
