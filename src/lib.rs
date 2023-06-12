use rand::prelude::*;
use rand_distr::StandardNormal;

pub struct Statistics {
    mean: f64,
    var: f64,
    skew: f64,
    ex_kurt: f64,
    n_sample: u32,
}

fn variance(samples: &[f64], bias: bool) -> f64 {
    /* biased estimator, as the same default value of numpy.var */
    let n = samples.len() as f64;
    let mu = samples.iter().sum::<f64>() / n;
    let mut var = samples.iter().map(|x| (*x - mu) * (*x - mu)).sum::<f64>();
    var /= match bias {
        true => n,
        false => n - 1.,
    };
    var
}

fn skewness(samples: &[f64], bias: bool) -> f64 {
    /* biased estimator, as the same default value of scipy.stats.kurtosis */
    let n = samples.len() as f64;
    let mu = samples.iter().sum::<f64>() / n;
    let (mut m3, mut s3) = (0f64, 0f64);

    for v in samples {
        let shift = (v - mu);
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

fn kurtosis(samples: &[f64], bias: bool) -> f64 {
    /* biased estimator, as the same default value of scipy.stats.kurtosis */
    let n = samples.len() as f64;
    let mu = samples.iter().sum::<f64>() / n;
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
        false => (n - 1.) / (n - 2.) / (n - 3.) * ((n + 1.) * m4 / m2 / m2 - 3. * (n - 1.)),
    }
}

pub fn cubic_transformation_sampling(tgt_stats: &Statistics, n_scenario: usize) -> Vec<f64> {
    // generating standard normal distribution samples
    let mut rng = thread_rng();
    let xs: Vec<f64> = StandardNormal
        .sample_iter(&mut rng)
        .take(n_scenario)
        .collect();

    let mut ex = [0f64; 12];
    for idx in 1..=ex.len() {
        ex[idx] = xs.iter().map(|x| x.powf(idx as f64)).sum::<f64>() / tgt_stats.n_sample as f64;
    }

    // to generate samples Z with zero mean, unit variance, the same skew and kurt with tgt.
    let mut z_moments = [0f64, 1f64, tgt_stats.skew, tgt_stats.ex_kurt + 3.];

    // using least square error algorithm to get params
    let mut params = [0f64, 1f64, 0f64, 0f64];

    let problem = CubicProblem {
        p: Vector4::new(0., 1., 0., 0.),
    };
    let (_result, report) = LevenbergMarquardt::new().minimize(problem);

    let mut scenarios = Vec::with_capacity(n_scenario);
    // scenarios
    //     .iter_mut()
    //     .map(|x| params[0] + params[1] * x + params[2] * x.powf(2.) + params[3] * x.powf(3.));

    scenarios
}

fn residuals(params: [f64; 4], ex: [f64; 12], ez: [f64; 4]) -> [f64; 4] {
    let (a, b, c, d) = (params[0], params[1], params[2], params[3]);
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
        + (12. * a * c * d * d + 6. * b * b * d * d + 12. * b * c * c * d + c * c * c * c) * ex[7]
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

    errors
}

struct CubicProblem {
    // holds current value of the n parameters
    p: Vector4<f64>,
}

impl LeastSquaresProblem<f64, U4, U4> for CubicProblem {
    type ParameterStorage = Owned<f64, U4>;
    type ResidualStorage = Owned<f64, U4>;
    type JacobianStorage = Owned<f64, U4, U4>;
    fn set_params(&mut self, p: &VectorN<f64, U4>) {
        self.p.copy_from(p);
        // do common calculations for residuals and the Jacobian here
    }

    fn params(&self) -> VectorN<f64, U2> {
        self.p
    }

    fn residuals(&self) -> Option<Vector2<f64>> {
        let [x, y] = [self.p.x, self.p.y];
        Some(Vector2::new(x * x + y - 11., x + y * y - 7.))
    }
    fn jacobian(&self) -> Option<Matrix2<f64>> {
        let [x, y] = [self.p.x, self.p.y];

        // first row of Jacobian, derivatives of first residual
        let d1_x = 2. * x;
        let d1_y = 1.;
        let d2_x = 1.;
        let d2_y = 2. * y;
        Some(Matrix2::new(d1_x, d1_y, d2_x, d2_y))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_variance() {
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
        assert!((variance(&samples, true) - 0.4333647219842476) < 1e-3);
        assert!((variance(&samples, false) - 0.45617339156236586) < 1e-3);
    }

    #[test]
    fn test_skewness() {
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
        assert!((skewness(&samples, true) - 0.3027475074365746) < 1e-3);
        assert!((skewness(&samples, false) - 0.327868632598646) < 1e-3)
    }

    #[test]
    fn test_kurtosis() {
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
        assert!((kurtosis(&samples, true) - -0.6516395832609172) < 1e-3);
        assert!((kurtosis(&samples, false) - -0.47713788797747014) < 1e-3);
    }
}
