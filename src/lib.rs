use rand::prelude::*;
use rand_distr::StandardNormal;

pub struct Statistics {
    mean: f64,
    var: f64,
    skew: f64,
    ex_kurt: f64,
    n_sample: u32,
}

fn variance(samples: &Vec<f64>) -> f64 {
    /* biased estimator, as the same default value of numpy.var */
    let n = samples.len();
    let mu = samples.iter().sum::<f64>() / n as f64;
    let var = samples.iter().map(|x| (*x - mu).powf(2.)).sum::<f64>() / n as f64;
    var
}

fn skewness(samples: &Vec<f64>) -> f64 {
    /* biased estimator, as the same default value of scipy.stats.kurtosis */
    let n = samples.len();
    let mu = samples.iter().sum::<f64>() / n as f64;
    let (mut m3, mut s3) = (0f64, 0f64);

    for idx in 0..n {
        m3 += (samples[idx] - mu).powf(3.);
        s3 += (samples[idx] - mu).powf(2.);
    }
    m3 /= n as f64;
    s3 /= n as f64;
    m3 / s3.powf(1.5)
}

fn kurtosis(samples: &Vec<f64>) -> f64 {
    /* biased estimator, as the same default value of scipy.stats.kurtosis */
    let n = samples.len();
    let mu = samples.iter().sum::<f64>() / n as f64;
    let (mut m4, mut s4) = (0f64, 0f64);

    for idx in 0..n {
        m4 += (samples[idx] - mu).powf(4.);
        s4 += (samples[idx] - mu).powf(2.);
    }
    m4 /= n as f64;
    s4 /= n as f64;
    m4 / s4.powf(2.) - 3.
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

    let mut scenarios = Vec::with_capacity(n_scenario);
    // scenarios
    //     .iter_mut()
    //     .map(|x| params[0] + params[1] * x + params[2] * x.powf(2.) + params[3] * x.powf(3.));

    scenarios
}

fn residuals(params: [f64; 4], EX: [f64; 12], EZ: [f64; 4]) -> [f64; 4] {
    let (a, b, c, d) = (params[0], params[1], params[2], params[3]);
    let (mut r1, mut r2, mut r3, mut r4) = (0f64, 0f64, 0f64, 0f64);

    r1 = a + b * EX[0] + c * EX[1] + d * EX[2] - EZ[0];

    r2 = d * d * EX[5]
        + 2. * c * d * EX[4]
        + (2. * b * d + c * c) * EX[3]
        + (2. * a * d + 2. * b * c) * EX[2]
        + (2. * a * c + b * b) * EX[1]
        + 2. * a * b * EX[0]
        + a * a
        - EZ[1];

    r3 = d * d * d * EX[8]
        + 3. * c * d * d * EX[7]
        + (3. * b * d * d + 3. * c * c * d) * EX[6]
        + (3. * a * d * d + 6. * b * c * d + c * c * c) * EX[5]
        + (6. * a * c * d + 3. * b * b * d + 3. * b * c * c) * EX[4]
        + (a * (6. * b * d + 3. * c * c) + 3. * b * b * c) * EX[3]
        + (3. * a * a * d + 6. * a * b * c + b * b * b) * EX[2]
        + (3. * a * a * c + 3. * a * b * b) * EX[1]
        + 3. * a * a * b * EX[0]
        + a * a * a
        - EZ[2];

    r4 = (d * d * d * d) * EX[11]
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
        - EZ[3];

    [r1, r2, r3, r4]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_variance() {
        let samples: Vec<f64> = (0..10).map(|x| x as f64).collect();
        let res = variance(&samples);
        assert!((res - 9.166666) < 1e-3)
    }

    #[test]
    fn test_skewness() {
        let samples: Vec<f64> = (0..10).map(|x| x as f64).collect();
        let res = skewness(&samples);
        assert!((res - 0.) < 1e-3)
    }

    #[test]
    fn test_kurtosis() {
        let samples: Vec<f64> = (0..10).map(|x| x as f64).collect();
        let res = kurtosis(&samples);
        assert!((res - -1.2242424242424244) < 1e-3)
    }
}
