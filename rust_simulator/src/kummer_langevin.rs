use rand::prelude::*;
use rand_distr::{Normal, Distribution};
use serde::{Serialize, Deserialize};
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

fn default_seed() -> u64 {
    42
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimConfig {
    pub grid_size: usize,
    pub t_max: f64,
    pub dt: f64,
    pub kappa: f64,
    pub gamma: f64,
    pub lambda: f64,
    pub v0: f64,
    pub t_crit: f64,
    pub t_init: f64,
    pub alpha_cooling: f64,
    pub epsilon_kummer: f64,
    pub lattice_a: f64,
    pub output_csv: String,
    pub output_json: String,
    #[serde(default)]
    pub use_supergravity: bool,
    #[serde(default = "default_seed")]
    pub seed: u64,
}

impl Default for SimConfig {
    fn default() -> Self {
        Self {
            grid_size: 48,
            t_max: 25.0,
            dt: 0.05,
            kappa: 1.0,
            gamma: 0.8,
            lambda: 1.0,
            v0: 1.5,
            t_crit: 1.0,
            t_init: 2.5,
            alpha_cooling: 0.35,
            epsilon_kummer: 0.25,
            lattice_a: 1.5,
            output_csv: "kummer_langevin_pointcloud.csv".to_string(),
            output_json: "kummer_langevin_summary.json".to_string(),
            use_supergravity: false,
            seed: 42,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SimulationSummary {
    pub grid_size: usize,
    pub total_steps: usize,
    pub final_time: f64,
    pub final_temperature: f64,
    pub final_string_count: usize,
    pub final_wall_pixel_count: usize,
    pub mean_field_norm: f64,
    pub mean_energy_density: f64,
    /// Tier X: Computed as final_temp < t_crit && mean_norm > 0.2; numerical diagnostic.
    pub symmetry_broken: bool,
    pub point_cloud_size: usize,
    pub seed: u64,
}

pub struct KummerLangevinSimulator {
    pub config: SimConfig,
    pub phi1: Vec<f64>,
    pub phi2: Vec<f64>,
    pub pi1: Vec<f64>,
    pub pi2: Vec<f64>,
    pub n: usize,
}

impl KummerLangevinSimulator {
    pub fn new(config: SimConfig) -> Self {
        let n = config.grid_size;
        let mut rng = StdRng::seed_from_u64(config.seed);
        let normal = Normal::new(0.0, 0.05).unwrap();

        let mut phi1 = vec![0.0; n * n];
        let mut phi2 = vec![0.0; n * n];
        let pi1 = vec![0.0; n * n];
        let pi2 = vec![0.0; n * n];

        // Initial symmetric random noise around zero
        for i in 0..(n * n) {
            phi1[i] = normal.sample(&mut rng);
            phi2[i] = normal.sample(&mut rng);
        }

        Self {
            config,
            phi1,
            phi2,
            pi1,
            pi2,
            n,
        }
    }

    #[inline(always)]
    fn idx(&self, i: usize, j: usize) -> usize {
        i * self.n + j
    }

    pub fn temperature(&self, t: f64) -> f64 {
        self.config.t_init / (1.0 + self.config.alpha_cooling * t)
    }

    pub fn vev_squared(&self, temp: f64) -> f64 {
        if temp < self.config.t_crit {
            self.config.v0 * self.config.v0 * (1.0 - temp / self.config.t_crit)
        } else {
            0.0
        }
    }

    pub fn potential_derivatives(&self, p1: f64, p2: f64, vev2: f64) -> (f64, f64, f64) {
        let rho2 = p1 * p1 + p2 * p2;
        let d_pot = rho2 - vev2;
        // Mexican hat component
        let v_mex = 0.25 * self.config.lambda * d_pot * d_pot;
        let dv1_mex = self.config.lambda * d_pot * p1;
        let dv2_mex = self.config.lambda * d_pot * p2;

        // Kummer 16-fold discrete symmetry breaking (4x4 orbifold fixed points)
        let k = 4.0 * PI / self.config.lattice_a;
        let sin_kp1 = (k * p1).sin();
        let sin_kp2 = (k * p2).sin();
        let cos_kp1 = (k * p1).cos();
        let cos_kp2 = (k * p2).cos();

        let v_kummer = self.config.epsilon_kummer * (2.0 - cos_kp1 - cos_kp2);
        let dv1_kummer = self.config.epsilon_kummer * k * sin_kp1;
        let dv2_kummer = self.config.epsilon_kummer * k * sin_kp2;

        (v_mex + v_kummer, dv1_mex + dv1_kummer, dv2_mex + dv2_kummer)
    }

    pub fn supergravity_potential_derivatives(&self, p1: f64, p2: f64) -> (f64, f64, f64) {
        // Upper half plane Poincare moduli coordinate: tau = x + i*y, y > 0
        // Map field p2 -> y (positive axio-dilaton imaginary part)
        let y = (p2.abs() + 0.2).max(0.05);
        let x = p1;

        // Tree-level GVW flux superpotential: W = a0 - tau * b0
        let a0 = 2.0;
        let b0 = 1.0;
        let w_re = a0 - x * b0;
        let w_im = -y * b0;
        let w_sq = w_re * w_re + w_im * w_im;

        // Covariant derivative D_tau W = -b0 + (i / 2y) * W
        let dw_re = -b0 - w_im / (2.0 * y);
        let dw_im = w_re / (2.0 * y);
        let dw_sq = dw_re * dw_re + dw_im * dw_im;

        // F-term potential V_F = e^K (K^{tau, taubar} |D_tau W|^2 - 3 |W|^2)
        let exp_k = 1.0 / (2.0 * y);
        let v_f = exp_k * (4.0 * y * y * dw_sq - 3.0 * w_sq);
        let v_lift = 0.5 / (y * y);
        let v_tot = v_f + v_lift;

        // Finite differences for forces
        let eps = 1e-5;
        let w_re_p = a0 - (x + eps) * b0;
        let dw_re_p = -b0 - w_im / (2.0 * y);
        let dw_im_p = w_re_p / (2.0 * y);
        let v_xp = (1.0 / (2.0 * y)) * (4.0 * y * y * (dw_re_p * dw_re_p + dw_im_p * dw_im_p) - 3.0 * (w_re_p * w_re_p + w_im * w_im)) + v_lift;

        let w_re_m = a0 - (x - eps) * b0;
        let dw_re_m = -b0 - w_im / (2.0 * y);
        let dw_im_m = w_re_m / (2.0 * y);
        let v_xm = (1.0 / (2.0 * y)) * (4.0 * y * y * (dw_re_m * dw_re_m + dw_im_m * dw_im_m) - 3.0 * (w_re_m * w_re_m + w_im * w_im)) + v_lift;

        let dv1 = (v_xp - v_xm) / (2.0 * eps);

        let y_p = y + eps;
        let y_m = (y - eps).max(0.01);
        let w_im_yp = -y_p * b0;
        let dw_re_yp = -b0 - w_im_yp / (2.0 * y_p);
        let dw_im_yp = w_re / (2.0 * y_p);
        let v_yp = (1.0 / (2.0 * y_p)) * (4.0 * y_p * y_p * (dw_re_yp * dw_re_yp + dw_im_yp * dw_im_yp) - 3.0 * (w_re * w_re + w_im_yp * w_im_yp)) + 0.5 / (y_p * y_p);

        let w_im_ym = -y_m * b0;
        let dw_re_ym = -b0 - w_im_ym / (2.0 * y_m);
        let dw_im_ym = w_re / (2.0 * y_m);
        let v_ym = (1.0 / (2.0 * y_m)) * (4.0 * y_m * y_m * (dw_re_ym * dw_re_ym + dw_im_ym * dw_im_ym) - 3.0 * (w_re * w_re + w_im_ym * w_im_ym)) + 0.5 / (y_m * y_m);

        let dv2 = (v_yp - v_ym) / (2.0 * eps);

        (v_tot, dv1, dv2)
    }

    pub fn identify_nearest_kummer_vacuum(&self, p1: f64, p2: f64) -> usize {
        let a = self.config.lattice_a;
        let mut best_id = 0;
        let mut min_dist_sq = f64::MAX;

        let mut id = 0;
        for i in 0..4 {
            for j in 0..4 {
                let center_x = (i as f64 - 1.5) * (a / 2.0);
                let center_y = (j as f64 - 1.5) * (a / 2.0);
                let d2 = (p1 - center_x).powi(2) + (p2 - center_y).powi(2);
                if d2 < min_dist_sq {
                    min_dist_sq = d2;
                    best_id = id;
                }
                id += 1;
            }
        }
        best_id
    }

    pub fn step(&mut self, t: f64, rng: &mut StdRng) {
        let n = self.n;
        let dt = self.config.dt;
        let gamma = self.config.gamma;
        let kappa = self.config.kappa;
        let temp = self.temperature(t);
        let vev2 = self.vev_squared(temp);

        let noise_std = (2.0 * gamma * temp / dt).sqrt();
        let normal = Normal::new(0.0, noise_std).unwrap();

        let mut new_phi1 = vec![0.0; n * n];
        let mut new_phi2 = vec![0.0; n * n];
        let mut new_pi1 = vec![0.0; n * n];
        let mut new_pi2 = vec![0.0; n * n];

        for i in 0..n {
            let ip1 = (i + 1) % n;
            let im1 = (i + n - 1) % n;
            for j in 0..n {
                let jp1 = (j + 1) % n;
                let jm1 = (j + n - 1) % n;
                let c = self.idx(i, j);

                let lap1 = self.phi1[self.idx(ip1, j)]
                    + self.phi1[self.idx(im1, j)]
                    + self.phi1[self.idx(i, jp1)]
                    + self.phi1[self.idx(i, jm1)]
                    - 4.0 * self.phi1[c];

                let lap2 = self.phi2[self.idx(ip1, j)]
                    + self.phi2[self.idx(im1, j)]
                    + self.phi2[self.idx(i, jp1)]
                    + self.phi2[self.idx(i, jm1)]
                    - 4.0 * self.phi2[c];

                let (_, dv1, dv2) = if self.config.use_supergravity {
                    self.supergravity_potential_derivatives(self.phi1[c], self.phi2[c])
                } else {
                    self.potential_derivatives(self.phi1[c], self.phi2[c], vev2)
                };

                let eta1 = normal.sample(rng);
                let eta2 = normal.sample(rng);

                let dpi1_dt = kappa * lap1 - dv1 - gamma * self.pi1[c] + eta1;
                let dpi2_dt = kappa * lap2 - dv2 - gamma * self.pi2[c] + eta2;

                new_pi1[c] = self.pi1[c] + dt * dpi1_dt;
                new_pi2[c] = self.pi2[c] + dt * dpi2_dt;

                new_phi1[c] = self.phi1[c] + dt * new_pi1[c];
                new_phi2[c] = self.phi2[c] + dt * new_pi2[c];
            }
        }

        self.phi1 = new_phi1;
        self.phi2 = new_phi2;
        self.pi1 = new_pi1;
        self.pi2 = new_pi2;
    }

    pub fn compute_defects(&self) -> (usize, usize, Vec<f64>, Vec<i32>) {
        let n = self.n;
        let mut string_count = 0;
        let mut wall_pixel_count = 0;
        let mut grad_sq = vec![0.0; n * n];
        let mut vorticity = vec![0; n * n];

        let wall_thresh = 0.65;

        for i in 0..n {
            let ip1 = (i + 1) % n;
            for j in 0..n {
                let jp1 = (j + 1) % n;
                let c = self.idx(i, j);

                let d1x = 0.5 * (self.phi1[self.idx(ip1, j)] - self.phi1[c]);
                let d1y = 0.5 * (self.phi1[self.idx(i, jp1)] - self.phi1[c]);
                let d2x = 0.5 * (self.phi2[self.idx(ip1, j)] - self.phi2[c]);
                let d2y = 0.5 * (self.phi2[self.idx(i, jp1)] - self.phi2[c]);

                let gsq = d1x * d1x + d1y * d1y + d2x * d2x + d2y * d2y;
                grad_sq[c] = gsq;
                if gsq > wall_thresh {
                    wall_pixel_count += 1;
                }

                let theta0 = self.phi2[c].atan2(self.phi1[c]);
                let theta1 = self.phi2[self.idx(ip1, j)].atan2(self.phi1[self.idx(ip1, j)]);
                let theta2 = self.phi2[self.idx(ip1, jp1)].atan2(self.phi1[self.idx(ip1, jp1)]);
                let theta3 = self.phi2[self.idx(i, jp1)].atan2(self.phi2[self.idx(i, jp1)]);

                let mut dtheta = 0.0;
                for (a, b) in [(theta0, theta1), (theta1, theta2), (theta2, theta3), (theta3, theta0)] {
                    let mut diff = b - a;
                    while diff > PI { diff -= 2.0 * PI; }
                    while diff <= -PI { diff += 2.0 * PI; }
                    dtheta += diff;
                }

                let wind = (dtheta / (2.0 * PI)).round() as i32;
                if wind != 0 {
                    string_count += 1;
                    vorticity[c] = wind;
                }
            }
        }

        (string_count, wall_pixel_count, grad_sq, vorticity)
    }

    pub fn run_simulation(&mut self) -> (SimulationSummary, Vec<Vec<f64>>) {
        let total_steps = (self.config.t_max / self.config.dt).ceil() as usize;
        // Derive dynamics RNG seed from config seed to avoid correlation between init and dynamics.
        // XOR with constant avoids zero if config.seed is zero.
        let dynamics_seed = self.config.seed ^ 0x9E37_79B9_7F4A_7C15u64;
        let mut rng = StdRng::seed_from_u64(dynamics_seed);

        println!("--- Starting Kummer Moduli Space Langevin Simulation ---");
        println!("Grid: {}x{}, Steps: {}, t_max: {:.2}, dt: {:.4}", self.n, self.n, total_steps, self.config.t_max, self.config.dt);

        let mut point_cloud_records = Vec::new();
        let snapshot_interval = total_steps.max(10) / 10;

        for step in 0..=total_steps {
            let t = step as f64 * self.config.dt;
            if step > 0 {
                self.step(t, &mut rng);
            }

            if step % snapshot_interval == 0 || step == total_steps {
                let temp = self.temperature(t);
                let vev2 = self.vev_squared(temp);
                let (_, _, grad_sq, vorticity) = self.compute_defects();

                for i in 0..self.n {
                    for j in 0..self.n {
                        let c = self.idx(i, j);
                        let p1 = self.phi1[c];
                        let p2 = self.phi2[c];
                        let (v, _, _) = self.potential_derivatives(p1, p2, vev2);
                        let gsq = grad_sq[c];
                        let e_dens = 0.5 * gsq + v;
                        let vort = vorticity[c] as f64;
                        let vac_id = self.identify_nearest_kummer_vacuum(p1, p2) as f64;

                        point_cloud_records.push(vec![
                            t,
                            i as f64,
                            j as f64,
                            p1,
                            p2,
                            gsq,
                            v,
                            temp,
                            e_dens,
                            vort,
                            vac_id,
                        ]);
                    }
                }
            }
        }

        let final_temp = self.temperature(self.config.t_max);
        let (final_strings, final_walls, final_gsq, _) = self.compute_defects();
        let vev2 = self.vev_squared(final_temp);

        let mut total_norm = 0.0;
        let mut total_energy = 0.0;
        for c in 0..(self.n * self.n) {
            let p1 = self.phi1[c];
            let p2 = self.phi2[c];
            let norm = (p1 * p1 + p2 * p2).sqrt();
            total_norm += norm;
            let (v, _, _) = self.potential_derivatives(p1, p2, vev2);
            total_energy += 0.5 * final_gsq[c] + v;
        }

        let mean_norm = total_norm / (self.n * self.n) as f64;
        let mean_energy = total_energy / (self.n * self.n) as f64;

        let summary = SimulationSummary {
            grid_size: self.n,
            total_steps,
            final_time: self.config.t_max,
            final_temperature: final_temp,
            final_string_count: final_strings,
            final_wall_pixel_count: final_walls,
            mean_field_norm: mean_norm,
            mean_energy_density: mean_energy,
            symmetry_broken: final_temp < self.config.t_crit && mean_norm > 0.2,
            point_cloud_size: point_cloud_records.len(),
            seed: self.config.seed,
        };

        (summary, point_cloud_records)
    }

    pub fn export_results(&self, summary: &SimulationSummary, records: &[Vec<f64>]) -> std::io::Result<()> {
        println!("Exporting {} points to CSV: {}", records.len(), self.config.output_csv);
        let csv_file = File::create(Path::new(&self.config.output_csv))?;
        let mut writer = BufWriter::new(csv_file);

        writeln!(writer, "t,x_pos,y_pos,phi1,phi2,grad_sq,potential,temperature,energy_density,vorticity,attractor_id")?;
        for r in records {
            writeln!(
                writer,
                "{:.4},{:.1},{:.1},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.0},{:.0}",
                r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10]
            )?;
        }
        writer.flush()?;

        println!("Exporting summary JSON: {}", self.config.output_json);
        let json_file = File::create(Path::new(&self.config.output_json))?;
        serde_json::to_writer_pretty(json_file, summary)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_same_seed_identical_summary() {
        // Positive control: same seed should produce bit-identical results
        let mut config1 = SimConfig::default();
        config1.seed = 42;
        config1.grid_size = 16; // Small grid for test speed
        config1.t_max = 0.1;
        config1.dt = 0.05;

        let mut config2 = SimConfig::default();
        config2.seed = 42;
        config2.grid_size = 16;
        config2.t_max = 0.1;
        config2.dt = 0.05;

        let mut sim1 = KummerLangevinSimulator::new(config1);
        let (summary1, _) = sim1.run_simulation();

        let mut sim2 = KummerLangevinSimulator::new(config2);
        let (summary2, _) = sim2.run_simulation();

        let s1 = serde_json::to_string(&summary1).unwrap();
        let s2 = serde_json::to_string(&summary2).unwrap();
        assert_eq!(s1, s2, "Same seed should produce identical summaries");
    }

    #[test]
    fn test_different_seed_differs() {
        // Negative control: different seed should produce different results
        let mut config1 = SimConfig::default();
        config1.seed = 42;
        config1.grid_size = 16;
        config1.t_max = 0.1;
        config1.dt = 0.05;

        let mut config2 = SimConfig::default();
        config2.seed = 123;
        config2.grid_size = 16;
        config2.t_max = 0.1;
        config2.dt = 0.05;

        let mut sim1 = KummerLangevinSimulator::new(config1);
        let (summary1, _) = sim1.run_simulation();

        let mut sim2 = KummerLangevinSimulator::new(config2);
        let (summary2, _) = sim2.run_simulation();

        assert_ne!(summary1.seed, summary2.seed, "Different seeds should be recorded");
        // With high probability, at least one physics field should differ
        assert!(
            summary1.mean_field_norm != summary2.mean_field_norm ||
            summary1.mean_energy_density != summary2.mean_energy_density ||
            summary1.final_string_count != summary2.final_string_count,
            "Different seeds should produce different simulation results"
        );
    }

    #[test]
    fn test_symmetry_broken_computed() {
        // Verify that symmetry_broken is computed, not hardcoded
        let config = SimConfig::default();
        let mut sim = KummerLangevinSimulator::new(config);
        let (summary, _) = sim.run_simulation();

        // symmetry_broken should be computed as final_temp < t_crit && mean_norm > 0.2
        let expected = summary.final_temperature < sim.config.t_crit && summary.mean_field_norm > 0.2;
        assert_eq!(summary.symmetry_broken, expected, "symmetry_broken should match computed condition");
    }
}
