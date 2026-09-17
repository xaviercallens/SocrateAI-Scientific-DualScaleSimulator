use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TachyonConfig {
    pub n_grid: usize,
    pub l_box: f64,
    pub t_max: f64,
    pub dt: f64,
    pub v0: f64,
    pub t0: f64,
    pub gamma: f64,
    pub gauge_coupling: f64,
    pub output_csv: String,
    pub output_json: String,
    #[serde(default = "default_seed")]
    pub seed: u64,
}

fn default_seed() -> u64 {
    42
}

impl Default for TachyonConfig {
    fn default() -> Self {
        Self {
            n_grid: 128,
            l_box: 10.0,
            t_max: 20.0,
            dt: 0.02,
            v0: 1.0,
            t0: 1.0,
            gamma: 0.5,
            gauge_coupling: 0.4,
            output_csv: "tachyon_condensation_telemetry.csv".to_string(),
            output_json: "tachyon_condensation_summary.json".to_string(),
            seed: 42,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TachyonPointRecord {
    pub t: f64,
    pub x: f64,
    pub tachyon_field: f64,
    pub tachyon_velocity: f64,
    pub tachyon_gradient: f64,
    pub gauge_field: f64,
    pub potential: f64,
    pub energy_density: f64,
    pub is_soliton_core: bool,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TachyonSummary {
    pub total_steps: usize,
    pub final_time: f64,
    pub soliton_width: f64,
    pub soliton_center: f64,
    pub soliton_peak_energy: f64,
    pub total_energy: f64,
    pub seed: u64,
    /// Tier X: Numerically computed as Q_initial == Q_final where Q = (sign(phi(L)) - sign(phi(-L)))/2.
    /// Not a K-theory class; a discrete-grid sign change diagnostic for the kink configuration.
    pub k_theory_charge_conserved: bool,
    /// Tier X: Count of adjacent sign changes in final t_field. Discrete-grid defect count, not Grothendieck rank.
    pub grothendieck_defect_rank: i64,
    /// Tier X: Cross-check: |Q_final| == defect_count. Both are numerical discretization diagnostics.
    pub rr_charge_match: bool,
}

pub struct TachyonSimulator {
    pub config: TachyonConfig,
    pub t_field: Vec<f64>,
    pub t_dot: Vec<f64>,
    pub a_field: Vec<f64>,
    pub a_dot: Vec<f64>,
    pub dx: f64,
}

/// Compute net kink charge: Q = (sign(phi[last]) - sign(phi[0])) / 2.
/// Tier X: discrete-grid diagnostic, not a true K-theory charge.
#[inline]
pub fn net_kink_charge(field: &[f64]) -> i64 {
    if field.is_empty() {
        return 0;
    }
    let sign_first = if field[0] > 0.0 { 1 } else { -1 };
    let sign_last = if field[field.len() - 1] > 0.0 { 1 } else { -1 };
    ((sign_last - sign_first) / 2) as i64
}

/// Count adjacent sign changes in field.
/// Tier X: discrete-grid defect count, not Grothendieck rank.
#[inline]
pub fn count_sign_changes(field: &[f64]) -> i64 {
    if field.len() < 2 {
        return 0;
    }
    let mut count = 0i64;
    for i in 0..field.len() - 1 {
        if (field[i] > 0.0) != (field[i + 1] > 0.0) {
            count += 1;
        }
    }
    count
}

impl TachyonSimulator {
    pub fn new(config: TachyonConfig) -> Self {
        let n = config.n_grid;
        let dx = config.l_box / (n as f64);

        let mut t_field = vec![0.0; n];
        let t_dot = vec![0.0; n];
        let mut a_field = vec![0.0; n];
        let a_dot = vec![0.0; n];

        // Initial anti-symmetric perturbation to seed the kink soliton (Sen conjecture)
        for i in 0..n {
            let x = (i as f64 - n as f64 * 0.5) * dx;
            t_field[i] = 0.1 * (x / 2.0).tanh();
            a_field[i] = 0.05 * (-0.5 * (x / 1.5).powi(2)).exp();
        }

        Self {
            config,
            t_field,
            t_dot,
            a_field,
            a_dot,
            dx,
        }
    }

    /// Sen tachyon potential V(T) = V0 / cosh(T / T0)
    #[inline(always)]
    pub fn potential(&self, t_val: f64) -> f64 {
        let arg = t_val / self.config.t0;
        self.config.v0 / arg.cosh().max(1e-12)
    }

    /// dV/dT = - (V0 / T0) * sinh(T / T0) / cosh^2(T / T0)
    #[inline(always)]
    pub fn potential_deriv(&self, t_val: f64) -> f64 {
        let arg = t_val / self.config.t0;
        let cosh_val = arg.cosh().max(1e-12);
        let sinh_val = arg.sinh();
        -(self.config.v0 / self.config.t0) * sinh_val / (cosh_val * cosh_val)
    }

    pub fn step(&mut self) {
        let n = self.config.n_grid;
        let dt = self.config.dt;
        let dx = self.dx;
        let gamma = self.config.gamma;
        let g2 = self.config.gauge_coupling * self.config.gauge_coupling;

        let mut new_t = vec![0.0; n];
        let mut new_t_dot = vec![0.0; n];
        let mut new_a = vec![0.0; n];
        let mut new_a_dot = vec![0.0; n];

        for i in 0..n {
            let ip1 = if i + 1 < n { i + 1 } else { n - 1 };
            let im1 = if i > 0 { i - 1 } else { 0 };

            // Spatial Laplacian
            let lap_t = (self.t_field[ip1] - 2.0 * self.t_field[i] + self.t_field[im1]) / (dx * dx);
            let lap_a = (self.a_field[ip1] - 2.0 * self.a_field[i] + self.a_field[im1]) / (dx * dx);

            let t_val = self.t_field[i];
            let a_val = self.a_field[i];

            let dv_dt = self.potential_deriv(t_val);

            // EOM: d^2 T / dt^2 = lap_T - dV/dT - gamma * d_t T - g^2 A^2 T
            let d2t_dt2 = lap_t - dv_dt - gamma * self.t_dot[i] - g2 * a_val * a_val * t_val;
            // Gauge EOM: d^2 A / dt^2 = lap_A - gamma_A * d_t A - g^2 T^2 A
            let d2a_dt2 = lap_a - 0.2 * self.a_dot[i] - g2 * t_val * t_val * a_val;

            new_t_dot[i] = self.t_dot[i] + dt * d2t_dt2;
            new_t[i] = self.t_field[i] + dt * new_t_dot[i];

            new_a_dot[i] = self.a_dot[i] + dt * d2a_dt2;
            new_a[i] = self.a_field[i] + dt * new_a_dot[i];
        }

        self.t_field = new_t;
        self.t_dot = new_t_dot;
        self.a_field = new_a;
        self.a_dot = new_a_dot;
    }

    pub fn run_simulation(&mut self) -> (TachyonSummary, Vec<TachyonPointRecord>) {
        println!("--- Starting Tachyon Condensation (Sen Soliton) Simulation ---");
        let n = self.config.n_grid;
        let total_steps = (self.config.t_max / self.config.dt).ceil() as usize;
        let snapshot_interval = total_steps.max(10) / 10;
        let mut records = Vec::new();

        // Capture initial charge before simulation starts
        let q_initial = net_kink_charge(&self.t_field);

        for step in 0..=total_steps {
            let t = step as f64 * self.config.dt;
            if step > 0 {
                self.step();
            }

            if step % snapshot_interval == 0 || step == total_steps {
                for i in 0..n {
                    let x = (i as f64 - n as f64 * 0.5) * self.dx;
                    let ip1 = if i + 1 < n { i + 1 } else { n - 1 };
                    let im1 = if i > 0 { i - 1 } else { 0 };
                    let grad_t = (self.t_field[ip1] - self.t_field[im1]) / (2.0 * self.dx);

                    let pot = self.potential(self.t_field[i]);
                    let e_dens = 0.5 * self.t_dot[i].powi(2)
                        + 0.5 * grad_t.powi(2)
                        + pot
                        + 0.5 * self.a_dot[i].powi(2);

                    let is_core = grad_t.abs() > 0.8 && e_dens > 0.5;

                    records.push(TachyonPointRecord {
                        t,
                        x,
                        tachyon_field: self.t_field[i],
                        tachyon_velocity: self.t_dot[i],
                        tachyon_gradient: grad_t,
                        gauge_field: self.a_field[i],
                        potential: pot,
                        energy_density: e_dens,
                        is_soliton_core: is_core,
                    });
                }
            }
        }

        // Analyze final soliton profile
        let mut max_grad = 0.0;
        let mut peak_x = 0.0;
        let mut peak_energy = 0.0;
        let mut total_e = 0.0;

        for i in 0..n {
            let x = (i as f64 - n as f64 * 0.5) * self.dx;
            let ip1 = if i + 1 < n { i + 1 } else { n - 1 };
            let im1 = if i > 0 { i - 1 } else { 0 };
            let grad_t = (self.t_field[ip1] - self.t_field[im1]) / (2.0 * self.dx);
            let pot = self.potential(self.t_field[i]);
            let e_dens = 0.5 * self.t_dot[i].powi(2) + 0.5 * grad_t.powi(2) + pot;
            total_e += e_dens * self.dx;

            if grad_t.abs() > max_grad {
                max_grad = grad_t.abs();
                peak_x = x;
                peak_energy = e_dens;
            }
        }

        // Compute charge diagnostics
        let q_final = net_kink_charge(&self.t_field);
        let defect_count = count_sign_changes(&self.t_field);

        let summary = TachyonSummary {
            total_steps,
            final_time: self.config.t_max,
            soliton_width: 2.0 * self.config.t0 / max_grad.max(0.1),
            soliton_center: peak_x,
            soliton_peak_energy: peak_energy,
            total_energy: total_e,
            seed: self.config.seed,
            k_theory_charge_conserved: q_initial == q_final,
            grothendieck_defect_rank: defect_count,
            rr_charge_match: q_final.abs() == defect_count,
        };

        println!(
            "Tachyon Condensation complete. Sen Kink Soliton Peak Energy: {:.4} at x={:.2}, Width={:.4}",
            peak_energy, peak_x, summary.soliton_width
        );

        (summary, records)
    }

    pub fn export_results(&self, summary: &TachyonSummary, records: &[TachyonPointRecord]) -> std::io::Result<()> {
        let csv_file = File::create(Path::new(&self.config.output_csv))?;
        let mut writer = BufWriter::new(csv_file);

        writeln!(writer, "t,x,tachyon_field,tachyon_velocity,tachyon_gradient,gauge_field,potential,energy_density,is_soliton_core")?;
        for r in records {
            writeln!(
                writer,
                "{:.4},{:.4},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{}",
                r.t, r.x, r.tachyon_field, r.tachyon_velocity, r.tachyon_gradient, r.gauge_field, r.potential, r.energy_density, r.is_soliton_core as i32
            )?;
        }
        writer.flush()?;

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
        let mut config1 = TachyonConfig::default();
        config1.seed = 42;
        config1.t_max = 0.1; // Short run for test speed
        config1.dt = 0.02;

        let mut config2 = TachyonConfig::default();
        config2.seed = 42;
        config2.t_max = 0.1;
        config2.dt = 0.02;

        let mut sim1 = TachyonSimulator::new(config1);
        let (summary1, _) = sim1.run_simulation();

        let mut sim2 = TachyonSimulator::new(config2);
        let (summary2, _) = sim2.run_simulation();

        let s1 = serde_json::to_string(&summary1).unwrap();
        let s2 = serde_json::to_string(&summary2).unwrap();
        assert_eq!(s1, s2, "Same seed should produce identical summaries");
    }

    #[test]
    fn test_different_seed_differs() {
        // Negative control: different seed would be a no-op for tachyon (no stochasticity)
        // but the seed field itself should differ
        let mut config1 = TachyonConfig::default();
        config1.seed = 42;
        config1.t_max = 0.1;
        config1.dt = 0.02;

        let mut config2 = TachyonConfig::default();
        config2.seed = 123;
        config2.t_max = 0.1;
        config2.dt = 0.02;

        let mut sim1 = TachyonSimulator::new(config1);
        let (summary1, _) = sim1.run_simulation();

        let mut sim2 = TachyonSimulator::new(config2);
        let (summary2, _) = sim2.run_simulation();

        assert_ne!(summary1.seed, summary2.seed, "Different seeds should be recorded");
    }

    #[test]
    fn test_charge_diagnostic_positive_control() {
        // Positive control: hand-constructed field with consistent boundary signs
        // should show charge conservation and matching defect count
        let field = vec![-0.1, -0.05, 0.0, 0.05, 0.1];
        let q = net_kink_charge(&field);
        let defects = count_sign_changes(&field);

        // sign(-0.1) = -1, sign(0.1) = +1 => Q = (+1 - (-1))/2 = +1
        assert_eq!(q, 1, "Expected charge +1 for boundary signs [-,-..+,+]");
        // Sign changes: -0.1→-0.05 (no), -0.05→0.0 (no), 0.0→0.05 (yes), 0.05→0.1 (no) = 1 change
        assert_eq!(defects, 1, "Expected 1 sign change in kink field");
        assert_eq!(q.abs(), defects, "Charge magnitude should match defect count");
    }

    #[test]
    fn test_charge_diagnostic_negative_control() {
        // Negative control: hand-constructed field whose boundary signs differ
        // at t=0 vs t=final (simulated by two different fields)
        // Field 1: boundary signs consistent (left negative)
        let field_initial = vec![-0.1, -0.05, 0.0, 0.05, 0.1];
        let q_initial = net_kink_charge(&field_initial);

        // Field 2: boundary signs flipped (right now negative)
        let field_final = vec![0.1, 0.05, 0.0, -0.05, -0.1];
        let q_final = net_kink_charge(&field_final);

        // Initial: sign(-0.1)=-1, sign(0.1)=+1 => Q=+1
        // Final: sign(0.1)=+1, sign(-0.1)=-1 => Q=-1
        assert_eq!(q_initial, 1);
        assert_eq!(q_final, -1);
        assert_ne!(q_initial, q_final, "Charge should change when boundary signs flip");
    }
}
