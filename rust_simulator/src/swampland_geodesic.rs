use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SwamplandConfig {
    pub s_max: f64,
    pub ds: f64,
    pub tau_init: (f64, f64),
    pub t_init: (f64, f64),
    pub v_tau: (f64, f64),
    pub v_t: (f64, f64),
    pub max_mode_n: i32,
    pub max_mode_m: i32,
    pub output_csv: String,
    pub output_json: String,
}

impl Default for SwamplandConfig {
    fn default() -> Self {
        Self {
            s_max: 5.0,
            ds: 0.01,
            tau_init: (0.0, 1.0),
            t_init: (0.0, 1.0),
            v_tau: (0.0, 1.0),
            v_t: (0.0, 1.0),
            max_mode_n: 4,
            max_mode_m: 4,
            output_csv: "swampland_geodesic_telemetry.csv".to_string(),
            output_json: "swampland_geodesic_summary.json".to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SwamplandStepRecord {
    pub s: f64,
    pub delta_d: f64,
    pub tau1: f64,
    pub tau2: f64,
    pub t1: f64,
    pub t2: f64,
    pub min_mass_gap: f64,
    pub theoretical_bound: f64,
    pub eft_cutoff: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SwamplandSummary {
    pub total_steps: usize,
    pub total_geodesic_distance: f64,
    pub initial_mass_gap: f64,
    pub final_mass_gap: f64,
    pub mass_gap_ratio: f64,
    pub sdc_alpha: f64,
    pub bound_satisfied: bool,
    pub eft_cutoff_broken: bool,
}

pub struct SwamplandSimulator {
    pub config: SwamplandConfig,
}

impl SwamplandSimulator {
    pub fn new(config: SwamplandConfig) -> Self {
        Self { config }
    }

    /// Computes the mass squared of mode (n, m) at (tau, T)
    #[inline(always)]
    pub fn mode_mass_sq(n: i32, m: i32, tau1: f64, tau2: f64, t2: f64) -> f64 {
        let n_f = n as f64;
        let m_f = m as f64;
        let numerator = (n_f - m_f * tau1).powi(2) + (m_f * tau2).powi(2);
        let denominator = t2 * tau2;
        numerator / denominator.max(1e-12)
    }

    /// Computes the lowest positive mass gap in the KK/winding spectrum
    pub fn compute_min_mass_gap(&self, tau1: f64, tau2: f64, t2: f64) -> f64 {
        let mut min_m = f64::MAX;
        let max_n = self.config.max_mode_n;
        let max_m = self.config.max_mode_m;

        for n in -max_n..=max_n {
            for m in -max_m..=max_m {
                if n == 0 && m == 0 {
                    continue;
                }
                let m_sq = Self::mode_mass_sq(n, m, tau1, tau2, t2);
                let mass = m_sq.max(0.0).sqrt();
                if mass > 1e-10 && mass < min_m {
                    min_m = mass;
                }
            }
        }
        min_m
    }

    /// Integrates the non-linear hyperbolic geodesic equations:
    /// d^2 x / ds^2 = 2/y (dx/ds)(dy/ds)
    /// d^2 y / ds^2 = 1/y ((dy/ds)^2 - (dx/ds)^2)
    pub fn run_simulation(&self) -> (SwamplandSummary, Vec<SwamplandStepRecord>) {
        println!("--- Starting Swampland Distance Geodesic Flow Simulation ---");
        let ds = self.config.ds;
        let n_steps = (self.config.s_max / ds).ceil() as usize;

        let mut tau1 = self.config.tau_init.0;
        let mut tau2 = self.config.tau_init.1;
        let mut u_tau1 = self.config.v_tau.0;
        let mut u_tau2 = self.config.v_tau.1;

        let mut t1 = self.config.t_init.0;
        let mut t2 = self.config.t_init.1;
        let mut u_t1 = self.config.v_t.0;
        let mut u_t2 = self.config.v_t.1;

        let mut delta_d = 0.0;
        let mut records = Vec::with_capacity(n_steps + 1);

        let alpha_sdc = 1.0 / std::f64::consts::SQRT_2; // 0.7071
        let initial_mass = self.compute_min_mass_gap(tau1, tau2, t2);
        let m0 = initial_mass;

        for step in 0..=n_steps {
            let s = step as f64 * ds;
            let mass_gap = self.compute_min_mass_gap(tau1, tau2, t2);
            let theoretical_bound = m0 * (-alpha_sdc * delta_d).exp();
            let eft_cutoff = 10.0 * (-alpha_sdc * delta_d).exp();

            records.push(SwamplandStepRecord {
                s,
                delta_d,
                tau1,
                tau2,
                t1,
                t2,
                min_mass_gap: mass_gap,
                theoretical_bound,
                eft_cutoff,
            });

            if step == n_steps {
                break;
            }

            // Metric velocity squared on H x H
            let v_sq_tau = (u_tau1 * u_tau1 + u_tau2 * u_tau2) / (tau2 * tau2);
            let v_sq_t = (u_t1 * u_t1 + u_t2 * u_t2) / (t2 * t2);
            let speed = (v_sq_tau + v_sq_t).max(0.0).sqrt();
            delta_d += speed * ds;

            // Geodesic acceleration (Christoffel symbols for Poincaré half-plane)
            let a_tau1 = 2.0 * u_tau1 * u_tau2 / tau2;
            let a_tau2 = (u_tau2 * u_tau2 - u_tau1 * u_tau1) / tau2;

            let a_t1 = 2.0 * u_t1 * u_t2 / t2;
            let a_t2 = (u_t2 * u_t2 - u_t1 * u_t1) / t2;

            // Stiff symplectic / Velocity-Verlet update
            tau1 += ds * u_tau1 + 0.5 * ds * ds * a_tau1;
            tau2 = (tau2 + ds * u_tau2 + 0.5 * ds * ds * a_tau2).max(1e-4);
            u_tau1 += ds * a_tau1;
            u_tau2 += ds * a_tau2;

            t1 += ds * u_t1 + 0.5 * ds * ds * a_t1;
            t2 = (t2 + ds * u_t2 + 0.5 * ds * ds * a_t2).max(1e-4);
            u_t1 += ds * a_t1;
            u_t2 += ds * a_t2;
        }

        let final_mass = records.last().map(|r| r.min_mass_gap).unwrap_or(0.0);
        let ratio = final_mass / initial_mass.max(1e-12);
        let bound_ok = records.iter().all(|r| r.min_mass_gap <= r.theoretical_bound * 1.5 + 1e-6);

        let summary = SwamplandSummary {
            total_steps: n_steps,
            total_geodesic_distance: delta_d,
            initial_mass_gap: initial_mass,
            final_mass_gap: final_mass,
            mass_gap_ratio: ratio,
            sdc_alpha: alpha_sdc,
            bound_satisfied: bound_ok,
            eft_cutoff_broken: final_mass < 0.1,
        };

        println!("Swampland Geodesic complete. Total Δd: {:.4}, Mass Gap Ratio: {:.6e}", delta_d, ratio);
        (summary, records)
    }

    pub fn export_results(&self, summary: &SwamplandSummary, records: &[SwamplandStepRecord]) -> std::io::Result<()> {
        let csv_file = File::create(Path::new(&self.config.output_csv))?;
        let mut writer = BufWriter::new(csv_file);

        writeln!(writer, "s,delta_d,tau1,tau2,t1,t2,min_mass_gap,theoretical_bound,eft_cutoff")?;
        for r in records {
            writeln!(
                writer,
                "{:.4},{:.6},{:.6},{:.6},{:.6},{:.6},{:.8},{:.8},{:.8}",
                r.s, r.delta_d, r.tau1, r.tau2, r.t1, r.t2, r.min_mass_gap, r.theoretical_bound, r.eft_cutoff
            )?;
        }
        writer.flush()?;

        let json_file = File::create(Path::new(&self.config.output_json))?;
        serde_json::to_writer_pretty(json_file, summary)?;

        Ok(())
    }
}
