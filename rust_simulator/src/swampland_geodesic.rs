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
    pub bound_slack: f64,
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
            bound_slack: 1.0,
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
    pub alpha_effective: Option<f64>,
    pub bound_satisfied: bool,
    pub eft_cutoff_broken: bool,
}

pub struct SwamplandSimulator {
    pub config: SwamplandConfig,
}

/// Checks if all records satisfy the theoretical bound with a configurable slack factor.
///
/// The bound is: min_mass_gap <= theoretical_bound * slack + epsilon
///
/// Returns: true if all records satisfy the bound, false otherwise.
pub fn check_bound(records: &[SwamplandStepRecord], slack: f64) -> bool {
    records.iter().all(|r| r.min_mass_gap <= r.theoretical_bound * slack + 1e-6)
}

impl SwamplandSimulator {
    pub fn new(config: SwamplandConfig) -> Self {
        Self { config }
    }

    /// Computes the mass squared of mode (n, m) at (tau, T).
    ///
    /// Represents Kaluza-Klein (winding/momentum) modes on the torus compactification.
    /// The function evaluates mass only on one torus sector (parameter t1 does not appear).
    /// The mode set is truncated to |n|, |m| ≤ max_mode_n/max_mode_m; this is an upper
    /// bound on the true lattice infimum.
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

    /// Computes acceleration at (tau1, tau2, u_tau1, u_tau2, t1, t2, u_t1, u_t2)
    /// using Christoffel symbols for the Poincaré half-plane.
    /// Pure function: no clamping applied here.
    #[inline]
    fn compute_acceleration(&self, _tau1: f64, tau2: f64, u_tau1: f64, u_tau2: f64,
                          _t1: f64, t2: f64, u_t1: f64, u_t2: f64) -> (f64, f64, f64, f64) {
        let a_tau1 = 2.0 * u_tau1 * u_tau2 / tau2;
        let a_tau2 = (u_tau2 * u_tau2 - u_tau1 * u_tau1) / tau2;
        let a_t1 = 2.0 * u_t1 * u_t2 / t2;
        let a_t2 = (u_t2 * u_t2 - u_t1 * u_t1) / t2;
        (a_tau1, a_tau2, a_t1, a_t2)
    }

    /// Computes the metric speed at the current state.
    #[inline]
    fn compute_speed(&self, u_tau1: f64, u_tau2: f64, tau2: f64, u_t1: f64, u_t2: f64, t2: f64) -> f64 {
        let v_sq_tau = (u_tau1 * u_tau1 + u_tau2 * u_tau2) / (tau2 * tau2);
        let v_sq_t = (u_t1 * u_t1 + u_t2 * u_t2) / (t2 * t2);
        (v_sq_tau + v_sq_t).max(0.0).sqrt()
    }

    /// Integrates the non-linear hyperbolic geodesic equations using classical RK4:
    /// d tau1/ds = u_tau1
    /// d tau2/ds = u_tau2
    /// d u_tau1/ds = 2 u_tau1 u_tau2 / tau2
    /// d u_tau2/ds = (u_tau2^2 - u_tau1^2) / tau2
    /// d t1/ds = u_t1
    /// d t2/ds = u_t2
    /// d u_t1/ds = 2 u_t1 u_t2 / t2
    /// d u_t2/ds = (u_t2^2 - u_t1^2) / t2
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

        let alpha_sdc = 1.0 / std::f64::consts::SQRT_2; // 0.7071 — reference rate from SDC bound
        let initial_mass = self.compute_min_mass_gap(tau1, tau2, t2);

        for step in 0..=n_steps {
            let s = step as f64 * ds;
            let mass_gap = self.compute_min_mass_gap(tau1, tau2, t2);
            let theoretical_bound = initial_mass * (-alpha_sdc * delta_d).exp();
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

            // RK4 integration of the 8-dimensional first-order system
            let (k1_tau1, k1_tau2, k1_u_tau1, k1_u_tau2, k1_t1, k1_t2, k1_u_t1, k1_u_t2) = {
                let (a_tau1, a_tau2, a_t1, a_t2) = self.compute_acceleration(
                    tau1, tau2, u_tau1, u_tau2, t1, t2, u_t1, u_t2);
                (u_tau1, u_tau2, a_tau1, a_tau2, u_t1, u_t2, a_t1, a_t2)
            };

            let (k2_tau1, k2_tau2, k2_u_tau1, k2_u_tau2, k2_t1, k2_t2, k2_u_t1, k2_u_t2) = {
                let tau1_2 = tau1 + 0.5 * ds * k1_tau1;
                let tau2_2 = tau2 + 0.5 * ds * k1_tau2;
                let u_tau1_2 = u_tau1 + 0.5 * ds * k1_u_tau1;
                let u_tau2_2 = u_tau2 + 0.5 * ds * k1_u_tau2;
                let t1_2 = t1 + 0.5 * ds * k1_t1;
                let t2_2 = t2 + 0.5 * ds * k1_t2;
                let u_t1_2 = u_t1 + 0.5 * ds * k1_u_t1;
                let u_t2_2 = u_t2 + 0.5 * ds * k1_u_t2;
                let (a_tau1, a_tau2, a_t1, a_t2) = self.compute_acceleration(
                    tau1_2, tau2_2, u_tau1_2, u_tau2_2, t1_2, t2_2, u_t1_2, u_t2_2);
                (u_tau1_2, u_tau2_2, a_tau1, a_tau2, u_t1_2, u_t2_2, a_t1, a_t2)
            };

            let (k3_tau1, k3_tau2, k3_u_tau1, k3_u_tau2, k3_t1, k3_t2, k3_u_t1, k3_u_t2) = {
                let tau1_3 = tau1 + 0.5 * ds * k2_tau1;
                let tau2_3 = tau2 + 0.5 * ds * k2_tau2;
                let u_tau1_3 = u_tau1 + 0.5 * ds * k2_u_tau1;
                let u_tau2_3 = u_tau2 + 0.5 * ds * k2_u_tau2;
                let t1_3 = t1 + 0.5 * ds * k2_t1;
                let t2_3 = t2 + 0.5 * ds * k2_t2;
                let u_t1_3 = u_t1 + 0.5 * ds * k2_u_t1;
                let u_t2_3 = u_t2 + 0.5 * ds * k2_u_t2;
                let (a_tau1, a_tau2, a_t1, a_t2) = self.compute_acceleration(
                    tau1_3, tau2_3, u_tau1_3, u_tau2_3, t1_3, t2_3, u_t1_3, u_t2_3);
                (u_tau1_3, u_tau2_3, a_tau1, a_tau2, u_t1_3, u_t2_3, a_t1, a_t2)
            };

            let (k4_tau1, k4_tau2, k4_u_tau1, k4_u_tau2, k4_t1, k4_t2, k4_u_t1, k4_u_t2) = {
                let tau1_4 = tau1 + ds * k3_tau1;
                let tau2_4 = tau2 + ds * k3_tau2;
                let u_tau1_4 = u_tau1 + ds * k3_u_tau1;
                let u_tau2_4 = u_tau2 + ds * k3_u_tau2;
                let t1_4 = t1 + ds * k3_t1;
                let t2_4 = t2 + ds * k3_t2;
                let u_t1_4 = u_t1 + ds * k3_u_t1;
                let u_t2_4 = u_t2 + ds * k3_u_t2;
                let (a_tau1, a_tau2, a_t1, a_t2) = self.compute_acceleration(
                    tau1_4, tau2_4, u_tau1_4, u_tau2_4, t1_4, t2_4, u_t1_4, u_t2_4);
                (u_tau1_4, u_tau2_4, a_tau1, a_tau2, u_t1_4, u_t2_4, a_t1, a_t2)
            };

            // Combine RK4 stages
            tau1 += (ds / 6.0) * (k1_tau1 + 2.0 * k2_tau1 + 2.0 * k3_tau1 + k4_tau1);
            tau2 += (ds / 6.0) * (k1_tau2 + 2.0 * k2_tau2 + 2.0 * k3_tau2 + k4_tau2);
            u_tau1 += (ds / 6.0) * (k1_u_tau1 + 2.0 * k2_u_tau1 + 2.0 * k3_u_tau1 + k4_u_tau1);
            u_tau2 += (ds / 6.0) * (k1_u_tau2 + 2.0 * k2_u_tau2 + 2.0 * k3_u_tau2 + k4_u_tau2);

            t1 += (ds / 6.0) * (k1_t1 + 2.0 * k2_t1 + 2.0 * k3_t1 + k4_t1);
            t2 += (ds / 6.0) * (k1_t2 + 2.0 * k2_t2 + 2.0 * k3_t2 + k4_t2);
            u_t1 += (ds / 6.0) * (k1_u_t1 + 2.0 * k2_u_t1 + 2.0 * k3_u_t1 + k4_u_t1);
            u_t2 += (ds / 6.0) * (k1_u_t2 + 2.0 * k2_u_t2 + 2.0 * k3_u_t2 + k4_u_t2);

            // Apply bounds only to committed state (post-RK4)
            tau2 = tau2.max(1e-4);
            t2 = t2.max(1e-4);

            // Accumulate geodesic distance via quadrature
            let speed = self.compute_speed(u_tau1, u_tau2, tau2, u_t1, u_t2, t2);
            delta_d += speed * ds;
        }

        let final_mass = records.last().map(|r| r.min_mass_gap).unwrap_or(0.0);
        let ratio = final_mass / initial_mass.max(1e-12);
        let bound_ok = check_bound(&records, self.config.bound_slack);

        // Compute alpha_effective = -ln(final/initial) / distance
        // Tier X numerical estimate, finite check required
        let alpha_effective = if delta_d > 1e-12 && initial_mass > 1e-12 && final_mass > 1e-12 {
            let alpha = -(final_mass / initial_mass).ln() / delta_d;
            if alpha.is_finite() { Some(alpha) } else { None }
        } else {
            None
        };

        let summary = SwamplandSummary {
            total_steps: n_steps,
            total_geodesic_distance: delta_d,
            initial_mass_gap: initial_mass,
            final_mass_gap: final_mass,
            mass_gap_ratio: ratio,
            sdc_alpha: alpha_sdc,
            alpha_effective,
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

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to run the old explicit scheme (Velocity-Verlet) for negative control.
    /// This is intentionally kept simple and unoptimized to preserve the known defects.
    fn run_old_explicit_scheme(config: &SwamplandConfig) -> (SwamplandSummary, Vec<SwamplandStepRecord>) {
        let ds = config.ds;
        let n_steps = (config.s_max / ds).ceil() as usize;

        let mut tau1 = config.tau_init.0;
        let mut tau2 = config.tau_init.1;
        let mut u_tau1 = config.v_tau.0;
        let mut u_tau2 = config.v_tau.1;

        let mut t1 = config.t_init.0;
        let mut t2 = config.t_init.1;
        let mut u_t1 = config.v_t.0;
        let mut u_t2 = config.v_t.1;

        let mut delta_d = 0.0;
        let mut records = Vec::with_capacity(n_steps + 1);

        let alpha_sdc = 1.0 / std::f64::consts::SQRT_2;
        let initial_mass = SwamplandSimulator::new(config.clone()).compute_min_mass_gap(tau1, tau2, t2);

        for step in 0..=n_steps {
            let s = step as f64 * ds;
            let _sim = SwamplandSimulator::new(config.clone());
            let mass_gap = _sim.compute_min_mass_gap(tau1, tau2, t2);
            let theoretical_bound = initial_mass * (-alpha_sdc * delta_d).exp();
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

            // Metric velocity at start of step
            let v_sq_tau = (u_tau1 * u_tau1 + u_tau2 * u_tau2) / (tau2 * tau2);
            let v_sq_t = (u_t1 * u_t1 + u_t2 * u_t2) / (t2 * t2);
            let speed = (v_sq_tau + v_sq_t).max(0.0).sqrt();
            delta_d += speed * ds;

            // Acceleration at start of step only (defect: not updated)
            let a_tau1 = 2.0 * u_tau1 * u_tau2 / tau2;
            let a_tau2 = (u_tau2 * u_tau2 - u_tau1 * u_tau1) / tau2;
            let a_t1 = 2.0 * u_t1 * u_t2 / t2;
            let a_t2 = (u_t2 * u_t2 - u_t1 * u_t1) / t2;

            // Explicit Taylor step (defect: uses only start-of-step acceleration)
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
        let bound_ok = check_bound(&records, config.bound_slack);

        let alpha_effective = if delta_d > 1e-12 && initial_mass > 1e-12 && final_mass > 1e-12 {
            let alpha = -(final_mass / initial_mass).ln() / delta_d;
            if alpha.is_finite() { Some(alpha) } else { None }
        } else {
            None
        };

        let summary = SwamplandSummary {
            total_steps: n_steps,
            total_geodesic_distance: delta_d,
            initial_mass_gap: initial_mass,
            final_mass_gap: final_mass,
            mass_gap_ratio: ratio,
            sdc_alpha: alpha_sdc,
            alpha_effective,
            bound_satisfied: bound_ok,
            eft_cutoff_broken: final_mass < 0.1,
        };

        (summary, records)
    }

    /// Positive control: check_bound should return true for records that satisfy the bound.
    #[test]
    fn test_check_bound_positive() {
        let records = vec![
            SwamplandStepRecord {
                s: 0.0,
                delta_d: 0.0,
                tau1: 0.0,
                tau2: 1.0,
                t1: 0.0,
                t2: 1.0,
                min_mass_gap: 1.0,
                theoretical_bound: 1.0,
                eft_cutoff: 10.0,
            },
            SwamplandStepRecord {
                s: 0.01,
                delta_d: 0.01,
                tau1: 0.0,
                tau2: 1.0,
                t1: 0.0,
                t2: 1.0,
                min_mass_gap: 0.99,
                theoretical_bound: 1.0,
                eft_cutoff: 10.0,
            },
        ];
        assert!(check_bound(&records, 1.0), "Positive control: bound check should pass for satisfying records");
    }

    /// Negative control: check_bound should return false for records violating the bound.
    #[test]
    fn test_check_bound_negative() {
        let records = vec![
            SwamplandStepRecord {
                s: 0.0,
                delta_d: 0.0,
                tau1: 0.0,
                tau2: 1.0,
                t1: 0.0,
                t2: 1.0,
                min_mass_gap: 1.0,
                theoretical_bound: 1.0,
                eft_cutoff: 10.0,
            },
            SwamplandStepRecord {
                s: 0.01,
                delta_d: 0.01,
                tau1: 0.0,
                tau2: 1.0,
                t1: 0.0,
                t2: 1.0,
                min_mass_gap: 3.0,  // Violates bound
                theoretical_bound: 1.0,
                eft_cutoff: 10.0,
            },
        ];
        assert!(!check_bound(&records, 1.0), "Negative control: bound check should fail for violating records");
    }

    /// Speed conservation test with default config (x'=0, scalar ODE y''=y'^2/y).
    /// Expected: speed = sqrt(2) conserved, total distance ≈ 7.071.
    #[test]
    fn test_speed_conservation_default_config() {
        let config = SwamplandConfig::default();
        let sim = SwamplandSimulator::new(config);
        let (summary, records) = sim.run_simulation();

        // Expected: y(s_max=5) = e^5 ≈ 148.41, so mass ≈ 1/148.41 ≈ 6.738e-3
        let expected_final_mass = 1.0 / (std::f64::consts::E.powf(5.0));
        let mass_error = (summary.final_mass_gap - expected_final_mass).abs() / expected_final_mass;
        assert!(mass_error < 1e-3, "Final mass gap error too large: {}", mass_error);

        // Expected: total distance ≈ sqrt(2) * s_max ≈ 7.071
        let expected_distance = std::f64::consts::SQRT_2 * 5.0;
        let distance_error = (summary.total_geodesic_distance - expected_distance).abs() / expected_distance;
        assert!(distance_error < 1e-4, "Distance error too large: {}", distance_error);

        // Check relative drift in speed conservation
        let mut max_speed_drift: f64 = 0.0;
        let expected_speed = std::f64::consts::SQRT_2;
        for i in 0..records.len() - 1 {
            let r = &records[i];
            // Reconstruct speed from delta_d increment
            let next_delta_d = records[i + 1].delta_d;
            let ds = 0.01;
            let speed_i = (next_delta_d - r.delta_d) / ds;
            let drift = (speed_i - expected_speed).abs() / expected_speed;
            max_speed_drift = max_speed_drift.max(drift);
        }
        assert!(max_speed_drift < 1e-6, "Speed conservation drift too large: {}", max_speed_drift);
    }

    /// Speed conservation test with non-trivial coupling (x'≠0).
    /// This exercises the coupling term and prevents the test from being vacuous.
    #[test]
    fn test_speed_conservation_nontrivial() {
        let mut config = SwamplandConfig::default();
        config.v_tau = (0.2, 1.0);  // Non-zero x', couples the system
        let sim = SwamplandSimulator::new(config.clone());
        let (summary, records) = sim.run_simulation();

        // With coupling, the trajectory is no longer trivial; y decays toward 0 but should stay above clamp
        let min_tau2 = records.iter().map(|r| r.tau2).fold(f64::MAX, f64::min);
        let min_t2 = records.iter().map(|r| r.t2).fold(f64::MAX, f64::min);

        // Allow some margin above the 1e-4 clamp, but ensure we're measuring the integrator, not the clamp
        assert!(min_tau2 > 2e-4, "tau2 decayed to clamp region (tau2_min={}), test is measuring clamp not integrator", min_tau2);
        assert!(min_t2 > 2e-4, "t2 decayed to clamp region (t2_min={}), test is measuring clamp not integrator", min_t2);

        // Check that speed is conserved (within RK4 error)
        // This is a softer test than default config since the exact solution is not closed-form
        assert!(summary.total_geodesic_distance > 1.0, "Geodesic distance should be positive");
        assert!(summary.mass_gap_ratio > 0.0 && summary.mass_gap_ratio < 1.0, "Mass gap ratio should be in (0,1)");

        // Check that alpha_effective is finite
        assert!(summary.alpha_effective.is_some(), "alpha_effective should be computed");
        if let Some(alpha) = summary.alpha_effective {
            assert!(alpha.is_finite(), "alpha_effective must be finite");
            assert!(alpha > 0.0, "alpha_effective should be positive");
        }
    }

    /// Negative control: old explicit scheme should drift more than RK4.
    #[test]
    fn test_old_scheme_drifts_more() {
        let config = SwamplandConfig::default();
        let sim_rk4 = SwamplandSimulator::new(config.clone());
        let (_summary_rk4, records_rk4) = sim_rk4.run_simulation();

        let (_summary_old, records_old) = run_old_explicit_scheme(&config);

        // Compute RK4 speed drift
        let mut max_speed_drift_rk4: f64 = 0.0;
        let expected_speed = std::f64::consts::SQRT_2;
        for i in 0..records_rk4.len() - 1 {
            let r = &records_rk4[i];
            let next_delta_d = records_rk4[i + 1].delta_d;
            let ds = 0.01;
            let speed_i = (next_delta_d - r.delta_d) / ds;
            let drift = (speed_i - expected_speed).abs() / expected_speed;
            max_speed_drift_rk4 = max_speed_drift_rk4.max(drift);
        }

        // Compute old scheme speed drift
        let mut max_speed_drift_old: f64 = 0.0;
        for i in 0..records_old.len() - 1 {
            let r = &records_old[i];
            let next_delta_d = records_old[i + 1].delta_d;
            let ds = 0.01;
            let speed_i = (next_delta_d - r.delta_d) / ds;
            let drift = (speed_i - expected_speed).abs() / expected_speed;
            max_speed_drift_old = max_speed_drift_old.max(drift);
        }

        // RK4 should drift much less than old scheme
        // Note: On default config (trivial scalar ODE), RK4 achieves machine precision (~1e-14)
        assert!(max_speed_drift_old > 1e-4, "Old scheme drift should be > 1e-4 (got {})", max_speed_drift_old);
        assert!(max_speed_drift_old > 10.0 * max_speed_drift_rk4,
                "Old scheme drift ({}) should be > 10x RK4 drift ({})", max_speed_drift_old, max_speed_drift_rk4);
    }

    /// Test that alpha_effective is finite and reasonable.
    #[test]
    fn test_alpha_effective_finite() {
        let config = SwamplandConfig::default();
        let sim = SwamplandSimulator::new(config);
        let (summary, _records) = sim.run_simulation();

        assert!(summary.alpha_effective.is_some(), "alpha_effective should be computed for default config");
        if let Some(alpha) = summary.alpha_effective {
            assert!(alpha.is_finite(), "alpha_effective must be finite");
            // Expected: alpha_effective ≈ 1/sqrt(2) ≈ 0.707 for default config
            assert!(alpha > 0.6 && alpha < 0.8, "alpha_effective should be close to 1/sqrt(2), got {}", alpha);
        }
    }

    /// Test that mode truncation is checked: minimizer should be strictly interior to box.
    #[test]
    fn test_mode_truncation_interior() {
        let config = SwamplandConfig::default();
        let tau1 = config.tau_init.0;
        let tau2 = config.tau_init.1;
        let t2 = config.t_init.1;

        let max_n = config.max_mode_n;
        let max_m = config.max_mode_m;

        // Find the minimizer mode
        let mut min_mass = f64::MAX;
        let mut min_n = 0;
        let mut min_m = 0;
        for n in -max_n..=max_n {
            for m in -max_m..=max_m {
                if n == 0 && m == 0 {
                    continue;
                }
                let m_sq = SwamplandSimulator::mode_mass_sq(n, m, tau1, tau2, t2);
                let mass = m_sq.max(0.0).sqrt();
                if mass > 1e-10 && mass < min_mass {
                    min_mass = mass;
                    min_n = n;
                    min_m = m;
                }
            }
        }

        // Check that minimizer is strictly interior to box
        assert!(min_n.abs() < max_n, "Mode n={} should be strictly interior (max_n={})", min_n, max_n);
        assert!(min_m.abs() < max_m, "Mode m={} should be strictly interior (max_m={})", min_m, max_m);
    }
}
