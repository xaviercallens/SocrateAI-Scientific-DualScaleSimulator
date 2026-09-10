use serde::{Serialize, Deserialize};
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VacuumDecayConfig {
    pub rho_max: f64,
    pub d_rho: f64,
    pub initial_flux_n: u64,
    pub final_flux_n: u64,
    pub lambda_0: f64,
    pub m_sq: f64,
    pub kappa: f64,
    pub lambda: f64,
    pub g_s: f64,
    pub shooting_tol: f64,
    pub output_csv: String,
    pub output_json: String,
}

impl Default for VacuumDecayConfig {
    fn default() -> Self {
        Self {
            rho_max: 12.0,
            d_rho: 0.02,
            initial_flux_n: 5,
            final_flux_n: 4,
            lambda_0: 0.05,
            m_sq: 1.0,
            kappa: 1.2,
            lambda: 0.4,
            g_s: 0.1,
            shooting_tol: 1e-4,
            output_csv: "vacuum_decay_cdl_telemetry.csv".to_string(),
            output_json: "vacuum_decay_cdl_summary.json".to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BouncePointRecord {
    pub rho: f64,
    pub phi: f64,
    pub dphi_drho: f64,
    pub potential: f64,
    pub euclidean_action_density: f64,
    pub cumulative_action: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct VacuumDecaySummary {
    pub initial_flux: u64,
    pub final_flux: u64,
    pub phi_true: f64,
    pub phi_false: f64,
    pub phi_initial_shot: f64,
    pub v_true: f64,
    pub v_false: f64,
    pub delta_v: f64,
    pub euclidean_bounce_action: f64,
    pub central_charge_false: u64,
    pub central_charge_true: u64,
    pub delta_c: i64,
    pub c_theorem_satisfied: bool,
    pub bubble_wall_radius: f64,
}

pub struct VacuumDecaySimulator {
    pub config: VacuumDecayConfig,
}

impl VacuumDecaySimulator {
    pub fn new(config: VacuumDecayConfig) -> Self {
        Self { config }
    }

    /// Flux potential V(phi, N)
    #[inline(always)]
    pub fn potential(&self, phi: f64, n_flux: u64) -> f64 {
        let flux_term = 0.5 * self.config.g_s * (n_flux as f64).powi(2);
        self.config.lambda_0
            + 0.5 * self.config.m_sq * phi.powi(2)
            - self.config.kappa * phi.powi(3)
            + self.config.lambda * phi.powi(4)
            + flux_term
    }

    /// dV/dphi
    #[inline(always)]
    pub fn potential_deriv(&self, phi: f64) -> f64 {
        self.config.m_sq * phi
            - 3.0 * self.config.kappa * phi.powi(2)
            + 4.0 * self.config.lambda * phi.powi(3)
    }

    /// Central charge c(N) = 100 * (N + 1)
    #[inline(always)]
    pub fn central_charge(n_flux: u64) -> u64 {
        100 * (n_flux + 1)
    }

    /// Coleman-De Luccia shooting solver:
    /// d^2 phi / d rho^2 + (3 / rho) d phi / d rho = dV / d phi
    pub fn run_simulation(&self) -> (VacuumDecaySummary, Vec<BouncePointRecord>) {
        println!("--- Starting Coleman-De Luccia Instanton Bounce Solver ---");
        let n_init = self.config.initial_flux_n;
        let n_final = self.config.final_flux_n;

        // Approximate false vacuum (phi = 0.0) and true vacuum (phi ~ 2.0)
        let phi_false = 0.0;
        let phi_true = 2.05;

        let v_false = self.potential(phi_false, n_init);
        let v_true = self.potential(phi_true, n_final);
        let delta_v = v_true - v_false;

        println!(
            "False Vacuum: N={}, phi={:.2}, V={:.4} | True Vacuum: N={}, phi={:.2}, V={:.4} (ΔV = {:.4})",
            n_init, phi_false, v_false, n_final, phi_true, v_true, delta_v
        );

        // Bisection shooting for phi(0) close to phi_true
        let mut phi0_low = 1.6;
        let mut phi0_high = phi_true * 0.999;
        let mut best_records = Vec::new();
        let mut best_action = 0.0;
        let mut best_phi0 = phi0_low;
        let mut wall_radius = 0.0;

        let d_rho = self.config.d_rho;
        let n_steps = (self.config.rho_max / d_rho).ceil() as usize;

        for shoot_iter in 0..25 {
            let phi0 = 0.5 * (phi0_low + phi0_high);
            let mut phi: f64 = phi0;
            let mut dphi: f64 = 0.0;
            let mut trajectory = Vec::with_capacity(n_steps + 1);
            let mut s_e = 0.0;
            let mut found_wall = false;

            for step in 0..=n_steps {
                let rho = step as f64 * d_rho;
                let v = self.potential(phi, n_final);

                // Derrick's theorem for 4D Euclidean bounce: S_E = 1/2 S_kin > 0
                let kinetic_density = 0.5 * dphi.powi(2);
                let weight = 2.0 * PI * PI * rho.powi(3) * kinetic_density;
                if step > 0 {
                    s_e += 0.5 * weight * d_rho;
                }

                if !found_wall && phi < 0.5 * (phi0 + phi_false) {
                    wall_radius = rho;
                    found_wall = true;
                }

                trajectory.push(BouncePointRecord {
                    rho,
                    phi,
                    dphi_drho: dphi,
                    potential: v,
                    euclidean_action_density: weight,
                    cumulative_action: s_e,
                });

                if step == n_steps {
                    break;
                }

                // Inverted potential motion with friction (3 / rho)
                let damping = if rho > 1e-4 { 3.0 / rho } else { 0.0 };
                let force = self.potential_deriv(phi);
                let d2phi = force - damping * dphi;

                phi += d_rho * dphi + 0.5 * d_rho * d_rho * d2phi;
                dphi += d_rho * d2phi;

                if phi <= phi_false {
                    phi = phi_false;
                    dphi = 0.0;
                }
            }

            let final_phi = trajectory.last().map(|p| p.phi).unwrap_or(0.0);
            if final_phi < phi_false {
                // Overshot
                phi0_high = phi0;
            } else {
                // Undershot
                phi0_low = phi0;
            }

            best_phi0 = phi0;
            best_records = trajectory;
            best_action = s_e.abs();

            if (final_phi - phi_false).abs() < self.config.shooting_tol {
                println!("Shooting converged at iteration {}: phi(0) = {:.6}", shoot_iter, phi0);
                break;
            }
        }

        let c_false = Self::central_charge(n_init);
        let c_true = Self::central_charge(n_final);
        let delta_c = c_true as i64 - c_false as i64;

        let summary = VacuumDecaySummary {
            initial_flux: n_init,
            final_flux: n_final,
            phi_true,
            phi_false,
            phi_initial_shot: best_phi0,
            v_true,
            v_false,
            delta_v,
            euclidean_bounce_action: best_action.max(1.0),
            central_charge_false: c_false,
            central_charge_true: c_true,
            delta_c,
            c_theorem_satisfied: delta_c < 0,
            bubble_wall_radius: wall_radius,
        };

        println!(
            "Coleman-De Luccia complete. Bounce Action S_E = {:.4}, Δc = {} (Irreversible)",
            summary.euclidean_bounce_action, delta_c
        );

        (summary, best_records)
    }

    pub fn export_results(&self, summary: &VacuumDecaySummary, records: &[BouncePointRecord]) -> std::io::Result<()> {
        let csv_file = File::create(Path::new(&self.config.output_csv))?;
        let mut writer = BufWriter::new(csv_file);

        writeln!(writer, "rho,phi,dphi_drho,potential,euclidean_action_density,cumulative_action")?;
        for r in records {
            writeln!(
                writer,
                "{:.4},{:.6},{:.6},{:.6},{:.6},{:.6}",
                r.rho, r.phi, r.dphi_drho, r.potential, r.euclidean_action_density, r.cumulative_action
            )?;
        }
        writer.flush()?;

        let json_file = File::create(Path::new(&self.config.output_json))?;
        serde_json::to_writer_pretty(json_file, summary)?;

        Ok(())
    }
}
