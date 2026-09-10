use rust_simulator::kummer_langevin::{KummerLangevinSimulator, SimConfig as KummerConfig};
use rust_simulator::swampland_geodesic::{SwamplandSimulator, SwamplandConfig};
use rust_simulator::tachyon_condensation::{TachyonSimulator, TachyonConfig};
use rust_simulator::vacuum_decay_cdl::{VacuumDecaySimulator, VacuumDecayConfig};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut mode = "all".to_string();

    let mut i = 1;
    while i < args.len() {
        if args[i] == "--mode" && i + 1 < args.len() {
            mode = args[i + 1].clone();
            i += 1;
        }
        i += 1;
    }

    println!("============================================================");
    println!("  LeanFlow / rusty-SUNDIALS Multi-Physics String Engine");
    println!("  Mode: {}", mode);
    println!("============================================================");

    if mode == "kummer" || mode == "all" {
        println!("\n>>> Running Loop 0: Kummer Orbifold Langevin Dynamics...");
        let config = KummerConfig::default();
        let mut sim = KummerLangevinSimulator::new(config);
        let (summary, records) = sim.run_simulation();
        sim.export_results(&summary, &records).expect("Failed to export Kummer telemetry");
        println!(">>> Kummer Simulation Complete. Wall pixels: {}, Strings: {}", summary.final_wall_pixel_count, summary.final_string_count);
    }

    if mode == "swampland" || mode == "all" {
        println!("\n>>> Running Loop 1: Swampland Distance Conjecture Geodesic Flow...");
        let config = SwamplandConfig::default();
        let sim = SwamplandSimulator::new(config);
        let (summary, records) = sim.run_simulation();
        sim.export_results(&summary, &records).expect("Failed to export Swampland telemetry");
        println!(">>> Swampland Complete. Geodesic distance: {:.4}, Mass gap collapse: {:.4e}", summary.total_geodesic_distance, summary.mass_gap_ratio);
    }

    if mode == "tachyon" || mode == "all" {
        println!("\n>>> Running Loop 2: Tachyon Condensation & Sen Soliton Formation...");
        let config = TachyonConfig::default();
        let mut sim = TachyonSimulator::new(config);
        let (summary, records) = sim.run_simulation();
        sim.export_results(&summary, &records).expect("Failed to export Tachyon telemetry");
        println!(">>> Tachyon Complete. Soliton peak energy: {:.4}, K-theory conserved: {}", summary.soliton_peak_energy, summary.k_theory_charge_conserved);
    }

    if mode == "vacuum-decay" || mode == "all" {
        println!("\n>>> Running Loop 3: Coleman-De Luccia Vacuum Decay & Holographic c-Theorem...");
        let config = VacuumDecayConfig::default();
        let sim = VacuumDecaySimulator::new(config);
        let (summary, records) = sim.run_simulation();
        sim.export_results(&summary, &records).expect("Failed to export Vacuum Decay telemetry");
        println!(">>> Vacuum Decay Complete. Bounce action S_E: {:.4}, Δc: {} (c-theorem: {})", summary.euclidean_bounce_action, summary.delta_c, summary.c_theorem_satisfied);
    }

    println!("\nAll requested simulations completed successfully.");
}
