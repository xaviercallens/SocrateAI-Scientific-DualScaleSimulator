use rust_simulator::kummer_langevin::{KummerLangevinSimulator, SimConfig as KummerConfig};
use rust_simulator::swampland_geodesic::{SwamplandSimulator, SwamplandConfig};
use rust_simulator::tachyon_condensation::{TachyonSimulator, TachyonConfig};
use rust_simulator::vacuum_decay_cdl::{VacuumDecaySimulator, VacuumDecayConfig};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut mode = "all".to_string();

    // The CLI flags below used to be PARSED NOWHERE: main() read only --mode and
    // then built every config with Config::default(), so the --grid-size / --t-max /
    // --output-csv / --output-json that workshopcosmo.py passes were silently
    // ignored and the documented invocation did not control the simulation.
    // Recorded as S1-F11 in audit/STREAM1_BRIDGE.md.
    let mut grid_size: Option<usize> = None;
    let mut t_max: Option<f64> = None;
    let mut seed: Option<u64> = None;
    let mut output_csv: Option<String> = None;
    let mut output_json: Option<String> = None;

    let mut i = 1;
    while i < args.len() {
        let next = |i: usize| -> Option<String> { args.get(i + 1).cloned() };
        match args[i].as_str() {
            "--mode" => { if let Some(v) = next(i) { mode = v; i += 1; } }
            "--grid-size" => { if let Some(v) = next(i) { grid_size = v.parse().ok(); i += 1; } }
            "--t-max" => { if let Some(v) = next(i) { t_max = v.parse().ok(); i += 1; } }
            "--seed" => { if let Some(v) = next(i) { seed = v.parse().ok(); i += 1; } }
            "--output-csv" => { if let Some(v) = next(i) { output_csv = Some(v); i += 1; } }
            "--output-json" => { if let Some(v) = next(i) { output_json = Some(v); i += 1; } }
            other => {
                if other.starts_with("--") {
                    eprintln!("warning: unrecognised flag {} ignored", other);
                }
            }
        }
        i += 1;
    }

    println!("============================================================");
    println!("  LeanFlow / rusty-SUNDIALS Multi-Physics String Engine");
    println!("  Mode: {}", mode);
    println!("============================================================");

    if mode == "kummer" || mode == "all" {
        println!("\n>>> Running Loop 0: Kummer Orbifold Langevin Dynamics...");
        let mut config = KummerConfig::default();
        if let Some(v) = grid_size { config.grid_size = v; }
        if let Some(v) = t_max { config.t_max = v; }
        if let Some(v) = seed { config.seed = v; }
        if let Some(ref v) = output_csv { config.output_csv = v.clone(); }
        if let Some(ref v) = output_json { config.output_json = v.clone(); }
        println!("    config: grid_size={} t_max={} seed={}", config.grid_size, config.t_max, config.seed);
        let mut sim = KummerLangevinSimulator::new(config);
        let (summary, records) = sim.run_simulation();
        sim.export_results(&summary, &records).expect("Failed to export Kummer telemetry");
        println!(">>> Kummer Simulation Complete. Wall pixels: {}, Strings: {}", summary.final_wall_pixel_count, summary.final_string_count);
    }

    if mode == "swampland" || mode == "all" {
        println!("\n>>> Running Loop 1: Swampland Distance Conjecture Geodesic Flow...");
        let mut config = SwamplandConfig::default();
        if let Some(v) = t_max { config.s_max = v; }
        if let Some(ref v) = output_csv { config.output_csv = v.clone(); }
        if let Some(ref v) = output_json { config.output_json = v.clone(); }
        println!("    config: s_max={}", config.s_max);
        let sim = SwamplandSimulator::new(config);
        let (summary, records) = sim.run_simulation();
        sim.export_results(&summary, &records).expect("Failed to export Swampland telemetry");
        println!(">>> Swampland Complete. Geodesic distance: {:.4}, Mass gap collapse: {:.4e}", summary.total_geodesic_distance, summary.mass_gap_ratio);
    }

    if mode == "tachyon" || mode == "all" {
        println!("\n>>> Running Loop 2: Tachyon Condensation & Sen Soliton Formation...");
        let mut config = TachyonConfig::default();
        if let Some(v) = t_max { config.t_max = v; }
        if let Some(v) = seed { config.seed = v; }
        if let Some(ref v) = output_csv { config.output_csv = v.clone(); }
        if let Some(ref v) = output_json { config.output_json = v.clone(); }
        println!("    config: t_max={} seed={}", config.t_max, config.seed);
        let mut sim = TachyonSimulator::new(config);
        let (summary, records) = sim.run_simulation();
        sim.export_results(&summary, &records).expect("Failed to export Tachyon telemetry");
        println!(">>> Tachyon Complete. Soliton peak energy: {:.4}, K-theory conserved: {}", summary.soliton_peak_energy, summary.k_theory_charge_conserved);
    }

    if mode == "vacuum-decay" || mode == "all" {
        println!("\n>>> Running Loop 3: Flat-Space Coleman Bounce & Holographic c-Theorem...");
        let mut config = VacuumDecayConfig::default();
        if let Some(ref v) = output_csv { config.output_csv = v.clone(); }
        if let Some(ref v) = output_json { config.output_json = v.clone(); }
        let sim = VacuumDecaySimulator::new(config);
        match sim.run_simulation() {
            Ok((summary, records)) => {
                sim.export_results(&summary, &records).expect("Failed to export Vacuum Decay telemetry");
                println!(">>> Vacuum Decay Complete. Bounce action S_E: {:.6}, Δc: {} (c-theorem: {})", summary.euclidean_bounce_action, summary.delta_c, summary.c_theorem_satisfied);
            }
            Err(e) => eprintln!("Vacuum Decay simulation failed: {}", e),
        }
    }

    println!("\nAll requested simulations completed successfully.");
}
