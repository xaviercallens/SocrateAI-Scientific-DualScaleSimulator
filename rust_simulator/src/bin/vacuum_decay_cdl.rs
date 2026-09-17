use rust_simulator::vacuum_decay_cdl::{VacuumDecaySimulator, VacuumDecayConfig};

fn main() {
    let mut config = VacuumDecayConfig::default();
    let args: Vec<String> = std::env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--rho-max" => {
                if i + 1 < args.len() {
                    config.rho_max = args[i + 1].parse().unwrap_or(12.0);
                    i += 1;
                }
            }
            "--output-csv" => {
                if i + 1 < args.len() {
                    config.output_csv = args[i + 1].clone();
                    i += 1;
                }
            }
            "--output-json" => {
                if i + 1 < args.len() {
                    config.output_json = args[i + 1].clone();
                    i += 1;
                }
            }
            _ => {}
        }
        i += 1;
    }

    let sim = VacuumDecaySimulator::new(config);
    match sim.run_simulation() {
        Ok((summary, records)) => {
            sim.export_results(&summary, &records).expect("Failed to export vacuum decay telemetry");
            println!("=== Flat-Space Coleman Bounce Simulation Success ===");
        }
        Err(e) => {
            eprintln!("Simulation failed: {}", e);
            std::process::exit(1);
        }
    }
}
