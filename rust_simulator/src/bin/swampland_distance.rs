use rust_simulator::swampland_geodesic::{SwamplandSimulator, SwamplandConfig};

fn main() {
    let mut config = SwamplandConfig::default();
    let args: Vec<String> = std::env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--s-max" => {
                if i + 1 < args.len() {
                    config.s_max = args[i + 1].parse().unwrap_or(5.0);
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

    let sim = SwamplandSimulator::new(config);
    let (summary, records) = sim.run_simulation();
    sim.export_results(&summary, &records).expect("Failed to export swampland telemetry");
    println!("=== Swampland Distance Simulation Success ===");
}
