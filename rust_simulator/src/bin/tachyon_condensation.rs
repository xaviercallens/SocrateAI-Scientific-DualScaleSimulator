use rust_simulator::tachyon_condensation::{TachyonSimulator, TachyonConfig};

fn main() {
    let mut config = TachyonConfig::default();
    let args: Vec<String> = std::env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--t-max" => {
                if i + 1 < args.len() {
                    config.t_max = args[i + 1].parse().unwrap_or(20.0);
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

    let mut sim = TachyonSimulator::new(config);
    let (summary, records) = sim.run_simulation();
    sim.export_results(&summary, &records).expect("Failed to export tachyon telemetry");
    println!("=== Tachyon Condensation Simulation Success ===");
}
