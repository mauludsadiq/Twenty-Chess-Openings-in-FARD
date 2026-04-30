use anyhow::{anyhow, Context, Result};
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::{env, fs};

#[derive(Debug, Deserialize, Clone)]
struct Edge {
    opening_id: String,
    endgame_id: String,
    #[serde(default)]
    motif_id: Option<String>,
    #[serde(default)]
    p_hat: Option<f64>,
    #[serde(default)]
    support_games: Option<u64>,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum EdgeFile {
    Wrapped { edges: Vec<Edge> },
    Raw(Vec<Edge>),
}

#[derive(Debug, Deserialize)]
struct OpeningsFile {
    openings: Vec<Opening>,
}

#[derive(Debug, Deserialize)]
struct Opening {
    id: String,
    name: String,
    #[serde(default)]
    eco_codes: Vec<String>,
    #[serde(default)]
    tags: Vec<String>,
    #[serde(default)]
    notes: String,
}

#[derive(Debug, Deserialize)]
struct EndgamesFile {
    endgames: Vec<Endgame>,
}

#[derive(Debug, Deserialize)]
struct Endgame {
    id: String,
    name: String,
    #[serde(default)]
    tags: Vec<String>,
    #[serde(default)]
    notes: String,
}

fn read_json<T: for<'de> Deserialize<'de>>(path: &str) -> Result<T> {
    let text = fs::read_to_string(path).with_context(|| format!("failed to read {path}"))?;
    serde_json::from_str(&text).with_context(|| format!("failed to parse {path}"))
}

fn load_edges(path: &str) -> Result<Vec<Edge>> {
    let parsed: EdgeFile = read_json(path)?;
    Ok(match parsed {
        EdgeFile::Wrapped { edges } => edges,
        EdgeFile::Raw(edges) => edges,
    })
}

fn shannon_entropy(probs: &[f64]) -> f64 {
    probs.iter().copied().filter(|p| *p > 0.0).map(|p| -p * p.log2()).sum()
}

fn usable_edges() -> Result<Vec<Edge>> {
    Ok(load_edges("data/Z_edges_fitted.json")?
        .into_iter()
        .filter(|e| e.p_hat.unwrap_or(0.0) > 0.0)
        .collect())
}

fn cmd_stats() -> Result<()> {
    let usable = usable_edges()?;
    if usable.is_empty() {
        println!("No usable edges in data/Z_edges_fitted.json");
        return Ok(());
    }

    let openings: Vec<String> = usable.iter().map(|e| e.opening_id.clone()).collect::<BTreeSet<_>>().into_iter().collect();
    let endgames: Vec<String> = usable.iter().map(|e| e.endgame_id.clone()).collect::<BTreeSet<_>>().into_iter().collect();

    let mut z_states = BTreeSet::new();
    for e in &usable {
        z_states.insert((e.opening_id.clone(), e.endgame_id.clone(), e.motif_id.clone().unwrap_or_else(|| "gamma_0".to_string())));
    }
    let z_states: Vec<(String, String, String)> = z_states.into_iter().collect();

    let opening_index: BTreeMap<String, usize> = openings.iter().enumerate().map(|(i, o)| (o.clone(), i)).collect();
    let z_index: BTreeMap<(String, String, String), usize> = z_states.iter().enumerate().map(|(i, z)| (z.clone(), i)).collect();

    let mut p = vec![vec![0.0_f64; z_states.len()]; openings.len()];
    let mut support_counts = vec![0_u64; openings.len()];

    for e in &usable {
        let i = opening_index[&e.opening_id];
        let z = (e.opening_id.clone(), e.endgame_id.clone(), e.motif_id.clone().unwrap_or_else(|| "gamma_0".to_string()));
        let m = z_index[&z];
        p[i][m] += e.p_hat.unwrap_or(0.0);
        support_counts[i] += e.support_games.unwrap_or(0);
    }

    for row in &mut p {
        let s: f64 = row.iter().sum();
        if s > 0.0 {
            for v in row {
                *v /= s;
            }
        }
    }

    let total_support: u64 = support_counts.iter().sum();
    if total_support == 0 {
        println!("Total support_games is zero; cannot build opening prior.");
        return Ok(());
    }

    let pi_openings: Vec<f64> = support_counts.iter().map(|c| *c as f64 / total_support as f64).collect();
    let h_o = shannon_entropy(&pi_openings);

    let mut mu_z = vec![0.0_f64; z_states.len()];
    for (i, pi_i) in pi_openings.iter().copied().enumerate() {
        for (m, p_im) in p[i].iter().copied().enumerate() {
            mu_z[m] += pi_i * p_im;
        }
    }

    let h_z = shannon_entropy(&mu_z);
    let endgame_index: BTreeMap<String, usize> = endgames.iter().enumerate().map(|(i, e)| (e.clone(), i)).collect();

    let mut mu_e = vec![0.0_f64; endgames.len()];
    for (m, p_m) in mu_z.iter().copied().enumerate() {
        mu_e[endgame_index[&z_states[m].1]] += p_m;
    }

    let h_e = shannon_entropy(&mu_e);

    println!("Chess opening→motif→endgame stats");
    println!("----------------------------------");
    println!("Openings: {}", openings.len());
    println!("Z states (opening,endgame,motif): {}", z_states.len());
    println!("Endgames: {}", endgames.len());
    println!("Total support_games: {}", total_support);
    println!();
    println!("H(O)               = {:.6} bits", h_o);
    println!("H(Z)               = {:.6} bits", h_z);
    println!("H(E)               = {:.6} bits", h_e);
    println!();
    println!("O → Z emanation (adding motif/endgame detail):");
    println!("  ΔH_emanation     = H(Z) - H(O) = {:.6} bits", h_z - h_o);
    println!("  ΔH_emanation/H(O)= {:.6}", (h_z - h_o) / h_o);
    println!();
    println!("Z → E collapse (forgetting motif, keeping only endgame):");
    println!("  ΔH_ZE            = H(Z) - H(E) = {:.6} bits", h_z - h_e);
    println!("  ΔH_ZE/H(Z)       = {:.6}", (h_z - h_e) / h_z);
    println!();

    let mut ranked: Vec<(&String, f64)> = endgames.iter().zip(mu_e.iter().copied()).collect();
    ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    println!("Endgame distribution under emanated prior:");
    for (e_id, prob) in ranked {
        println!("  {}: {:.6}", e_id, prob);
    }

    Ok(())
}

fn cmd_kernel() -> Result<()> {
    let usable = usable_edges()?;
    let endgames: Vec<String> = usable.iter().map(|e| e.endgame_id.clone()).collect::<BTreeSet<_>>().into_iter().collect();

    let mut rows: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    for e in usable {
        *rows.entry(e.opening_id).or_default().entry(e.endgame_id).or_default() += e.p_hat.unwrap_or(0.0);
    }

    println!("Opening -> Closing (endgame) kernel");
    println!("-----------------------------------");
    println!("Endgames: {}", endgames.join(", "));
    println!();

    for (opening, mut row) in rows {
        let sum: f64 = row.values().sum();
        if sum > 0.0 {
            for v in row.values_mut() {
                *v /= sum;
            }
        }
        let mut ranked: Vec<(String, f64)> = row.into_iter().collect();
        ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        println!("{opening}:");
        for (endgame, p) in ranked {
            if p > 0.0 {
                println!("  {endgame}: {p:.6}");
            }
        }
        println!();
    }

    Ok(())
}

fn cmd_catalog() -> Result<()> {
    let openings: OpeningsFile = read_json("data/openings.json")?;
    let endgames: EndgamesFile = read_json("data/endgames.json")?;

    println!("Opening catalog");
    println!("---------------");
    for o in openings.openings {
        println!("{} | {} | ECO: {} | tags: {}", o.id, o.name, o.eco_codes.join(","), o.tags.join(","));
        if !o.notes.is_empty() {
            println!("  {}", o.notes);
        }
    }

    println!();
    println!("Endgame catalog");
    println!("---------------");
    for e in endgames.endgames {
        println!("{} | {} | tags: {}", e.id, e.name, e.tags.join(","));
        if !e.notes.is_empty() {
            println!("  {}", e.notes);
        }
    }

    Ok(())
}

fn usage() {
    println!("Usage:");
    println!("  cargo run -- stats");
    println!("  cargo run -- kernel");
    println!("  cargo run -- catalog");
}

fn main() -> Result<()> {
    let cmd = env::args().nth(1).unwrap_or_else(|| "stats".to_string());
    match cmd.as_str() {
        "stats" => cmd_stats(),
        "kernel" => cmd_kernel(),
        "catalog" => cmd_catalog(),
        "-h" | "--help" | "help" => {
            usage();
            Ok(())
        }
        other => Err(anyhow!("unknown command: {other}")),
    }
}
