use bincode;
use clap::{Parser, Subcommand};
use flate2::read::GzDecoder;
use futures::stream::StreamExt;
use reqwest::Url;
use scraper::{Html, Selector};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

// -----------------------------------------------------------------------------
// UTILITIES & CONSTANTS
// -----------------------------------------------------------------------------

fn get_reader(filename: &str) -> Box<dyn BufRead> {
    let file = File::open(filename).unwrap_or_else(|_| panic!("Cannot open {}", filename));
    if filename.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

fn prob2score(p: f64) -> f64 {
    if p == 0.0 { -100.0 } else { (p / 0.25).log2() }
}

fn anti(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| match c {
        b'A'|b'a' => b'T', b'C'|b'c' => b'G', b'G'|b'g' => b'C', b'T'|b't' => b'A',
        b'R'|b'r' => b'Y', b'Y'|b'y' => b'R', b'K'|b'k' => b'M', b'M'|b'm' => b'K',
        b'W'|b'w' => b'W', b'S'|b's' => b'S', b'B'|b'b' => b'V', b'V'|b'v' => b'B',
        b'D'|b'd' => b'H', b'H'|b'h' => b'D',
        _ => c,
    }).collect()
}

// -----------------------------------------------------------------------------
// DATA STRUCTURES & CACHING
// -----------------------------------------------------------------------------

#[derive(Serialize, Deserialize, Clone, Debug)]
struct Feature {
    seqid: String, typ: String, beg: usize, end: usize, strand: char, fid: String, pid: String,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
struct Transcript {
    tx_feat: Feature, exons: Vec<Feature>, introns: Vec<Feature>, cdss: Vec<Feature>,
    utr5s: Vec<Feature>, utr3s: Vec<Feature>, is_coding: bool,
}

impl Transcript {
    fn finalize(&mut self) {
        self.is_coding = !self.cdss.is_empty();
    }
}

#[derive(Serialize, Deserialize, Default)]
struct GenomeData {
    sequences: HashMap<String, Vec<u8>>,
    genes: Vec<(Feature, Vec<Transcript>)>,
}

impl GenomeData {
    fn load_fasta(&mut self, path: &str, v: bool) {
        if v { println!("  [Verbose] Parsing FASTA: {}", path); }
        let reader = get_reader(path);
        let mut name = String::new();
        let mut seq = Vec::new();

        for line in reader.lines() {
            let line = line.unwrap();
            let line = line.trim();
            if line.starts_with('>') {
                if !name.is_empty() {
                    self.sequences.insert(name.clone(), seq.clone());
                    seq.clear();
                }
                name = line[1..].split_whitespace().next().unwrap().trim_start_matches("chr").to_string();
            } else {
                seq.extend_from_slice(line.as_bytes());
            }
        }
        if !name.is_empty() { self.sequences.insert(name, seq); }
    }

    fn load_gff3(&mut self, path: &str, v: bool) {
        if v { println!("  [Verbose] Parsing GFF3: {}", path); }
        let reader = gegt_reader(path);
        let mut features_by_pid: HashMap<String, Vec<Feature>> = HashMap::new();
        let mut genes_raw = Vec::new();

        for line in reader.lines() {
            let line = line.unwrap();
            if line.starts_with('#') || line.is_empty() { continue; }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 9 { continue; }

            let seqid = fields[0].trim_start_matches("chr").to_string();
            let typ = fields[2].to_string();
            let beg = fields[3].parse::<usize>().unwrap();
            let end = fields[4].parse::<usize>().unwrap();
            let strand = fields[6].chars().next().unwrap_or('.');
            
            let mut fid = String::new();
            let mut pids = vec![String::new()];

            for attr in fields[8].trim_end_matches(';').split(';') {
                if let Some((k, v)) = attr.split_once('=') {
                    if k == "ID" { fid = v.to_string(); }
                    if k == "Parent" { pids = v.split(',').map(|s| s.to_string()).collect(); }
                }
            }

            for pid in pids {
                let feat = Feature { seqid: seqid.clone(), typ: typ.clone(), beg, end, strand, fid: fid.clone(), pid: pid.clone() };
                if typ == "gene" { genes_raw.push(feat); } 
                else { features_by_pid.entry(pid).or_default().push(feat); }
            }
        }

        for gene_feat in genes_raw {
            let mut txs = Vec::new();
            if let Some(children) = features_by_pid.get(&gene_feat.fid) {
                for tx_feat in children {
                    let mut tx = Transcript {
                        tx_feat: tx_feat.clone(), exons: Vec::new(), introns: Vec::new(),
                        cdss: Vec::new(), utr5s: Vec::new(), utr3s: Vec::new(), is_coding: false,
                    };
                    if let Some(tx_children) = features_by_pid.get(&tx_feat.fid) {
                        for child in tx_children {
                            match child.typ.as_str() {
                                "exon" => tx.exons.push(child.clone()),
                                "intron" => tx.introns.push(child.clone()),
                                "CDS" => tx.cdss.push(child.clone()),
                                "five_prime_UTR" => tx.utr5s.push(child.clone()),
                                "three_prime_UTR" => tx.utr3s.push(child.clone()),
                                _ => {}
                            }
                        }
                    }
                    tx.finalize();
                    txs.push(tx);
                }
            }
            self.genes.push((gene_feat, txs));
        }
    }
}

fn build_or_load_cache(db_path: &str, fasta_path: Option<&str>, gff3_path: Option<&str>, force: bool, v: bool) -> GenomeData {
    let cache_file = format!("{}.bin", db_path);
    if !force && Path::new(&cache_file).exists() {
        if v { println!("  [Verbose] Loading cached database from {}", cache_file); }
        let file = File::open(&cache_file).unwrap();
        
        // FIXED: BufReader makes deserialization lightning fast!
        let reader = BufReader::new(file); 
        return bincode::deserialize_from(reader).expect("Failed to read cache");
    }

    let mut gd = GenomeData::default();
    gd.load_fasta(fasta_path.expect("FASTA required to build DB"), v);
    gd.load_gff3(gff3_path.expect("GFF3 required to build DB"), v);

    if v { println!("  [Verbose] Serializing database to {}", cache_file); }
    let file = File::create(&cache_file).unwrap();
    let writer = BufWriter::new(file);
    bincode::serialize_into(writer, &gd).unwrap();
    gd
}

// -----------------------------------------------------------------------------
// MARKOV MODEL
// -----------------------------------------------------------------------------

#[derive(Default)]
struct Model { k: usize, probs: HashMap<Vec<u8>, HashMap<u8, f64>>, log_odds_cache: HashMap<Vec<u8>, f64> }

impl Model {
    fn create(&mut self, seqs: &[Vec<u8>], k: usize) {
        self.k = k;
        let mut counts: HashMap<Vec<u8>, [u32; 4]> = HashMap::new();
        let nt_to_idx = |b: u8| -> Option<usize> {
            match b { b'A'|b'a'=>Some(0), b'C'|b'c'=>Some(1), b'G'|b'g'=>Some(2), b'T'|b't'=>Some(3), _=>None }
        };

        for seq in seqs {
            if seq.len() <= k { continue; }
            for i in k..seq.len() {
                if let Some(idx) = nt_to_idx(seq[i]) {
                    let kmer = seq[i - k..i].to_vec();
                    let entry = counts.entry(kmer).or_insert([1, 1, 1, 1]);
                    entry[idx] += 1;
                }
            }
        }

        let idx_to_nt = [b'A', b'C', b'G', b'T'];
        for (kmer, cnts) in counts {
            let total: u32 = cnts.iter().sum();
            let mut nt_probs = HashMap::new();
            for i in 0..4 {
                let p = cnts[i] as f64 / total as f64;
                nt_probs.insert(idx_to_nt[i], p);
                let mut cache_key = kmer.clone();
                cache_key.push(idx_to_nt[i]);
                self.log_odds_cache.insert(cache_key, prob2score(p));
            }
            self.probs.insert(kmer, nt_probs);
        }
    }

    fn score(&self, seq: &[u8]) -> f64 {
        if seq.len() <= self.k { return 0.0; }
        let mut total = 0.0;
        for window in seq.windows(self.k + 1) {
            total += self.log_odds_cache.get(window).unwrap_or(&-2.0);
        }
        total
    }

    fn export_file(&self, path: &str) {
        let mut file = BufWriter::new(File::create(path).unwrap());
        writeln!(file, "% MM {} {}", path, self.probs.len() * 4).unwrap();
        let mut sorted_keys: Vec<_> = self.probs.keys().collect();
        sorted_keys.sort();
        for kmer in sorted_keys {
            if let Some(nts) = self.probs.get(kmer) {
                let mut nts_sorted: Vec<_> = nts.iter().collect();
                nts_sorted.sort_by_key(|&(k, _)| k);
                for (nt, p) in nts_sorted {
                    writeln!(file, "{}{} {:.6}", String::from_utf8_lossy(kmer), *nt as char, p).unwrap();
                }
                writeln!(file).unwrap();
            }
        }
    }

    fn import_file(&mut self, path: &str) {
        let reader = get_reader(path);
        let mut k_found = None;
        for line in reader.lines() {
            let line = line.unwrap();
            if line.starts_with('%') || line.trim().is_empty() { continue; }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 2 { continue; }
            let key = parts[0].as_bytes();
            if key.len() < 2 { continue; }
            let (kmer, nt) = key.split_at(key.len() - 1);
            let val: f64 = parts[1].parse().unwrap();
            if k_found.is_none() { k_found = Some(kmer.len()); }
            let entry = self.probs.entry(kmer.to_vec()).or_default();
            entry.insert(nt[0], val);
            self.log_odds_cache.insert(key.to_vec(), prob2score(val));
        }
        self.k = k_found.unwrap_or(0);
    }
}

// -----------------------------------------------------------------------------
// ASYNC BULK PIPELINE (The Auto Crawler)
// -----------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct GenomeJob {
    prefix: String,
    fasta_url: String,
    gff_url: String,
}

/// Crawls a directory URL and heuristically pairs FASTA and GFF3 files
async fn crawl_database(base_url_str: &str) -> Result<Vec<GenomeJob>, Box<dyn std::error::Error>> {
    println!("Crawling {} ...", base_url_str);
    let base_url = Url::parse(base_url_str)?;
    let html_content = reqwest::get(base_url.clone()).await?.text().await?;
    
    let document = Html::parse_document(&html_content);
    let selector = Selector::parse("a").unwrap();
    
    let mut fastas = Vec::new();
    let mut gffs = Vec::new();

    for element in document.select(&selector) {
        if let Some(href) = element.value().attr("href") {
            let full_url = base_url.join(href)?.to_string();
            if href.ends_with(".fna.gz") || href.ends_with(".fa.gz") {
                fastas.push(full_url);
            } else if href.ends_with(".gff.gz") || href.ends_with(".gff3.gz") {
                gffs.push(full_url);
            }
        }
    }

    let mut jobs = Vec::new();
    
    for fasta_url in fastas {
        let f_name = fasta_url.split('/').last().unwrap();
        let prefix = f_name.replace(".fna.gz", "").replace(".fa.gz", "");
        if let Some(gff_url) = gffs.iter().find(|g| g.contains(&prefix)) {
            jobs.push(GenomeJob {
                prefix: prefix.clone(), fasta_url: fasta_url.clone(), gff_url: gff_url.clone(),
            });
        }
    }
    
    Ok(jobs)
}

/// Streams a URL directly to the disk
async fn download_file(url: &str, dest: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut response = reqwest::get(url).await?;
    if !response.status().is_success() {
        return Err(format!("Failed to download {}: Status {}", url, response.status()).into());
    }
    let mut file = std::fs::File::create(dest)?;
    while let Some(chunk) = response.chunk().await? {
        std::io::Write::write_all(&mut file, &chunk)?;
    }
    Ok(())
}

/// Worker that processes a single genome end-to-end
async fn process_genome_worker(job: GenomeJob, k: usize, v: bool) -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    println!(">> Starting Genome: {}", job.prefix);
    
    std::fs::create_dir_all("bulk_out")?;
    
    let fasta_path = format!("{}_temp.fa.gz", job.prefix);
    let gff_path = format!("{}_temp.gff.gz", job.prefix);
    let db_path = format!("{}_db", job.prefix);

    // 1. Download
    download_file(&job.fasta_url, &fasta_path).await?;
    download_file(&job.gff_url, &gff_path).await?;
    
    // 2. Offload processing to background thread
    let prefix_clone = job.prefix.clone();
    let res = tokio::task::spawn_blocking(move || {
        let gd = build_or_load_cache(&db_path, Some(&fasta_path), Some(&gff_path), true, v);
        
        let mut mmbase = Model::default();
        let seqs: Vec<Vec<u8>> = gd.sequences.values().cloned().collect();
        mmbase.create(&seqs, k);

        let mut models = HashMap::new();
        let configs: [(&str, fn(&Transcript) -> &Vec<Feature>); 4] = [
            ("exon", |tx| &tx.exons), ("intron", |tx| &tx.introns),
            ("utr5", |tx| &tx.utr5s), ("utr3", |tx| &tx.utr3s),
        ];

        for (name, accessor) in configs {
            let mut feat_seqs = Vec::new();
            for (_, txs) in &gd.genes {
                for tx in txs {
                    for feat in accessor(tx) {
                        if let Some(seq) = gd.sequences.get(&feat.seqid) {
                            if feat.beg > 0 && feat.end <= seq.len() {
                                let f_seq = &seq[feat.beg - 1..feat.end];
                                feat_seqs.push(if feat.strand == '-' { anti(f_seq) } else { f_seq.to_vec() });
                            }
                        }
                    }
                }
            }
            if !feat_seqs.is_empty() {
                let mut mm = Model::default();
                mm.create(&feat_seqs, k);
                models.insert(name, mm);
            }
        }

        let mut badgenes = Vec::new();
        for (gene, txs) in &gd.genes {
            for tx in txs {
                if !tx.is_coding { continue; }
                
                let (mut t_model_sc, mut t_model_pos, mut t_base_sc, mut t_base_pos) = (0.0, 0, 0.0, 0);

                for (m_key, accessor) in configs {
                    if let Some(model) = models.get(m_key) {
                        for feat in accessor(tx) {
                            if let Some(seq) = gd.sequences.get(&feat.seqid) {
                                if feat.beg > 0 && feat.end <= seq.len() {
                                    let mut f_seq = seq[feat.beg - 1..feat.end].to_vec();
                                    if feat.strand == '-' { f_seq = anti(&f_seq); }

                                    t_model_sc += model.score(&f_seq);
                                    t_model_pos += f_seq.len().saturating_sub(model.k);
                                    t_base_sc += mmbase.score(&f_seq);
                                    t_base_pos += f_seq.len().saturating_sub(mmbase.k);
                                }
                            }
                        }
                    }
                }

                if t_model_pos > 0 && t_base_pos > 0 {
                    let diff = (t_base_sc / t_base_pos as f64) - (t_model_sc / t_model_pos as f64);
                    if diff > 0.0 { badgenes.push((gene.fid.clone(), tx.tx_feat.fid.clone(), diff)); }
                }
            }
        }

        badgenes.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
        let out_file = format!("bulk_out/{}_badgenes.txt", prefix_clone);
        let mut out = BufWriter::new(File::create(&out_file).unwrap());
        writeln!(out, "GeneID\tTranscriptID\tScoreDiff").unwrap();
        for (gid, tid, diff) in &badgenes {
            writeln!(out, "{}\t{}\t{:.4}", gid, tid, diff).unwrap();
        }
        println!("<< Analysis complete. Found {} anomalous transcripts. Results written to {}", badgenes.len(), out_file);

        // 3. Cleanup temp files inside the closure! 
        // (Since they were moved in here, the thread owns them and can delete them safely)
        let _ = std::fs::remove_file(&fasta_path);
        let _ = std::fs::remove_file(&gff_path);
        let _ = std::fs::remove_file(format!("{}.bin", db_path));

    }).await;

    res?;
    println!("<< Finished Genome: {} in {:.2?}", job.prefix, start.elapsed());
    Ok(())
}

// -----------------------------------------------------------------------------
// CLI ENTRY POINT
// -----------------------------------------------------------------------------

#[derive(Parser)]
#[command(name = "gene_qc")]
#[command(about = "Optimized Gene Quality Control Pipeline (Local & Auto)")]
struct Cli {
    #[arg(short, long, global = true)]
    verbose: bool,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Create a local binary cache from a FASTA and GFF3 file
    Createdb {
        db: String,
        fasta: String,
        gff3: String,
    },
    /// Train Markov models from a local database
    Train {
        db: String,
        #[arg(short, long, default_value_t = 3)]
        k: usize,
    },
    /// Score and find bad genes using local models
    Findbgs {
        db: String,
    },
    /// Crawl a URL, download, parse, and score multiple genomes concurrently
    Auto {
        url: String,
        #[arg(short, long, default_value_t = 3)]
        k: usize,
    },
}

#[tokio::main]
async fn main() {
    let cli = Cli::parse();
    let v = cli.verbose;

    match cli.command {
        // --- LOCAL: Create DB ---
        Commands::Createdb { db, fasta, gff3 } => {
            println!("Building database...");
            let start = Instant::now();
            build_or_load_cache(&db, Some(&fasta), Some(&gff3), true, v);
            println!("Database saved to {}.bin in {:.2?}", db, start.elapsed());
        }
        
        // --- LOCAL: Train Models ---
        Commands::Train { db, k } => {
            println!("Loading database...");
            let gd = build_or_load_cache(&db, None, None, false, v);
            std::fs::create_dir_all("models").unwrap();

            println!("Training baseline model...");
            let mut mmbase = Model::default();
            let seqs: Vec<Vec<u8>> = gd.sequences.values().cloned().collect();
            mmbase.create(&seqs, k);
            mmbase.export_file("models/baseline.mm");

            // FIXED: Typed closures as function pointers
            let configs: [(&str, &str, fn(&Transcript) -> &Vec<Feature>); 4] = [
                ("exon", "exon.mm", |tx| &tx.exons),
                ("intron", "intron.mm", |tx| &tx.introns),
                ("5' UTR", "5utr.mm", |tx| &tx.utr5s),
                ("3' UTR", "3utr.mm", |tx| &tx.utr3s),
            ];

            for (name, filename, accessor) in configs {
                println!("Training {} model...", name);
                let mut feat_seqs = Vec::new();
                for (_, txs) in &gd.genes {
                    for tx in txs {
                        for feat in accessor(tx) {
                            if let Some(seq) = gd.sequences.get(&feat.seqid) {
                                if feat.beg > 0 && feat.end <= seq.len() {
                                    let f_seq = &seq[feat.beg - 1..feat.end];
                                    feat_seqs.push(if feat.strand == '-' { anti(f_seq) } else { f_seq.to_vec() });
                                }
                            }
                        }
                    }
                }
                if !feat_seqs.is_empty() {
                    let mut mm = Model::default();
                    mm.create(&feat_seqs, k);
                    mm.export_file(&format!("models/{}", filename));
                }
            }
            println!("Training complete! Models saved in 'models/' directory.");
        }

        // --- LOCAL: Find Bad Genes ---
        Commands::Findbgs { db } => {
            println!("Loading database...");
            let gd = build_or_load_cache(&db, None, None, false, v);
            
            println!("Loading models...");
            let mut models = HashMap::new();
            let model_files = [("base", "baseline.mm"), ("exon", "exon.mm"), ("intron", "intron.mm"), ("utr5", "5utr.mm"), ("utr3", "3utr.mm")];
            
            for (key, file) in model_files {
                let path = format!("models/{}", file);
                if Path::new(&path).exists() {
                    let mut m = Model::default();
                    m.import_file(&path);
                    models.insert(key, m);
                }
            }

            let base_model = models.get("base").expect("Missing baseline model");
            let mut badgenes = Vec::new();

            println!("Scoring genes...");
            let checks: [(&str, fn(&Transcript) -> &Vec<Feature>); 4] = [
                ("exon", |t| &t.exons), ("intron", |t| &t.introns),
                ("utr5", |t| &t.utr5s), ("utr3", |t| &t.utr3s)
            ];

            for (gene, txs) in &gd.genes {
                for tx in txs {
                    if !tx.is_coding { continue; }
                    
                    let (mut t_model_sc, mut t_model_pos, mut t_base_sc, mut t_base_pos) = (0.0, 0, 0.0, 0);

                    for (m_key, accessor) in checks {
                        if let Some(model) = models.get(m_key) {
                            for feat in accessor(tx) {
                                if let Some(seq) = gd.sequences.get(&feat.seqid) {
                                    if feat.beg > 0 && feat.end <= seq.len() {
                                        let mut f_seq = seq[feat.beg - 1..feat.end].to_vec();
                                        if feat.strand == '-' { f_seq = anti(&f_seq); }

                                        t_model_sc += model.score(&f_seq);
                                        t_model_pos += f_seq.len().saturating_sub(model.k);
                                        t_base_sc += base_model.score(&f_seq);
                                        t_base_pos += f_seq.len().saturating_sub(base_model.k);
                                    }
                                }
                            }
                        }
                    }

                    if t_model_pos > 0 && t_base_pos > 0 {
                        let diff = (t_base_sc / t_base_pos as f64) - (t_model_sc / t_model_pos as f64);
                        if diff > 0.0 { badgenes.push((gene.fid.clone(), tx.tx_feat.fid.clone(), diff)); }
                    }
                }
            }

            badgenes.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
            let mut out = BufWriter::new(File::create("badgenes.txt").unwrap());
            writeln!(out, "GeneID\tTranscriptID\tScoreDiff").unwrap();
            
            // FIXED: &badgenes prevents borrow-checker panic
            for (gid, tid, diff) in &badgenes {
                writeln!(out, "{}\t{}\t{:.4}", gid, tid, diff).unwrap();
            }
            println!("Analysis complete. Found {} anomalous transcripts. Results written to badgenes.txt", badgenes.len());
        }

        // --- ASYNC AUTO: Online Web Crawler ---
        Commands::Auto { url, k } => {
            println!("Initializing bulk pipeline...");
            
            let jobs = match crawl_database(&url).await {
                Ok(j) => j,
                Err(e) => { eprintln!("Failed to crawl directory: {}", e); return; }
            };
            
            println!("Found {} valid FASTA/GFF3 pairs to process.", jobs.len());
            
            let fetches = futures::stream::iter(jobs)
                .map(|job| {
                    tokio::spawn(async move {
                        if let Err(e) = process_genome_worker(job.clone(), k, v).await {
                            eprintln!("Error processing {}: {:?}", job.prefix, e);
                        }
                    })
                })
                .buffer_unordered(4); // limit to 4 concurrent downloads/runs
                
            fetches.collect::<Vec<_>>().await;
            println!("Bulk pipeline complete! Check the 'bulk_out/' directory.");
        }
    }
}