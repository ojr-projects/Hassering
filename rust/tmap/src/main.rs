// %% tmap: exact sparse-pattern spectra for the genome tendency map (docs/GENOME_TENDENCY_MAP_REVIEW.md, Rev. 2)
// Window-free, position-anchored and exact. Every base anchors every pattern; N (including masked coding sequence
// in the non-coding mode) breaks any pattern whose footprint touches it.
//
//   tmap spectrum <genomes.txt> <out_prefix> [--core 4] [--gmax 96] [--threads N]
//
// genomes.txt: one FASTA path per line (one genome per strain, pooled). Outputs, little-endian:
//   <prefix>.cores.f32   per gapped shape (k = 4, 5; span <= 13; first offset 0; itertools order) and cell:
//                        [observed, expected] pairs (both strands; expected from the pooled order-1 Markov chain)
//   <prefix>.pairs.u32   offset-pair spectrum: count[g][a][b] = #(core a at i, core b at i + c + g), g = 0..gmax,
//                        contiguous cores of length c, both strands, whole footprint N-free
//   <prefix>.json        dimensions, shapes, model, positions
use rayon::prelude::*;
use std::io::{BufRead, BufReader, Write};

const W: usize = 13;

fn read_fasta(path: &str) -> Vec<Vec<u8>> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut cur: Option<Vec<u8>> = None;
    for line in BufReader::new(f).lines() {
        let line = line.unwrap();
        if line.starts_with('>') {
            if let Some(s) = cur.take() {
                out.push(s);
            }
            cur = Some(Vec::new());
        } else if let Some(s) = cur.as_mut() {
            s.extend(line.trim().bytes().map(|b| match b {
                b'A' | b'a' => 0u8,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            }));
        }
    }
    if let Some(s) = cur.take() {
        out.push(s);
    }
    out
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&c| if c < 4 { 3 - c } else { 4 }).collect()
}

/// Both orientations of every record of every genome.
fn strands(genomes: &[Vec<Vec<u8>>]) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    for g in genomes {
        for r in g {
            out.push(r.clone());
            out.push(revcomp(r));
        }
    }
    out
}

#[inline]
fn pext(x: u32, m: u32) -> u32 {
    #[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))]
    unsafe {
        return core::arch::x86_64::_pext_u32(x, m);
    }
    #[allow(unreachable_code)]
    {
        let (mut r, mut bb, mut mm) = (0u32, 1u32, m);
        while mm != 0 {
            let low = mm & mm.wrapping_neg();
            if x & low != 0 {
                r |= bb;
            }
            bb <<= 1;
            mm &= mm - 1;
        }
        r
    }
}

struct Shape {
    offs: Vec<usize>,
    k: usize,
    bm: u32,
    pm: u16,
}

fn shape(offs: &[usize]) -> Shape {
    let (mut bm, mut pm) = (0u32, 0u16);
    for &o in offs {
        bm |= 3u32 << (2 * (W - 1 - o));
        pm |= 1u16 << (W - 1 - o);
    }
    Shape { offs: offs.to_vec(), k: offs.len(), bm, pm }
}

fn all_shapes(k: usize) -> Vec<Shape> {
    fn rec(start: usize, left: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if left == 0 {
            out.push(cur.clone());
            return;
        }
        for v in start..W {
            cur.push(v);
            rec(v + 1, left - 1, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    rec(1, k - 1, &mut vec![0], &mut out);
    out.iter().map(|o| shape(o)).collect()
}

/// 26-bit window codes and 13-bit N masks of one strand.
fn window_codes(a: &[u8]) -> (Vec<u32>, Vec<u16>) {
    let n = a.len();
    let get = |i: usize| if i < n { a[i] } else { 4 };
    let (mut c, mut m) = (0u32, 0u16);
    for j in 0..W {
        let b = get(j);
        c = (c << 2) | (b & 3) as u32;
        m = (m << 1) | (b == 4) as u16;
    }
    let cmask: u32 = (1u32 << (2 * W)) - 1;
    let mmask: u16 = (1u16 << W) - 1;
    let mut code = vec![0u32; n];
    let mut nm = vec![0u16; n];
    for t in 0..n {
        code[t] = c;
        nm[t] = m;
        let b = get(t + W);
        c = ((c << 2) | (b & 3) as u32) & cmask;
        m = ((m << 1) | (b == 4) as u16) & mmask;
    }
    (code, nm)
}

/// Contiguous core codes of length c starting at each position, with validity (no N in the core).
fn core_codes(a: &[u8], c: usize) -> (Vec<u32>, Vec<bool>) {
    let n = a.len();
    let mut code = vec![0u32; n];
    let mut ok = vec![false; n];
    let mask: u32 = if 2 * c >= 32 { u32::MAX } else { (1u32 << (2 * c)) - 1 };
    let (mut x, mut run) = (0u32, 0usize);
    for i in 0..n {
        let b = a[i];
        if b == 4 {
            run = 0;
            x = 0;
        } else {
            x = ((x << 2) | b as u32) & mask;
            run += 1;
        }
        if run >= c {
            let s = i + 1 - c;
            code[s] = x;
            ok[s] = true;
        }
    }
    (code, ok)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str, d: &str| -> String {
        args.iter().position(|a| a == k).map(|i| args[i + 1].clone()).unwrap_or_else(|| d.to_string())
    };
    if args.len() < 4 || args[1] != "spectrum" {
        eprintln!("usage: tmap spectrum <genomes.txt> <out_prefix> [--core 4] [--gmax 96] [--threads N]");
        std::process::exit(2);
    }
    let threads: usize = get("--threads", "0").parse().unwrap();
    if threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    }
    let c: usize = get("--core", "4").parse().unwrap();
    let gmax: usize = get("--gmax", "96").parse().unwrap();
    let list = std::fs::read_to_string(&args[2]).unwrap();
    let paths: Vec<&str> = list.lines().map(|l| l.trim()).filter(|l| !l.is_empty()).collect();
    let t0 = std::time::Instant::now();
    let genomes: Vec<Vec<Vec<u8>>> = paths.par_iter().map(|p| read_fasta(p)).collect();
    let st = strands(&genomes);
    let valid: u64 = st.iter().map(|s| s.iter().filter(|&&b| b < 4).count() as u64).sum();
    eprintln!("{} genomes, {} strands, {} valid bases (both strands) read in {:.1}s", genomes.len(), st.len(), valid,
              t0.elapsed().as_secs_f64());

    // ------------------------------------------------ cores: every gapped shape k = 4, 5 (span <= 13)
    let t1 = std::time::Instant::now();
    let shapes: Vec<Shape> = all_shapes(4).into_iter().chain(all_shapes(5)).collect();
    let wc: Vec<(Vec<u32>, Vec<u16>)> = st.par_iter().map(|s| window_codes(s)).collect();
    let obs: Vec<Vec<u32>> = shapes
        .par_iter()
        .map(|sh| {
            let mut h = vec![0u32; 1 << (2 * sh.k)];
            for (code, nm) in &wc {
                for (cc, m) in code.iter().zip(nm.iter()) {
                    if m & sh.pm == 0 {
                        h[pext(*cc, sh.bm) as usize] += 1;
                    }
                }
            }
            h
        })
        .collect();
    drop(wc);
    // order-1 Markov from the dinucleotides of the same strands
    let mut di = [[0.5f64; 4]; 4];
    for s in &st {
        for w in s.windows(2) {
            if w[0] < 4 && w[1] < 4 {
                di[w[0] as usize][w[1] as usize] += 1.0;
            }
        }
    }
    let tot: f64 = di.iter().flatten().sum();
    let mut pi = [0.0; 4];
    let mut t = [[0.0; 4]; 4];
    for a in 0..4 {
        let rs: f64 = di[a].iter().sum();
        pi[a] = rs / tot;
        for b in 0..4 {
            t[a][b] = di[a][b] / rs;
        }
    }
    let mut tp = vec![[[0.0f64; 4]; 4]; W + 1];
    for i in 0..4 {
        tp[0][i][i] = 1.0;
    }
    for g in 1..=W {
        for i in 0..4 {
            for j in 0..4 {
                tp[g][i][j] = (0..4).map(|k| tp[g - 1][i][k] * t[k][j]).sum();
            }
        }
    }
    let mut fh = std::io::BufWriter::new(std::fs::File::create(format!("{}.cores.f32", args[3])).unwrap());
    for (sh, h) in shapes.iter().zip(&obs) {
        let n: f64 = h.iter().map(|&x| x as f64).sum();
        let mut v: Vec<f64> = pi.to_vec();
        for i in 1..sh.k {
            let g = &tp[sh.offs[i] - sh.offs[i - 1]];
            let mut nv = vec![0.0; v.len() * 4];
            for (idx, &x) in v.iter().enumerate() {
                for b in 0..4 {
                    nv[idx * 4 + b] = x * g[idx & 3][b];
                }
            }
            v = nv;
        }
        for (o, p) in h.iter().zip(v.iter()) {
            fh.write_all(&(*o as f32).to_le_bytes()).unwrap();
            fh.write_all(&((p * n) as f32).to_le_bytes()).unwrap();
        }
    }
    fh.flush().unwrap();
    let t_cores = t1.elapsed().as_secs_f64();

    // ------------------------------------------------ offset-pair spectrum of contiguous cores
    let t2 = std::time::Instant::now();
    let cc = 1usize << (2 * c);
    let cells = (gmax + 1) * cc * cc;
    let pairs: Vec<u32> = st
        .par_iter()
        .fold(
            || vec![0u32; cells],
            |mut h, s| {
                let (code, ok) = core_codes(s, c);
                let n = s.len();
                // next N at or after each position: the pair's whole footprint (cores and gap) must be N-free, so
                // in the non-coding mode no pair straddles masked coding sequence
                let mut nn = vec![n; n + 1];
                for i in (0..n).rev() {
                    nn[i] = if s[i] == 4 { i } else { nn[i + 1] };
                }
                for i in 0..n {
                    if !ok[i] {
                        continue;
                    }
                    let a = code[i] as usize * cc;
                    for g in 0..=gmax {
                        let j = i + c + g;
                        if j + c > nn[i] {
                            break;
                        }
                        if ok[j] {
                            h[g * cc * cc + a + code[j] as usize] += 1;
                        }
                    }
                }
                h
            },
        )
        .reduce(|| vec![0u32; cells], |mut a, b| {
            for (x, y) in a.iter_mut().zip(b) {
                *x += y;
            }
            a
        });
    let mut fp = std::io::BufWriter::new(std::fs::File::create(format!("{}.pairs.u32", args[3])).unwrap());
    for x in &pairs {
        fp.write_all(&x.to_le_bytes()).unwrap();
    }
    fp.flush().unwrap();
    let t_pairs = t2.elapsed().as_secs_f64();
    let shp: Vec<String> = shapes.iter().map(|s| format!("{:?}", s.offs)).collect();
    let js = format!(
        "{{\"genomes\": {}, \"valid_bases\": {}, \"core\": {}, \"gmax\": {}, \"n_shapes\": {}, \"shapes\": [{}], \
         \"pi\": {:?}, \"t\": {:?}, \"time_s\": {{\"cores\": {:.2}, \"pairs\": {:.2}}}}}",
        genomes.len(), valid, c, gmax, shapes.len(), shp.join(", "), pi, t, t_cores, t_pairs);
    std::fs::write(format!("{}.json", args[3]), js).unwrap();
    eprintln!("cores {:.1}s, pairs {:.1}s", t_cores, t_pairs);
}
