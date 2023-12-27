use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;


fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn tajd_cstes(nb: u64) -> (f64, f64) {
	let mut a1: f64 = 0.0;
	let mut a2: f64 = 0.0;
	for n in 1..nb {
		a1 += 1.0 / n as f64;
		a2 += 1.0 / (n.pow(2) as f64);
	}

	let b1: f64 = (nb + 1) as f64 / ((3 * (nb - 1)) as f64);
	let b2: f64 = ((2 * (nb.pow(2) + nb + 3)) as f64 )/ ((9 * nb * (nb - 1)) as f64);
	let c1: f64 = b1 - (1.0 / a1);
	let c2: f64 = b2 - ((nb + 2) as f64 / (a1 * nb as f64)) + (a2 / a1.powi(2));
	let e1: f64 = c1 / a1;
	let e2: f64 = c2 / (a1.powi(2) + a2);
	return (e1, e2);
}

fn main() {
  let all_nb: u64 = 242824;
  let mut pi_tmp: u64 = 0;
  let mut h_tmp: u64 = 0;
  let mut seg: u64 = 0;
  let mut tajd: f64 = 0.0;
  let det: f64 = 8276.0;

  // calculate a_1
  let mut a_1: f64 = 0.0;
  for m in 1..(all_nb+1) {
    a_1 += 1.0 / m as f64;
  }
  println!("a_1 calculated");
  // Read in pi values
  let mut linnum: u8 = 0;
  if let Ok(lines) = read_lines("./pykeys_new.txt") {
    println!("reading file");
	for line in lines {
        linnum += 1;
		if let Ok(ip) = line {
			let mut spl = ip.split_whitespace();
			let snp = spl.next().unwrap().parse::<u64>().unwrap();
			let nb = spl.next().unwrap().parse::<u32>().unwrap();
			let freq = spl.next().unwrap().parse::<u64>().unwrap();
	        pi_tmp += freq * (all_nb - freq) * snp;
	        h_tmp += freq.pow(2) * snp;
	        seg += snp;
		}
	}
  }
  println!("pi calculated");

  let tipi: u64 = all_nb;
  let pi = pi_tmp as f64 / ((tipi as f64 * (tipi as f64 - 1.0)) / 2.0) / det ;
  let pi_corrected = 0.75 * (1.0 - ( 4.8 * pi / 3.0)).ln();
  let theta = seg as f64 / a_1 / det;
  let thetaW_corrected = ( det / (det - seg as f64)).ln() / a_1;

  let mut d: f64 = 0.0;
  let (e1, e2) = tajd_cstes(all_nb);
  let p1p2: f64 = e1 * seg as f64 + e2 * seg as f64 * (seg as f64 - 1.0);
  if pi > 0.0 && theta > 0.0 {
	d = pi - theta;
  }

  if p1p2 > 0.0 {
	tajd = d * det / p1p2.sqrt();
  }

  let hache = h_tmp as f64/ ((tipi as f64 * (tipi as f64 - 1.0)) / 2.0);
  let final_h = pi - hache;

  println!("pi:\t{}", pi);
  println!("pi Corrected:\t{}", pi_corrected);
  println!("thetaW:\t{}", theta);
  println!("thetaW Corrected:\t{}", thetaW_corrected);
  println!("D:\t{}", d);
  println!("Tajima's D:\t{}", tajd);
  println!("H:\t{}", final_h);
}
