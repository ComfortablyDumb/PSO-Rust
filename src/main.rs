use rand::Rng;

fn min(a: f64, b: f64) -> f64 {
    if a > b {
        return b;
    }
    return a;
}

fn obj_func(p: &Vec<f64>) -> f64 {
    let mut res: f64 = 10. * (p.len() as f64);
    for xi in p {
        res += xi.powf(2.) - 10. * (2. * std::f64::consts::PI * xi).cos();
    }

    return res;
}

fn norm(p:&Vec<f64>)->f64{
    let mut res:f64  = 0.;

    for xi in p{
        res += xi.powf(2.);
    }
    return res;

}

fn join(p: &Vec<f64>, c: &str) -> String {
    let mut res = String::from("(");
    for xi in p {
        res = res + &xi.to_string() + &c;
    }

    res.pop();
    res += ")";

    return res;
}

fn pso<F: Fn(&Vec<f64>) -> f64>(
    f: F,
    d: usize,
    n_particles: usize,
    phi1: f64,
    phi2: f64,
    n_max: f64,
    eps: f64,
    w_max: f64,
    w_min: f64,
    ub: f64,
    lb: f64,
) {
    let mut rng = rand::thread_rng();
    let dw = (w_max - w_min) / n_max;
    let b = min(ub.abs(), lb.abs());
    let mut particles = vec![vec![0.0f64; d]; n_particles];
    let mut velocities = vec![vec![0.0f64; d]; n_particles];

    for i in 0..n_particles {
        for j in 0..particles[i].len() {
            particles[i][j] = rng.gen_range(-b, b);
            velocities[i][j] = rng.gen_range(-b, b);
        }
    }

    let mut pb = particles.clone();
    let mut gb = pb[0].clone();

    for i in 0..n_particles {
        if f(&pb[i]) < f(&gb) {
            gb = pb[i].clone();
        }
    }

    let mut w = w_max;

    let mut n: f64 = 1.;
    let mut u1: f64;
    let mut u2: f64;

    while n < n_max && f(&gb).abs() > eps {
        for p in 0..n_particles {
            u1 = rng.gen();
            u2 = rng.gen();

            for i in 0..d {
                velocities[p][i] = velocities[p][i] * w
                    + phi1 * u1 * (pb[p][i] - particles[p][i])
                    + phi2 * u2 * (gb[i] - particles[p][i]);
                particles[p][i] += velocities[p][i]
            }

            if f(&pb[p]) > f(&particles[p]) {
                pb[p] = particles[p].clone();
                if f(&gb) > f(&pb[p]) {
                    gb = pb[p].clone();
                }
            }
        }
        w -= dw;
        n += 1.;
    }

    println!(
        "gb:{}\nf(gb):{}\nn:{}",
        join(&gb, ","),
        f(&gb),
        n
    );
}

fn main() {

    pso(obj_func,2,1000,1.,1.,1000.,1E-10,0.9,0.3,1.,-1.);
}
