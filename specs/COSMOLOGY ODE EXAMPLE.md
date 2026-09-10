# TASK: CREATE COSMOLOGY ODE EXAMPLE IN RUSTY-SUNDIALS

I need you to write a new example file `examples/cosmology_quintessence.rs` for my `Rusty-SUNDIALS` solver. We are integrating the highly stiff Einstein-Klein-Gordon equations for a scalar field in a hyperbolic target space over cosmological time.

## 1. Setup & ODE System
The state vector `y_vec` has 5 components: `[a, x, y, u, v]`, where `tau = x + i y`, and `u = dx/dt`, `v = dy/dt`.
The metric is Poincaré: $ds^2 = (dx^2 + dy^2)/(2y^2)$.

Create a function `compute_potential(x: f64, y: f64) -> (f64, f64, f64)` returning `(V, dV/dx, dV/dy)`.
- Model a smooth potential with a saddle at `x=0, y = 1/sqrt(12)` (Fricke point, high energy) and a global minimum at `x=0.5, y = sqrt(3)/2` (Orbifold point, low energy).

## 2. RHS Function (cosmology_rhs)
Match our CVODE RHS trait. Inside:
1. Clamp `y` to be strictly > 0 (e.g., `max(1e-8)` to prevent singularity).
2. Calculate kinetic energy: `T = (u^2 + v^2) / (2 * y^2)`.
3. Background densities: `rho_m = 0.315 / a^3`, `rho_r = 9.2e-5 / a^4`.
4. Total density: `rho_tot = T + V + rho_m + rho_r`.
5. Friedmann Hubble parameter: `H = sqrt(rho_tot / 3.0)`.
6. `ydot[0] = a * H`
7. `ydot[1] = u`
8. `ydot[2] = v`
9. `ydot[3] = (2.0 / y) * u * v - 3.0 * H * u - y^2 * dVdx`
10. `ydot[4] = (u^2 - v^2) / y - 3.0 * H * v - y^2 * dVdy`

## 3. Solver Configuration
- Use `CVode::new(Bdf, Newton)`. This is a hyper-stiff problem.
- Initial conditions (post-inflation): `a = 1e-10`, `x = 0.001`, `y = 1/sqrt(12) + 0.001`, `u=0`, `v=0`.
- Tolerances: Rel = 1e-8, Abs = 1e-10.
- Integrate step-by-step and output CSV data `(t, a, H, x, y)` to prove the modulus flows safely to the `(0.5, sqrt(3)/2)` attractor.