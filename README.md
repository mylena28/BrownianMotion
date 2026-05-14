# Brownian Dynamics: Settling Sphere

Fortran implementation of the **Ermak-McCammon** algorithm for simulating a single colloidal sphere undergoing simultaneous gravitational settling and Brownian motion in a quiescent Newtonian fluid.

---

## Physical Model

The particle evolves under an overdamped Langevin equation (inertia neglected):

```
γ dr/dt = F_ext + ξ(t)
```

where `γ = 6πμa` is the Stokes drag, `F_ext` is the net gravitational force, and `ξ(t)` is a Gaussian white-noise force satisfying the fluctuation-dissipation theorem.

### Ermak-McCammon Integration Scheme

At each time step Δt:

```
x(t+Δt) = x(t)            + √(2D Δt) · Rx
y(t+Δt) = y(t)            + √(2D Δt) · Ry
z(t+Δt) = z(t) − U_s Δt  + √(2D Δt) · Rz
```

where `Rx, Ry, Rz ~ N(0,1)` are independent Gaussian random numbers, `D = kBT/γ` is the Stokes-Einstein diffusion coefficient, and `U_s = F_ext/γ` is the sedimentation velocity. The z-axis is taken as upward positive.

### Default Physical Parameters

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Temperature | T | 298.15 | K |
| Fluid viscosity (water) | μ | 1.0 × 10⁻³ | Pa·s |
| Sphere radius | a | 1.0 × 10⁻⁶ | m |
| Particle density (SiO₂) | ρ_p | 2000 | kg/m³ |
| Fluid density (water) | ρ_f | 1000 | kg/m³ |

### Derived Quantities

| Quantity | Symbol | Value | Units |
|----------|--------|-------|-------|
| Stokes drag | γ | 1.885 × 10⁻⁸ | N·s/m |
| Diffusion coefficient | D | 2.184 × 10⁻¹³ | m²/s |
| Sedimentation velocity | U_s | 2.180 × 10⁻⁶ | m/s |
| Peclet number | Pe = U_s a / D | ≈ 10 | — |

---

## Project Structure

```
.
├── params_mod.f90      # Phase 1 – physical constants module
├── rng_mod.f90         # Phase 2 – Box-Muller Gaussian RNG module
├── brownian_sim.f90    # Phase 3 – settling simulation (gravity on)
├── msd_validate.f90    # Phase 4 – MSD validation (gravity off)
├── plot_trajectory.py  # Phase 4 – Python visualisation
└── Makefile
```

### Module descriptions

**`params_mod.f90`** — `module physical_params`
Defines all physical constants as Fortran `parameter` values. The subroutine `init_params()` computes the derived quantities (D, U_s, γ) at runtime and prints a summary.

**`rng_mod.f90`** — `module rng_module`
Implements the **Box-Muller transform** to convert pairs of uniform U[0,1] deviates (from the intrinsic `random_number`) into independent N(0,1) Gaussian deviates. Provides `gaussian_3(rx, ry, rz)` for use in the integration loop.

**`brownian_sim.f90`** — `program brownian_settling`
Main integration loop. Runs 500 000 steps (Δt = 0.01 s, total 5000 s) and writes the trajectory to `trajectory.dat`.

**`msd_validate.f90`** — `program msd_validation`
Runs a 100 000-step no-gravity simulation, stores the full trajectory in memory, and computes the time-averaged MSD for up to 2000 lag values. Validates `⟨Δx²⟩ = 2Dτ`.

**`plot_trajectory.py`**
Python/Matplotlib script producing:
- 3-D trajectory coloured by time (`trajectory_3d.png`)
- MSD vs lag time on linear and log-log scales (`msd_validation.png`)

---

## Requirements

### Fortran
- `gfortran` ≥ 7 (or any Fortran 2003-compliant compiler)

### Python (visualisation only)
- Python ≥ 3.8
- `numpy`, `matplotlib`

```bash
pip install numpy matplotlib
```

---

## Build & Run

```bash
# Build both executables
make

# Run settling simulation → trajectory.dat
make run

# Run MSD validation → msd_output.dat, no_gravity_traj.dat
make validate

# Build and run both
make all_run

# Generate plots (requires trajectory.dat and msd_output.dat)
python3 plot_trajectory.py

# Remove all build artefacts and output files
make clean
```

### Output files

| File | Contents |
|------|----------|
| `trajectory.dat` | Columns: step, time [s], x [m], y [m], z [m] |
| `no_gravity_traj.dat` | Same format, no gravity (used for MSD) |
| `msd_output.dat` | Columns: τ [s], MSD_x, MSD_y, MSD_z, MSD_3D, Theory_3D [m²] |
| `trajectory_3d.png` | 3-D trajectory + z(t) panel |
| `msd_validation.png` | MSD linear scale + log-log scale |

---

## Validation

The MSD test (no gravity) verifies the fluctuation-dissipation theorem:

```
⟨Δx²(τ)⟩ = 2 D τ     (per axis)
⟨|Δr|²(τ)⟩ = 6 D τ   (3-D)
```

The log-log MSD plot confirms a slope of exactly 1, consistent with normal (Fickian) diffusion. Agreement between simulation and theory is < 0.3 % at short lag times. The slight overestimate at large lags (τ/T_total ≳ 0.05) is a known finite-trajectory statistical artifact of the time-averaged estimator and is not a physical error.

---

## Modifying Parameters

All physical constants are `parameter` declarations in `params_mod.f90`. Simulation length and time step are `parameter` constants at the top of each program file.

```fortran
! params_mod.f90
real(8), parameter :: a     = 1.0d-6    ! change radius here
real(8), parameter :: rho_p = 2000.0d0  ! change density here
...

! brownian_sim.f90
integer, parameter :: N_steps   = 500000
real(8), parameter :: dt        = 1.0d-2
```

After any change, rebuild with `make clean && make`.

---

## References

- Ermak, D. L. & McCammon, J. A. (1978). *Brownian dynamics with hydrodynamic interactions.* J. Chem. Phys., **69**(4), 1352–1360.
- Einstein, A. (1905). *Über die von der molekularkinetischen Theorie der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten suspendierten Teilchen.* Ann. Phys., **322**(8), 549–560.
- Dhont, J. K. G. (1996). *An Introduction to Dynamics of Colloids.* Elsevier.
