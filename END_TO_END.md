# Complete Python-led calculation workflow

This project implements the stages needed for an independent reproduction study
of [Thapa et al., PRL 128, 236402](https://doi.org/10.1103/PhysRevLett.128.236402).
It runs real calculations. The small examples verify software execution; they do
not reproduce the paper's formation times, electronic properties or conductivity
ratio. The original paper uses VASP/PAW; this freely available implementation
uses periodic all-electron PySCF with Gaussian orbitals. That substitution requires
convergence and cross-validation.

## Run the workflow

From the repository root:

```sh
.venv/bin/python end_to_end.py --config configs/end_to_end_smoke.json --out results/my_complete_test
```

Use a new output directory. To continue after an interruption:

```sh
.venv/bin/python end_to_end.py --config configs/end_to_end_smoke.json --out results/my_complete_test --resume
```

Completed stages are reused. An incomplete stage is preserved under
`failed_attempts/` and restarted. In particular, interrupted MD restarts from its
seed; the full Nose-Hoover thermostat state is not checkpointed. This is stage
recovery, not a bitwise continuation of an interrupted trajectory. The configuration
must be unchanged. Each new invocation saves source copies/hashes; completed
stages retain their original provenance when code changes between invocations.

Installation on another compatible machine:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

The GAP model and companion sparse files are external dependencies and must also be present. The README records
the download, checksum-verifying extractor and model provenance. No proprietary
VASP executable or PAW files are needed. ASE drives calculations; compiled QUIP,
PySCF, NumPy/SciPy and DFT-D3 libraries perform expensive numerical operations.
The workflow and analysis are Python, not a pure-Python implementation of DFT.

## What happens, from input to report

1. Generate periodic random carbon at the chosen density and seed. Run actual
   constant-volume Nose-Hoover-chain MD with the archived 2017 carbon GAP. Save
   coordinates, momenta, thermodynamics and structural time series.
2. Relax the final configuration at fixed volume with GAP/FIRE. Record whether
   forces actually converged. Construct an independent AB graphite reference at
   the same requested density. An optional adjacent-vacancy reference is relaxed;
   a 5-8-5 reconstruction is not assumed or imposed.
3. Measure coordination, local planar geometry, reciprocal density modulation,
   layer spacing, periodic chordless rings and 3D RDF. Compute a 2D in-plane RDF
   only when the layer checks pass. Translation, rotation and periodic images
   are included. These are documented independent order parameters, not the
   authors' unavailable exact planarity statistic.
4. Solve periodic PBE with a Gaussian all-electron basis and uniform k-point
   mesh. Save energies, occupations, coefficients, chemical potential, SCF logs
   and settings. Refuse unconverged SCF results. Compare energy per atom using
   consistently calculated sample and reference energies.
5. Evaluate total charge density, four orbital probability densities near the
   chemical potential at the first k-point, and the density of states. Export
   numerical arrays, Gaussian cube files and plots. Orbital identity is not
   automatically labeled pi or pi-star. HSE06 independently supplies an optional
   second charge/orbital calculation.
6. Use PBE momentum matrix elements to calculate the conductivity tensor at
   multiple photon energies and broadenings. Calculate complex/signed and
   positive space-projected conductivity (SPC). Check its integral against the
   conductivity and compare real-space matrix elements with analytic integrals.
7. Compare graphite and sample; report isotropic ratios, and in-plane ratios
   only when both structures have identified layers. A numerically zero reference
   conductivity produces an undefined ratio, not a large misleading number.
8. Optionally optimize cell and atoms at zero external pressure with PBE+D3.
   Finite differences of the electronic free energy plus D3 provide consistent
   forces/stress. This is expensive but implemented with actual DFT calculations.
   The separate zero-pressure geometry provides a spacing measurement; the
   fixed-density electronic comparison remains at its original volume.
9. Save `REPORT.md`, `summary.json`, `comparison.json`, provenance, source
   snapshots, raw states and visualization files.

The small default enables PBE, HSE06, charge/orbitals and transport for eight-atom
sample/reference cells, but disables the costly cell optimization and vacancy
case. `zero_pressure.enabled` and `defect_reference` enable those implemented
branches. The electronic integration test separately exercises real PBE+D3 cell
optimization on two atoms.

## Transport definition and units

Inside `agworkflow/transport.py`, energies are Hartree and positions are Bohr.
For per-spin Fermi occupations f and Gaussian delta broadening eta, the weight is

```
g_ij(k,w) = 2*pi/(V*Nk) * (f_i-f_j)/(E_j-E_i)
            * Gaussian_eta(E_j-E_i-w) * e^2/(hbar*a0)
```

The factor 2 counts spin. The degenerate limit of the divided difference is
`f*(1-f)/sigma_FD`. All ordered pairs, including diagonal matrix elements, are
included. Conductivity in S/m is `sum g_ij Re[p_ij,a conj(p_ij,b)]`.
The frequency input is photon energy in eV, not frequency in hertz.

This divided-difference Greenwood expression agrees with the usual 1/omega
expression on the exact delta-function support. At finite Gaussian broadening,
it is an explicit regularization; the value at zero photon energy is a
finite-broadening estimate, not a measured scattering-limited DC conductivity.
Basis, k-mesh, system size, temperature, active-band window and broadening must be
converged. See the [Kubo-Greenwood formulation](https://kubo-greenwood-formula-bd0633.gitlab.io/theory/).

Following [Subedi, Prasai and Drabold's SPC construction](https://people.ohio.edu/drabold/pubs/238.pdf),
let `xi_ij,a(r)=psi_i*(r) [-i d_a psi_j(r)]`. The row sum per unit volume is

```
row(r) = (1/3) sum_kij,a g_ij xi_ij,a(r) conj(p_ij,a)
```

Its complex integral recovers the isotropic conductivity. The positive SPC field
is `abs(row(r))`; its integral generally does not equal the conductivity. This
factorization avoids constructing a grid-by-grid matrix. Uniform quadrature
errors are measured and exported, not corrected by artificial renormalization.

Canonical momentum is used only for the local, all-electron PBE Hamiltonian.
Transport from hybrid states or pseudopotentials is rejected because the required
nonlocal velocity terms are not implemented. HSE06 is used for charge/orbitals.
Its default screened exchange uses `erfc(w*r)/r = 1/r - erf(w*r)/r`, reusing
full-Coulomb Gaussian density fitting and PySCF's long-range Fourier integrator.
Set `electronic.screened_exchange` to `direct` to use the slower direct screened
integrals. Finite fitting errors still require convergence. The implementation
uses [PySCF periodic DFT](https://pyscf.org/user/pbc/scf.html).

## Results and inspection

- `sample_md/trajectory.extxyz`: actual trajectory, viewable in atomistic viewers.
- `*_structure/summary.json`: layer checks and periodic ring counts.
- `*_structure/rdf_in_plane.csv`: present only when justified by detected layers.
- `*_pbe/states.npz`, `cell.json`, `settings.json`: reusable electronic states.
- `*_pbe/density.cube`, `orbital_k0_band*.cube`: volumetric charge/orbitals.
- `*_pbe/conductivity.csv` and `.npz`: full tensor in the NPZ; diagonals and
  isotropic average in CSV, resolved by photon energy and broadening.
- `*_pbe/spc_signed.cube`, `spc_absolute.cube`, `spc.npz`: conductivity fields.
- `*_pbe/quality.json`: charge, momentum, overlap and SPC integration errors.
- `*_pbe/*.png`, `*_hse06/*.png`: computed density, orbital, DOS and transport plots.
- `*_zero_pressure/summary.json`: convergence, stress, density and layer spacing.
- `*_in_plane_comparison.npz`: tensor projections and ratios when layers exist.

The all-electron core density needs a sufficiently fine field grid. The SCF
Becke grid and the exported uniform field grid are different numerical controls.
`require_field_accuracy=true` makes the workflow stop on failed configured field
checks. A converged SCF alone does not establish a converged basis or field grid.

## Tests and larger studies

Fast analytic/geometry tests:

```sh
.venv/bin/python -m unittest -v test_structure.py test_workflow.py
```

Actual optional-stage integration tests (use fresh directories):

```sh
.venv/bin/python validate_electronic.py cell --out results/my_cell_test
.venv/bin/python validate_electronic.py kpoints --out results/my_kpoint_test
.venv/bin/python validate_electronic.py exchange --out results/my_exchange_test
```

The cell test deliberately runs one optimizer step; it verifies forces/stress and
the optimizer path, not a converged zero-pressure material. The k-point test uses
a 3×1×1 mesh, including complex Bloch states away from time-reversal-invariant
points. The exchange test compares both HSE06 evaluations on one periodic C2 cell.

Electronic convergence on one unchanged structure:

```sh
.venv/bin/python electronic_convergence.py --structure results/end_to_end_smoke/reference/final.extxyz --config configs/end_to_end_smoke.json --plan configs/electronic_convergence.json --out results/my_electronic_convergence
```

The example plan increases field grid, k mesh, then basis. Extend it until the
observables of interest stabilize; its final entry is not automatically adequate.
Broadenings already form a sweep within every electronic analysis. Alter the
smearing, active window, SCF grid and other parameters in the same JSON plan.

Prepare independent-seed jobs without launching them:

```sh
.venv/bin/python run_sweep.py --config configs/research_start.json --grid configs/seed_sweep.json --out results/my_seed_study
```

Add `--run` to execute them sequentially, or `--run --resume` to reuse completed
stages. The sweep JSON accepts dotted paths such as `structural.temperature_K` and
`structural.density_g_cm3`. Atom-count studies should also choose appropriate
`reference_repeat` values; a Cartesian parameter grid does not enforce paired
sample/reference sizes.

`configs/research_start.json` is an expensive 80-atom, 100-ps starting template,
with a larger basis, k mesh, vacancy reference and zero-pressure relaxation.
It is not a prevalidated production recipe. Larger 400/1000-atom studies and
multiple seeds remain necessary to examine finite-size effects. An earlier
400-atom GAP benchmark on this Mac extrapolated to roughly ten days per 100 ps;
all-electron DFT, especially numerical cell relaxation, adds substantial cost.
No large production run is launched automatically.

## What remains before claiming the paper is reproduced

The code is available, but the scientific study still needs long trajectories,
actual formation of disordered graphitic layers, independent seeds, converged
electronic settings, and comparison of the paper's observables. The archived GAP
variant's exact identity, original configurations and missing supplementary
settings remain unresolved. The precise original VASP/PAW result cannot be
promised from this alternative backend. These are research and validation steps,
not curves that should be supplied by fitting the desired answer.
