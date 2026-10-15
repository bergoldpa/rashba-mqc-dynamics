# Mixed quantum-classical dynamics: Ehrenfest, koopmons and exact (SOFT) reference

All MATLAB code lives in `core/`. The five test cases differ only in their
parameters (`caseN/case_parameters.m`); their results are stored in the case folders.

## Authors and citation

Paul Bergold, Giovanni Manfredi, Cesare Tronci

This code accompanies the paper

> P. Bergold, G. Manfredi, C. Tronci, *Quantum-classical dynamics of Rashba spin-orbit
> coupling*, arXiv:2603.23758 (2026), https://arxiv.org/abs/2603.23758

If you use this code, please cite the paper and the code (see `CITATION.cff`).

## License

MIT License, see `LICENSE`.

## Results

The repository contains only the code and the LaTeX sources of the figures. All
results (`.mat`, `.dat`, `.png`, `.fig`, compiled figures) are created by running
`run_project.m` (see below) and are written into the case folders.

## Folder structure

```
<project folder>/
├── setup.m                                 adds core/ to the MATLAB path
├── run_project.m                           driver script: comment lines in/out, then Run
├── core/                                   shared code (only one copy)
│   ├── get_case_parameters.m               reads case_parameters.m of the current case
│   ├── run_case.m                          postprocessing of one case (all or selected parts)
│   ├── run_simulation.m                    runs one simulation of one case in its folder
│   ├── simulation/
│   │   ├── ehrenfest_simulation.m          Ehrenfest (mean-field) dynamics
│   │   ├── ehrenfest_input.m               its input parameters
│   │   ├── koopmons_simulation.m           koopmon dynamics (with correlation terms)
│   │   ├── koopmons_input.m                its input parameters
│   │   ├── soft_simulation.m               exact quantum reference (split-operator FFT)
│   │   └── max_energy_drift.m              energy conservation check of a finished run
│   └── postprocessing/
│       ├── plot_phase_space_soft.m         Wigner function snapshots (run first)
│       ├── plot_phase_space_ehrenfest_koopmons.m
│       ├── phase_space_plot_settings.m     shared settings of both phase-space scripts
│       ├── plot_bloch_purity.m / plot_bloch_purity_soft.m
│       ├── compute_spin_correlation.m / compute_spin_correlation_soft.m
│       └── convert_bloch_purity_time.m / convert_spin_correlation_time.m
├── case1 … case5/
│   ├── case_parameters.m                   all parameters of this case (model parameters in atomic
│   │                                       units, physical values in the comments; numerical settings)
│   ├── soft/                               soft_simulation.mat
│   ├── ehrenfest/n500/                     ehrenfest_simulation.mat
│   ├── koopmons/n500_a05/                  koopmons_simulation.mat
│   ├── postprocessing/                     working folders of the postprocessing scripts
│   │   └── phase_space/ bloch_purity/ spin_correlation/
│   └── figures/                            LaTeX (pgfplots) figures, data in input/
│       └── phase_space/ bloch_purity/ spin_correlation/
└── potential_energy_surfaces/              energy surfaces of the MomC model vs momentum
    ├── momentum_surfaces.m                 computes the surfaces of all cases -> figures_momentum/MOM*.dat
    └── figures_momentum/                   LaTeX figures: cases 1–2 and 3–5 side by side
```

Each run, postprocessing and figure-input folder contains a short `README.md` describing
what runs there and which files it produces.

Case 3 additionally contains `ehrenfest/n5000/` (N = 5000, set in its
`local_parameters.m`) and `ehrenfest/n500_vs_n5000/compare_n500_vs_n5000.m`.

## Requirements

- MATLAB R2022a or newer (tested with R2023b and R2026a)
- Statistics and Machine Learning Toolbox (`sobolset`, `norminv` for the initial sampling
  of the Ehrenfest and koopmon simulations)
- For the figures: a LaTeX distribution with `pgfplots` and the `standalone` class

## Setup

In MATLAB, open the project folder and run

```matlab
setup       % adds core/ to the MATLAB path
savepath    % optional: keep the path for future sessions
```

## Reproducing the results

**Easiest:** open `run_project.m`, set `cases = [...]`, comment in the lines you want
(simulations, postprocessing parts, `momentum_surfaces`) and press **Run**.
Each simulation runs in its own folder; afterwards MATLAB returns to the project folder.

**By hand,** for each test case N = 1, …, 5 (in MATLAB, after `setup`):

1. Run the simulations, each in its own folder:
   ```matlab
   cd caseN/soft;              soft_simulation
   cd ../ehrenfest/n500;       ehrenfest_simulation
   cd ../../koopmons/n500_a05; koopmons_simulation
   ```
   For case 3, also run `ehrenfest_simulation` in `case3/ehrenfest/n5000`, then
   `compare_n500_vs_n5000` in `case3/ehrenfest/n500_vs_n5000`.
2. Run the complete postprocessing: `run_case(N)`.
3. Compile the LaTeX figures in `caseN/figures/phase_space`, `…/bloch_purity` and
   `…/spin_correlation` (e.g. `pdflatex phase_space_summary.tex`).

Energy surfaces: run `momentum_surfaces` (any folder), then compile the two `.tex`
files in `potential_energy_surfaces/figures_momentum`.

## Running a simulation

Change to the run folder and call the script, e.g.

```matlab
cd case2/ehrenfest/n500
ehrenfest_simulation
```

`get_case_parameters` searches upwards from the current folder for `case_parameters.m`.
Results (`<script name>.mat`, `plot_energy.fig`, `time_seconds.dat`) are written to the current folder.

- **Override a value for one folder only:** add a `local_parameters.m`, e.g. `P.n_b = 5000;`.
- **New test case:** create a folder, copy a `case_parameters.m` into it and adapt the values.
- **Harmonic potential:** `P.harmonic = 1` (on) or `0` (off).
- **Time unit of exported data:** `P.time_scale = "PICO"` (ps) or `"NANO"` (ns).

## Postprocessing

**Everything at once:** `run_case(n)` runs the complete postprocessing of case n
(phase space, Bloch/purity, spin correlation, time conversion) from the existing
simulation results and copies the input files of the LaTeX figures into
`caseN/figures/*/input/`. Then compile the `.tex` files.

**Selected parts:** `run_case(n, "phase_space")`, `run_case(n, "bloch_purity")`,
`run_case(n, "spin_correlation")`, or several, e.g. `run_case(n, ["phase_space" "bloch_purity"])`.

**Simulations from any folder:** `run_simulation(n, method)` with method
`"soft"`, `"ehrenfest"`, `"ehrenfest_n5000"` (case 3 only) or `"koopmons"`.

**Single scripts:**


1. Change to the matching folder `caseN/postprocessing/<topic>` and run the script.
   The scripts read the simulation results directly from the case folders
   (e.g. `../../soft/soft_simulation.mat`) and never modify or delete them.
   In the Ehrenfest/koopmon scripts, choose the result with the `string_name` line.
2. The scripts don't ask questions. The method (`K` = koopmons, `E` = Ehrenfest) follows from the
   loaded file; the smoothing width `alpha_opt` (0.2) is set at the top of
   `plot_phase_space_ehrenfest_koopmons.m`.
3. For the time-series plots, run `convert_*_time.m` in the `output/` folder
   (ps or ns, set by `P.time_scale`). Then copy the `*_N.dat` files to
   `caseN/figures/…/input/` and compile the `.tex` file.

For the phase-space plots, run `plot_phase_space_soft.m` first. It creates
`WIGNER_VAL_REF.mat`, which defines the colour scale for the Ehrenfest/koopmon plots.

## Conventions

- SOFT stores 1001 snapshots at t = 0, τ, …, t_final; Ehrenfest/koopmons store the
  states at t = 0, dt, 2·dt, … In all output, snapshot l is at time (l−1)·step.
- If an Ehrenfest/koopmon run stops early (relative energy error > `tol_energy`),
  only the computed steps are saved and a warning is shown.
- Atomic time units are converted to ps/ns with 1 a.u. of time = 2.4188843265857e-17 s (CODATA).
- Output files of the postprocessing start with the method: `soft_`, `ehrenfest_`, `koopmons_`
  (e.g. `ehrenfest_n1.dat`, `koopmons_bloch_full.fig`), so the methods never overwrite each other.
- All `.dat` files are written with `writematrix` at full precision.
- Only the MomC model is implemented in the simulations.
- Spin correlation coefficients are NaN where a standard deviation vanishes
  (e.g. at t = 0 with the spin along z); these points are left out of the `_N.dat` files.

## Stored data

- `caseN/soft/soft_simulation.mat` holds the wave functions `psi1_val`, `psi2_val`.
  The densities `abs(psi1_val).^2`, `abs(psi2_val).^2` and their sum are not stored;
  they can be recomputed exactly from these.
- `caseN/postprocessing/phase_space/WIGNER_VAL_REF.mat` holds only the colour scale and
  settings (a few MB). To get the Wigner functions back exactly, rerun `plot_phase_space_soft.m`
  on `soft_simulation.mat`; they are then in the workspace as `wigner_fft_val`.
