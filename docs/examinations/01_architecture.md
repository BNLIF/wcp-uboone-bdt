# 01 — Architecture and Data Flow

**wcp-uboone-bdt** is the post-reconstruction analysis framework for the
MicroBooNE Low-Energy Excess (LEE) search and related cross-section measurements.
It takes output ntuples from the Wire-Cell Pattern Recognition (WCP) reconstruction,
applies BDT-based selection, builds multi-channel systematic covariance matrices,
and performs frequentist fits (chi-squared, goodness-of-fit, Feldman-Cousins) for
the LEE signal strength. The framework is written in C++17 with ROOT and built via
the `wcb`/`waf` build system.

---

## Table of Contents

1. [Pipeline overview](#1-pipeline-overview)
2. [Component dependency diagram](#2-component-dependency-diagram)
3. [Configuration files](#3-configuration-files)
4. [Hardcoded paths baked into source](#4-hardcoded-paths-baked-into-source)
5. [Build system](#5-build-system)
6. [Source tree layout](#6-source-tree-layout)

---

## 1. Pipeline overview

The analysis proceeds in five stages:

### Stage 1 — Event-level filtering and BDT scoring
**App:** `bdt_convert.cxx`

Reads five TTrees from the WCP selection output (`wcpselection/` directory per file):
- `T_BDTvars` — ~hundreds of BDT tagger input variables (struct `TaggerInfo`)
- `T_eval` — reco/truth summary (struct `EvalInfo`)
- `T_PFeval` — particle-flow + kinematics (struct `PFevalInfo`)
- `T_KINEvars` — reconstructed kinematic variables (struct `KineInfo`)
- `T_pot` — proton-on-target counters (struct `POTInfo`)

Applies ~20 TMVA and XGBoost BDT scores, filters by a hardcoded good-runs list
(~565 lines of inline integers in the source), and writes a condensed "checkout"
TTree per input file.

### Stage 2 — Per-file histogram filling
**Apps:** `convert_checkout_hist.cxx` and `convert_checkout_hist_xs.cxx`

Reads the condensed checkout TTrees. For each configured analysis channel (defined
in `configurations/cov_input.txt`), evaluates the channel selection cuts
(`LEEana::get_cut_pass` from `cuts.h`) and applies analysis weights
(`LEEana::get_weight`). Fills per-channel histograms and writes one ROOT file per
input file.

The `_xs` variant additionally fills 2D truth-signal × reco histograms for the
cross-section unfolding response matrix.

### Stage 3 — Merging and combined histogram production
**Apps:** `merge_hist.cxx` and `merge_hist_xs.cxx`

Opens every per-file histogram ROOT file and sums them into a single combined
`merge.root`. Uses `LEEana::CovMatrix` (`src/master_cov_matrix.cxx`) to manage
the channel/file mapping, and `LEEana::Bayes` to compute per-bin Bayesian MC-stat
uncertainties for the MC-stat covariance seed. Output filename is hardcoded as
`"merge.root"` (`apps/merge_hist.cxx`).

### Stage 4 — Systematic covariance matrix construction
Five dedicated apps (thin drivers over `LEEana::CovMatrix`):

| App | Covariance type | Main generator called |
|-----|----------------|-----------------------|
| `det_cov_matrix.cxx` | Detector systematics | `gen_det_cov_matrix` (mcm_1.h) |
| `xf_cov_matrix.cxx` | Flux + cross-section reweight | `gen_xf_cov_matrix` (mcm_2.h) |
| `xs_cov_matrix.cxx` | XS unfolding response | `gen_xs_cov_matrix` |
| `stat_cov_matrix.cxx` | Data statistical | `gen_data_stat_cov_matrix` (mcm_data_stat.h) |
| `stat_pred_cov_matrix.cxx` | MC prediction statistical | `gen_pred_stat_cov_matrix` (mcm_pred_stat.h) |

Each writes one or more fractional or absolute covariance ROOT files consumed
by the fitter in Stage 5.

The detector covariance construction uses an intermediate step: `merge_det.cxx`
first pairs CV and detector-variation NTuples by `(run, event)` matching into
a combined file, which `det_cov_matrix.cxx` then reads.

### Stage 5 — Fitting, GoF, and Feldman-Cousins
**App:** `read_TLee_v20.cxx`

Instantiates `TLee` (`src/TLee.cxx`), reads the merged histogram file and all
covariance ROOT files via paths set in `inc/WCPLEEANA/Configure_Lee.h`, assembles
the full covariance matrix (flux/XS + detector + MC-stat + data-stat), and runs:
- Best-fit LEE signal strength (Minuit2 MIGRAD via `Minimization_Lee_strength_FullCov`)
- Goodness-of-fit with conditional constraint (Schwartz decomposition)
- Feldman-Cousins Δchi-squared confidence intervals

### Wiener-SVD unfolding path (cross-section only)
`src/WienerSVD.cxx` is invoked from `xs_cov_matrix.cxx` (and demonstrated in
`wiener_example.cxx`). It unfolds the measured spectrum using the response matrix
from the `_xs` histogram files. See `07_algorithms.md#wiener-svd-unfolding`.

### GP smoothing (within Stage 4)
`src/GPRegressor.cxx` + `src/GPKernel.cxx` are called inside `gen_det_cov_matrix`
(`src/mcm_1.h`) to smooth the bias vector of the detector covariance bootstrap
before the amplification step. They are not separate pipeline stages; they are
called internally. See `07_algorithms.md#gaussian-process-regression`.

---

## 2. Component dependency diagram

```
Raw WCP ntuples per file
  (wcpselection/{T_BDTvars, T_eval, T_PFeval, T_KINEvars, T_pot})
           │
           ▼
    ┌─────────────────────┐
    │   bdt_convert       │  ← hardcoded good_run_list_vec
    │   (applies BDT      │  ← TMVA/XGBoost weight files (external)
    │    scores + filter) │  ← inc/WCPLEEANA/{tagger,eval,kine,pfeval,pot}.h
    └─────────────────────┘
           │  condensed checkout TTrees (one per input file)
           ▼
  ┌────────────────────────────────────────────┐
  │   convert_checkout_hist[_xs]               │
  │   (fills per-channel histograms per file)  │  ← configurations/cov_input.txt
  │   uses: cuts.h, eval.h, kine.h, pfeval.h   │  ← configurations/cv_input.txt
  └────────────────────────────────────────────┘
           │  per-file histogram ROOT files
           ▼
  ┌────────────────────────────────────────────┐
  │   merge_hist[_xs]                          │
  │   (sums histograms, runs Bayes MC-stat)    │  → merge.root  (hardcoded name)
  │   uses: LEEana::CovMatrix, LEEana::Bayes   │
  └────────────────────────────────────────────┘
           │  merge.root + Bayes error graphs
           ├──────────────────────┐
           ▼                      ▼
  Raw checkout TTrees     ┌──────────────────────────────────────┐
  (CV + detvar pairs)     │ Systematic covariance drivers        │
           │              │                                      │
           ▼              │  xf_cov_matrix  → cov_{N}.root       │
  ┌────────────────┐      │  det_cov_matrix → cov_det_{N}.root   │
  │  merge_det     │      │  xs_cov_matrix  → cov_xs_{N}.root    │
  │  (pair CV +    │      │  stat_cov_matrix → cov_data_stat.root│
  │   detvar TTrees│      │  stat_pred_cov  → cov_mc_stat/       │
  │   by run,event)│      │                                      │
  └────────────────┘      │  All use: LEEana::CovMatrix          │
           │              │  det uses: GPSmoothing (GPRegressor) │
           └──────────────┴──────────────────────────────────────┘
                                    │  covariance ROOT files
                                    ▼
                      ┌─────────────────────────────┐
                      │   read_TLee_v20              │
                      │   (TLee class)               │  ← Configure_Lee.h (baked-in config)
                      │   - best-fit Lee strength    │  ← merge.root (spectra)
                      │   - GoF (conditional constr.)│  ← flux/XS cov files
                      │   - Feldman-Cousins Δchi²    │  ← detector cov files
                      │                             │  ← MC-stat log tables
                      └─────────────────────────────┘
                                    │
                               Fit results
                         (histograms, FC scan trees)
```

---

## 3. Configuration files

Runtime behavior is controlled by plain-text config files expected under
`./configurations/` relative to the working directory from which the apps are run.
None of these files exist in the repository; they are created per-analysis.

| File | Used by | Content |
|------|---------|---------|
| `configurations/cov_input.txt` | `CovMatrix` ctor (`master_cov_matrix.cxx:133`) | Channel list: name, variable, nbins, range, file-type flag, syst flags, MC-stat group, LEE flag |
| `configurations/cv_input.txt` | `CovMatrix` ctor (`master_cov_matrix.cxx:215`) | File-type → input ROOT file → output ROOT file, POT, run period |
| `configurations/file_ch.txt` | `CovMatrix` ctor (`master_cov_matrix.cxx:230`) | Input file → cut name mapping |
| `configurations/rw_cv_input.txt` | `CovMatrix` ctor (`master_cov_matrix.cxx:246`) | Reweight config: weight name, truth variable, binning |
| `configurations/xs_ch.txt` | XS mode | XS channel mapping |
| `configurations/xs_real_bin.txt` | XS mode | True-to-reco bin correspondence |
| `configurations/gp_input.txt` | GP smoothing (`mcm_1.h:152`) | GP training data and kernel config |
| `configurations/alt_var_xbins.txt` | `alt_var_index` (`cuts.h:4034`) | Alternative variable bin edges (lazy-loaded) |

**TLee configuration** is not a runtime file; it is a set of global variables
defined in `inc/WCPLEEANA/Configure_Lee.h` and compiled into `read_TLee_v20`.
This includes paths to the spectra file, flux/XS covariance directory, detector
covariance directory, and MC-stat log directory. Changing the input dataset
requires recompiling. See `05_bugs.md#B-05` for the ODR/guard concern.

**dict/LinkDef.h** — ROOT dictionary generation descriptor for the package.
Contains ROOT `ClassDef` linkage entries. Not a runtime config.

---

## 4. Hardcoded paths baked into source

Several absolute filesystem paths are compiled into the binaries:

| Location | Hardcoded path | Notes |
|----------|---------------|-------|
| `apps/applyNuMIGeomtryWeights.cxx:34–37` | `/home/xqian/wire-cell/wcp-uboone-bdt/scripts/NuMI_Geometry_Weights_Histograms.root` | Default argument; overrideable at runtime |
| `apps/applyNuMIGeomtryWeights.cxx:35` | `/data1/xqian/...` | Second hardcoded example path |
| `src/WienerSVD.cxx:107–113` | `/uboone/data/users/lcoopert/...` | Inside a **commented-out** TFile write block; not active |
| `inc/WCPLEEANA/Configure_Lee.h:5–59` | Multiple `/home/xji/...`, `/home/lee/...`, `/home/xqian/...` paths | All **commented-out** alternative configurations |
| `apps/merge_hist.cxx` | `"merge.root"` | Hardcoded output filename (relative path) |

No currently-active path in the compiled code references a machine-specific
absolute path. The `applyNuMIGeomtryWeights` paths are default argument values
and can be overridden at invocation.

---

## 5. Build system

**File:** `wscript` (31 lines, Python/waf)

- `bld.smplpkg('WCPuBooNE_BDT_APP', use=['ROOTSYS'])` — the `smplpkg` helper
  provided by `wcb` globs all `*.cxx` files in the **flat** `apps/` directory.
  Files in `apps/old/` are **not** included (they are in a subdirectory). There
  is no explicit app enumeration; adding a new `.cxx` to `apps/` automatically
  builds a new executable.
- `HAVE_VLNEVAL` — if the VLNEVAL library is present (external), Boost
  `program_options` is also linked, enabling `eval_vlne.cxx` which wraps the
  variable-length-neutrino energy estimator.
- Compiler flags: `-Wall -Wno-unused-local-typedefs -Wno-unused-function`.
  `-Wpedantic` and `-Werror` are commented out (`wscript:19`).
- ROOT is the only required external library in a standard build.

---

## 6. Source tree layout

```
wcp-uboone-bdt/
├── wscript                     waf build definition
├── dict/
│   └── LinkDef.h               ROOT dictionary descriptor
├── inc/WCPLEEANA/              Header-only event-record + selection layer
│   ├── Configure_Lee.h         TLee fitter runtime settings (global variables)
│   ├── cuts.h                  Selection engine, weights, signal predicates (4067 L)
│   ├── tagger.h                BDT tagger variables struct (2994 L)
│   ├── bdt.h                   cal_*_bdt implementations (650 L)
│   ├── eval.h, kine.h,         Event record structs
│   │   pfeval.h, weights.h,
│   │   pot.h
│   ├── master_cov_matrix.h     CovMatrix class declaration
│   ├── TLee.h                  TLee class declaration
│   ├── WienerSVD.h             WienerSVD function declaration
│   ├── bayes.h                 Bayes class declaration
│   ├── GPKernel.h, GPRegressor.h, GPPoint.h
│   └── Util.h
├── src/                        Class implementations
│   ├── TLee.cxx                Fit / GoF / Feldman-Cousins (2848 L)
│   ├── master_cov_matrix.cxx   Systematic covariance builder (2285 L)
│   ├── mcm_1.h                 Detector cov + GP smoothing (helper, #included)
│   ├── mcm_2.h                 Flux/XS cov (helper, #included)
│   ├── mcm_data_stat.h         Data-stat cov (helper, #included)
│   ├── mcm_pred_stat.h         MC-stat cov (helper, #included)
│   ├── WienerSVD.cxx           Wiener-SVD unfolding (236 L)
│   ├── WienerSVD_3D.C          3D stitched regularization matrices (#included)
│   ├── bayes.cxx               Bayesian MC-stat posterior (415 L)
│   ├── GPKernel.cxx            RBF + RQ kernel functions (141 L)
│   ├── GPRegressor.cxx         GP posterior inference (202 L)
│   ├── GPSmoothing.C           GP smoothing driver (#included in mcm_1.h)
│   ├── Util.cxx                Matrix/histogram utility functions (152 L)
│   └── draw.icc                ROOT style helpers (#included in TLee.cxx)
├── apps/                       ~30 executables (all built by wscript glob)
│   ├── read_TLee_v20.cxx       LEE fit / FC driver
│   ├── bdt_convert.cxx         Event-level converter + BDT scorer (3242 L)
│   ├── merge_hist[_xs].cxx     Histogram merger
│   ├── merge_det.cxx           CV + detvar TTree pairer
│   ├── {det,xf,xs,stat,        Covariance matrix producers
│   │   stat_pred}_cov_matrix
│   ├── convert_checkout_hist   Per-file histogram filler
│   │   [_xs].cxx
│   └── old/                    Dead code (not built)
│       ├── nueCC_convert.cxx
│       └── numuCC_convert.cxx
├── scripts/                    Standalone ROOT macros (not built)
│   ├── plot_FC_new.cc          FC plot display
│   ├── plot_systematics.cc     Systematics breakdown plots
│   ├── extractNuMIGeometry     NuMI geometry weight builder
│   │   Weights.C
│   └── NuMI_Geometry_Weights   Prebuilt output of above
│       _Histograms.root
└── docs/
    └── examinations/           This examination (read-only documentation)
```
