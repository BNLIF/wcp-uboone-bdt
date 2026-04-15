# 02 — Core Framework Technical Examination

*MicroBooNE Low-Energy-Excess Analysis — wcp-uboone-bdt*
*Examined: 2026-04-14*

---

## Table of Contents

1. [Util — Matrix / Histogram Utilities](#1-util--matrix--histogram-utilities)
2. [WienerSVD — Spectral Unfolding](#2-wienersvd--spectral-unfolding)
3. [GPKernel / GPRegressor / GPSmoothing — Gaussian-Process Smoothing](#3-gpkernel--gpregressor--gpsmoothing--gaussian-process-smoothing)
4. [bayes — Bayesian Error Propagation](#4-bayes--bayesian-error-propagation)
5. [master_cov_matrix (CovMatrix) — Covariance Matrix Construction](#5-master_cov_matrix-covmatrix--covariance-matrix-construction)
6. [TLee — Top-Level Statistical Engine](#6-tlee--top-level-statistical-engine)

---

## 1. Util — Matrix / Histogram Utilities

### Purpose

`Util.cxx` / `inc/WCPLEEANA/Util.h` provide a small set of interactive and converter utilities used during development and debugging. They convert between ROOT histogram objects and the `TMatrixD`/`TVectorD` objects used throughout the analysis, and contain diagnostic print routines for SVD decomposition verification.

### Key Classes and Methods

| Function | File : Line | Role |
|---|---|---|
| `Matrix(row, col)` | `Util.cxx:12` | Interactively fills and returns a `TMatrixD` from `cin` |
| `Vector(row)` | `Util.cxx:40` | Interactively fills and returns a `TVectorD` from `cin` |
| `MatrixMatirx(M1, M2)` | `Util.cxx:29` | Prints M1, M2, and their product M1×M2 |
| `MatrixVector(M, V)` | `Util.cxx:55` | Prints M, V, and M×V |
| `SVD(M)` | `Util.cxx:67` | Performs ROOT `TDecompSVD`, prints U, V, D, and the reconstruction M≈UDV |
| `H2M(histo, mat, rowcolumn)` | `Util.cxx:109` | Copies a `TH2D` into a `TMatrixD`; `rowcolumn` selects transposition convention |
| `H2V(histo, vec)` | `Util.cxx:123` | Copies a `TH1D` into a `TVectorD` |
| `M2H(mat, histo)` | `Util.cxx:132` | Copies a `TMatrixD` into a `TH2D` |
| `V2H(vec, histo)` | `Util.cxx:144` | Copies a `TVectorD` into a `TH1D` |

### Inputs Consumed / Outputs Produced

- Inputs: `TH1D*`, `TH2D*`, `TMatrixD`, `TVectorD` — all provided by the caller in-memory.
- Outputs: modified `TMatrixD`/`TVectorD` passed by reference, or printed diagnostics to `stdout`.
- No ROOT files are opened or written.

### Algorithm Sketch

`H2M` iterates histogram bins and assigns `mat(i,j) = histo(i+1,j+1)` (or transposed when `rowcolumn=false`). `SVD` uses `TDecompSVD` and reconstructs M from the U, D, V factors to verify factorization quality. The interactive `Matrix` and `Vector` routines read values from `stdin` one entry at a time, making them unsuitable for automated pipelines.

### Notes and Red Flags

- **`Util.cxx:29`** — The function is named `MatrixMatirx` (typo: "Matirx"), but the header `Util.h:17` declares it as `MatrixMatrix`. The definition name and the declaration name disagree; calling `MatrixMatrix` as declared will produce a linker error or resolve to the wrong symbol depending on how the header is consumed. → see 05_bugs.md#B-01
- **`Util.cxx:12`** — `Matrix()` and `Vector()` (`Util.cxx:40`) read from `cin` inside a library function. Any code path that calls these during a batch job will block indefinitely waiting for terminal input.
- **`Util.cxx:105`** — `SVD()` calls `MM.Draw("colz")` at the end; calling `Draw` outside a `TCanvas` context is a no-op at best, an error at worst, and indicates leftover debug code.
- **`Util.h:17`** — The header declares `MatrixMatrix` but the implementation spelling is `MatrixMatirx`. This is a compilation / linkage mismatch that will silently fail or be unreachable.

---

## 2. WienerSVD — Spectral Unfolding

### Purpose

`WienerSVD.cxx` / `inc/WCPLEEANA/WienerSVD.h` implement the Wiener-filter variant of SVD-based unfolding described in arXiv:1710.09757. Given a detector response matrix, a measured spectrum, a signal prior, and a covariance matrix, it returns the unfolded signal spectrum, the additional-smearing matrix (`AddSmear`), the Wiener-filter weights (`WF`), and the full covariance of the unfolded result. The companion file `WienerSVD_3D.C` provides the `C3_3D()` helper that builds the 3D smoothness regularisation matrix.

### Key Classes and Methods

| Function | File : Line | Role |
|---|---|---|
| `Matrix_C(n, type)` | `WienerSVD.cxx:29` | Constructs the n×n smoothness matrix C for types 0 (unit), 1 (1st deriv.), 2 (2nd deriv.), 22/23/32/33 (multi-dim variants) |
| `WienerSVD(...)` | `WienerSVD.cxx:119` | Main unfolding function; returns unfolded `TVectorD` and fills `AddSmear`, `WF`, `UnfoldCov`, `covRotation`, `covRotation_t` |
| `is_on_edge(edges, element)` | `WienerSVD.cxx:16` | Helper: true if `element` is in `edges` vector |
| `get_slice(edges, index)` | `WienerSVD.cxx:22` | Helper: returns the slice index for a bin |
| `C3_3D(version)` | `WienerSVD_3D.C` | Builds third-derivative regularisation matrix for 3D unfolding |

### Inputs Consumed / Outputs Produced

- **Inputs (all in-memory):** `Response` (m×n TMatrixD), `Signal` (n-element TVectorD prior), `Measure` (m-element TVectorD), `Covariance` (m×m TMatrixD), `C_type` (int), `Norm_type` (float), `flag_WienerFilter` (float, default 1.0).
- **Outputs (by reference):** `AddSmear` (n×n additional-smearing matrix), `WF` (n-element Wiener-filter weights), `UnfoldCov` (n×n covariance of unfolded result), `covRotation`, `covRotation_t`.
- No ROOT files are opened or written by this function itself.

### Algorithm Sketch

The covariance matrix is decomposed via `TDecompSVD` to obtain an orthogonal rotation Q that whitens the measurement space. The rotated response R = Q·Response is then factored as U·D·Vᵀ. The smoothness matrix C (chosen by `C_type`) is applied to the signal prior and normalised by `Norm_type`. The Wiener filter W is built bin-by-bin as W(i,i) = S(i)²/(D(i)²·S(i)²+1), weighting each SVD mode by its signal-to-noise ratio. The unfolded spectrum is C⁻¹·V·W·Dᵀ·Uᵀ·Q·Measure, and the covariance propagates via covRotation·Covariance·covRotationᵀ. → see 07_algorithms.md#wiener-svd-unfolding

### Notes and Red Flags

- **`WienerSVD.cxx:190`** — `bool test_wiener = false;` is hardcoded. The alternative Wiener filter formula at line 209 (using `signal_ratio=2.0` instead of the prior) is permanently disabled. This is dead code that cannot be activated without source changes. → see 05_bugs.md#B-02
- **`WienerSVD.cxx:40`** — For `type==22` or `type==32`, `dim_edges` is hardcoded as `{ 0, 3, 7, 11, 14, 18, 22, 26, 31, 36}`. For `type==23` or `type==33`, a different hardcoded set is used at line 47–50. These bin-boundary tables are not derived from the runtime dimensions `n`; if the binning changes, the tables must be updated manually.
- **`WienerSVD.cxx:52`** — An unconditional `std::cout` debug print executes on every call, regardless of type. In batch production this produces verbose output on every unfolding call.
- **`WienerSVD.cxx:107–113`** — A commented-out block writes a debug ROOT file to a hardcoded `/uboone/data/users/lcoopert/...` path. This is developer-specific and would silently fail on any other system.
- **`WienerSVD.cxx:163`** — `C0.Invert()` is called on `C0` in-place (making it `C_inv`), then `Signal = C*Signal` uses the original `C`. This relies on `C` being a copy of `C0` made at line 162 (`TMatrixD C = C0`), which is correct, but the flow is subtle and easy to break if `C0` is reused.

---

## 3. GPKernel / GPRegressor / GPSmoothing — Gaussian-Process Smoothing

### Purpose

These files implement a Gaussian Process regressor used to smooth the detector-systematic covariance matrix after bootstrapping. `GPKernel.cxx` / `GPKernel.h` define the kernel functions (Radial Basis Function and Rational Quadratic). `GPRegressor.cxx` / `GPRegressor.h` implement fitting via marginal-likelihood maximisation and posterior prediction. `GPPoint.h` defines the 5-dimensional input point. `GPSmoothing.C` is the top-level driver called from `mcm_1.h`.

### Key Classes and Methods

| Class / Function | File : Line | Role |
|---|---|---|
| `GPPoint(double x[5])` | `GPPoint.h:9` | 5-dimensional input point; stores coordinates in `double x[5]` |
| `GPKernel::operator()(pts)` | `GPKernel.cxx:7` | Builds symmetric kernel matrix K(pts, pts) |
| `GPKernel::KernelBlock(pts1, pts2)` | `GPKernel.cxx:21` | Builds rectangular cross-kernel block K(pts1, pts2) |
| `RBFKernel::Mag(p1, p2)` | `GPKernel.cxx:43` | Computes scaled squared distance used in the RBF kernel |
| `RBFKernel::Element(pt1, pt2, dpar_idx)` | `GPKernel.cxx:66` | Evaluates RBF kernel element or its derivative w.r.t. hyperparameter dpar_idx |
| `RBFKernel::SetThetas(thetas)` | `GPKernel.cxx:93` | Sets hyperparameters from log-space optimiser output via `exp(theta[i])` |
| `GPRegressor::Fit(X, y, solveHP)` | `GPRegressor.cxx:14` | Normalises y, optionally solves hyperparameters, computes Cholesky K, fAlpha, fKInv |
| `GPRegressor::Predict(X)` | `GPRegressor.cxx:94` | Computes posterior mean and covariance at new points |
| `GPRegressor::SolveHyperParameters()` | `GPRegressor.cxx:65` | Maximises marginal likelihood using ROOT BFGS2 minimiser |
| `MarginalLikelihood::DoEval(par)` | `GPRegressor.cxx:178` | Evaluates negative log marginal likelihood |
| `MarginalLikelihood::DoDerivative(par, ipar)` | `GPRegressor.cxx:191` | Computes gradient of negative log marginal likelihood |
| `GPSmoothing(vec_mean, cov_mat, filename, flag)` | `GPSmoothing.C:50` | Driver: reads config, constructs GP, replaces bootstrapped mean and covariance with GP posterior |

### Inputs Consumed / Outputs Produced

- **GPSmoothing inputs:** `vec_mean` (TVectorD*, bootstrapped mean differences), `cov_mat_bootstrapping` (TMatrixD*, bootstrapped covariance), `./configurations/gp_input.txt` (tab-separated: log_scales, 6 kernel params, 5 lines of bin centres).
- **GPSmoothing outputs:** overwrites `*vec_mean` and `*cov_mat_bootstrapping` in-place with the GP posterior mean and covariance (when `smoothing_par % 2 == 1`).
- Debug output written to `./hist_rootfiles/DetVar/debug_smoothing.root` when `smoothing_par >= 2`.

### Algorithm Sketch

`GPSmoothing` reads a configuration file specifying the 5D bin-centre grid and 6 RBF hyperparameters, then constructs one `GPPoint` per bin. `GPRegressor::Fit` normalises the target vector, builds the Cholesky decomposition of K+noise, and computes `fAlpha = K⁻¹y`. `GPRegressor::Predict` evaluates the posterior mean as K(X*,X)·fAlpha and posterior covariance as K(X*,X*) − K(X*,X)·K⁻¹·K(X,X*). The smoothed posterior covariance replaces the raw bootstrapped covariance, reducing statistical noise in the detector-systematic matrix. → see 07_algorithms.md#gp-smoothing

### Notes and Red Flags

- **`GPKernel.cxx:50–52`** — `RBFKernel::Mag` receives `GPPoint` arguments by value, then directly mutates the raw coordinate pointers obtained via `p1x = p1.X()` and writes `p1x[i] = log(p1x[i])` (and likewise `p2x[i]`). Because `GPPoint::X()` returns a pointer to the internal `double x[5]` array, this permanently modifies the point's coordinates every time `Mag` is called in log-scale mode. Any subsequent kernel evaluation using the same `GPPoint` objects will use corrupted coordinates. → see 05_bugs.md#B-03
- **`GPRegressor.cxx:191–200`** — `MarginalLikelihood::DoDerivative` uses `fAlpha` and `fKInv`, which are populated by `MarginalLikelihood::Solve`. If the BFGS minimiser calls `DoDerivative` with a `par` that differs from the last `DoEval` call, `fAlpha` and `fKInv` are stale (computed at the previous parameter point). The derivative is then evaluated at an inconsistent point. → see 05_bugs.md#B-04
- **`GPRegressor.cxx:71–80`** — `SolveHyperParameters` allocates `double* par = new double[n]` at line 71 but never frees it. This is a heap memory leak per call to `SolveHyperParameters`. → see 05_bugs.md#B-05
- **`GPRegressor.cxx:79`** — All hyperparameter bounds are set to `(-5, 5)` in log-space (i.e., `exp(-5)` to `exp(5)` ≈ 0.0067 to 148). There is no mechanism to widen these bounds, which may prevent the optimiser from finding solutions for very large or very small length scales.
- **`GPSmoothing.C:114`** — `reg.Fit(gp_points, (*vec_mean), false)` passes `false` for `solveHyperParams`, meaning hyperparameters from the config file are used without optimisation. This is by design but means the kernel parameters must be manually tuned.

---

## 4. bayes — Bayesian Error Propagation

### Purpose

`bayes.cxx` / `inc/WCPLEEANA/bayes.h` implement a `Bayes` class that computes the Bayesian posterior PDF for a Poisson-distributed signal rate by convolving multiple "measurement components". Each component represents a sample (e.g. from a different run period) with its own expected count and uncertainty. The class is used by `CovMatrix` to compute per-bin Bayesian uncertainties (via `get_bayes_errors`).

### Key Classes and Methods

| Class / Function | File : Line | Role |
|---|---|---|
| `Bayes::Bayes()` | `bayes.cxx:8` | Constructor; zeroes `num_component`, `mean`, `acc_sigma2`; nulls `f_conv`, `f_conv_num`, `g1` |
| `Bayes::~Bayes()` | `bayes.cxx:20` | Deletes all `TF1*` in `meas_pdf_vec`, `f1_num_vec`, `f_test_vec`; deletes `f_conv_num`, `f_conv`, `g1` |
| `Bayes::add_meas_component(meas, sigma2, weight, flag)` | `bayes.cxx:42` | Adds one measurement component; creates a test TF1 for MC sampling |
| `Bayes::do_convolution()` | `bayes.cxx:104` | Performs sequential FFT convolution of all component PDFs; samples result into TGraph `g1` |
| `Bayes::get_covariance()` | `bayes.cxx:294` | Numerical integral over `f_conv_num` to compute variance around `mean` |
| `Bayes::get_covariance_mc()` | `bayes.cxx:333` | MC sampling (100 000 throws) for cross-check of `get_covariance()` |
| `Bayes::calculate_lower_upper(nsigma)` | `bayes.cxx:183` | Bisection search for the lower and upper bounds enclosing `nsigma` credible interval |
| `Prop_Poisson_Pdf(x, par)` | `bayes.cxx:377` | TF1 callback implementing the effective Poisson PDF for one component |

### Inputs Consumed / Outputs Produced

- **Inputs:** per-component `(meas, sigma2, weight)` tuples added via `add_meas_component`.
- **Outputs:** `f_conv` (pointer to convolved TF1), `f_conv_num` (TGraph-backed TF1 for integration), `mean` (sum of component means), `get_covariance()` / `calculate_lower_upper()` return values.
- No ROOT files opened or written.

### Algorithm Sketch

Each call to `add_meas_component` accumulates a component mean and variance. `do_convolution` iterates over components: the first component directly initialises `f_conv`; each subsequent component is convolved with the running `f_conv` using `TF1Convolution` (FFT, 10 000 points). The final `f_conv` is sampled on a 5001-point grid into `TGraph g1`, which is then wrapped in `f_conv_num` for fast integration. Credible intervals are found by bisection on the cumulative integral of `f_conv_num`. → see 07_algorithms.md#bayesian-credible-intervals

### Notes and Red Flags

- **`bayes.cxx:138–151`** — Inside `do_convolution`, each new convolved `TF1` is appended to `f_conv_vec` at line 149. The destructor (line 30–33) explicitly comments out the deletion loop for `f_conv_vec`. All intermediate convolution `TF1*` objects are therefore leaked when `Bayes` is destroyed. → see 05_bugs.md#B-06
- **`bayes.cxx:37–38`** — The destructor checks `if (f_conv != (TF1*)0) delete f_conv`. However, `f_conv` is also the last element of `f_conv_vec` (line 151: `f_conv = new_f_conv`; line 149: `f_conv_vec.push_back(new_f_conv)`). If the `f_conv_vec` deletion loop were ever re-enabled, `f_conv` would be deleted twice. → see 05_bugs.md#B-07
- **`bayes.cxx:46–49`** — The `weight` argument to `add_meas_component` is ignored for non-zero `meas`: the stored weight is always 1 when `meas != 0`. The `flag` parameter is accepted but unused entirely. API callers may not realise their weight is discarded.
- **`bayes.cxx:342–345`** — `get_covariance_mc` draws 100 000 random samples from each `f_test_vec[i]`, which are independent component draws. The joint sum is used to estimate the combined variance. With O(10) components each requiring 100 000 evaluations of a 60 000-point TF1, this is extremely slow.

---

## 5. master_cov_matrix (CovMatrix) — Covariance Matrix Construction

### Purpose

`master_cov_matrix.cxx` (with included fragments `mcm_1.h`, `mcm_2.h`, `mcm_data_stat.h`, `mcm_pred_stat.h`) implements `LEEana::CovMatrix`. It reads configuration files describing analysis channels, MC/data files, and systematic weights, then constructs the full multi-systematic covariance matrix used downstream by `TLee`. It also fills prediction and data histograms from event-by-event TTree information.

### Key Classes and Methods

| Class / Method | File : Lines | Role |
|---|---|---|
| `CovMatrix(cov_filename, cv_filename, file_filename, rw_filename)` | `master_cov_matrix.cxx:46` | Parses four config text files; initialises channel maps, file lists, reweighting info |
| `gen_xf_cov_matrix(run, ...)` | `mcm_2.h:4` | Computes flux+XS fractional covariance by iterating over reweight universes |
| `gen_det_cov_matrix(run, ..., flag_gp)` | `mcm_1.h:10` | Bootstraps detector-systematic covariance (1000 throws); calls `GPSmoothing` |
| `gen_data_stat_cov_matrix(run, ...)` | `mcm_data_stat.h:1` | Bootstraps data-statistical covariance (5000 throws via `TPrincipal`) |
| `gen_pred_stat_cov_matrix(run, ...)` | `mcm_pred_stat.h:1` | Propagates MC-prediction statistical uncertainty |
| `fill_pred_histograms(run, ...)` | `master_cov_matrix.h:91` | Fills prediction histograms from TTrees, applying POT normalisation and optional LEE weight |
| `fill_data_histograms(run, ...)` | `master_cov_matrix.h:90` | Fills data observation histograms |
| `get_events_weights(input_filename, ...)` | `master_cov_matrix.h:106` | Event-loop over TTree; populates map of passed events with weights |
| `get_events_info(input_filename, ...)` | `master_cov_matrix.h:125` | Event-loop for detector systematics; reads detector-variation weights |

### Inputs Consumed / Outputs Produced

- **Config files read:** `configurations/cov_input.txt`, `configurations/cv_input.txt`, `configurations/file_ch.txt`, `configurations/rw_cv_input.txt`, optionally `configurations/gp_input.txt`, `configurations/xs_ch.txt`, `configurations/osc_parameter.txt`.
- **ROOT files read:** MC and data ntuples listed in `cv_input.txt` / `file_ch.txt`; TTree names are constructed from channel-cut configuration.
- **Outputs written / filled:** `TMatrixD* cov_xf_mat`, `TMatrixD* cov_det_mat`, `TMatrixD* cov_mat_bootstrapping`, and all histogram maps passed by reference; optionally `debug_smoothing.root`.

### Algorithm Sketch

The constructor reads text configuration files to build internal maps relating channels to file types, histogram binning, and systematic categories. `gen_xf_cov_matrix` iterates over each reweight universe (one column of the weight vector per event), accumulates the outer product of per-universe spectrum shifts, and sums to form the fractional covariance. `gen_det_cov_matrix` bootstraps 1000 detector-variation samples using `TPrincipal`; the resulting bootstrapped covariance is passed through `GPSmoothing` before a second `TPrincipal` step that converts the (possibly non-PD) smoothed matrix into a proper covariance. `gen_data_stat_cov_matrix` bootstraps 5000 data-statistical samples similarly. → see 07_algorithms.md#covariance-matrix-construction

### Notes and Red Flags

- **`mcm_2.h:1990`** (approximately line 1990 of `master_cov_matrix.cxx` as included) — In the `"reweight"` branch of the universe-loop, `gRandom->SetSeed(j*reweight*77777)` is called inside the innermost loop over 1000 universe indices `j`. Because `gRandom` is a global `TRandom` (not `TRandom3`), and the seed is `j * reweight * 77777` (a floating-point product cast implicitly to `UInt_t`), the same seed will be reused whenever `reweight` rounds to the same integer multiple. This breaks the independence of universe throws. → see 05_bugs.md#B-08
- **`master_cov_matrix.cxx:41`** — The four implementation fragments (`mcm_1.h`, `mcm_2.h`, `mcm_data_stat.h`, `mcm_pred_stat.h`) are `#include`d directly into the `.cxx` translation unit rather than being compiled as separate units. This means all four "header" files are actually implementation files, creating a single monolithic translation unit that is difficult to build incrementally or test in isolation.
- **`mcm_1.h:57`** — `double data_pot = 5e19;` is hardcoded as a local initialisation before being overwritten from `map_inputfile_info`. If the lookup fails (e.g. due to a missing run period), the fallback silently uses 5×10¹⁹ POT rather than signalling an error.
- **`mcm_1.h:77`** — The bootstrapping loop (`for (int qx = 0; qx != 1000; qx++)`) uses `TPrincipal` to accumulate rows, then retrieves the covariance matrix. `TPrincipal` is designed for PCA, not covariance estimation; while it does store the covariance matrix, this is an unusual and fragile choice compared to direct accumulation.
- **`master_cov_matrix.cxx` constructor, line 75`** — `flag_reweight` is set to `true` if any single reweight entry has `flag_reweight_i==1`, but `rw_type` is never updated from its default of 0. Code that branches on `rw_type` (e.g. the `rw_type==3` block at `mcm_2.h:2026`) will therefore never execute in a straightforward configuration.

---

## 6. TLee — Top-Level Statistical Engine

### Purpose

`TLee.cxx` / `inc/WCPLEEANA/TLee.h` implement the `TLee` class, the outermost statistical engine that ingests pre-built spectra and covariance matrices and performs the full LEE analysis: collapse to observation bins, goodness-of-fit tests with and without sideband constraint, Minuit2 minimisation of the LEE signal strength, and Feldman-Cousins confidence interval generation. It is the last object to run in the pipeline, consuming outputs from `CovMatrix`.

### Key Classes and Methods

| Method | File : Lines | Role |
|---|---|---|
| `Set_config_file_directory(...)` | `TLee.cxx:2410` | Stores paths to the spectra ROOT file and covariance directories |
| `Set_Spectra_MatrixCov()` | `TLee.cxx:2426` | Reads spectra and all covariance matrices from ROOT/text files; populates all internal maps |
| `Set_POT_implement()` | `TLee.cxx:2341` | Scales all spectra and covariances by `scaleF_POT²` |
| `Set_TransformMatrix()` | `TLee.cxx:2331` | Placeholder — body is empty |
| `Set_Collapse()` | `TLee.cxx:2218` | Applies LEE scale factor and collapse transform matrix; builds `matrix_pred_newworld` and `matrix_absolute_cov_newworld` |
| `GetChi2(pred, meas, cov)` | `TLee.cxx:537` | Computes CNP-corrected chi² for a given prediction/measurement pair |
| `Exe_Goodness_of_fit(num_Y, num_X, ...)` | `TLee.cxx:1084` | Full GoF with optional sideband constraint; produces plots and fills `val_GOF_noConstrain` / `val_GOF_wiConstrain` |
| `Exe_Goodness_of_fit(vc_target, vc_support, index)` | `TLee.cxx:984` | Channel-level wrapper for the low-level GoF function |
| `Exe_Goodness_of_fit_detailed(...)` | `TLee.cxx:949` | Bin-level wrapper for GoF |
| `Set_Variations(num_toy)` | `TLee.cxx:499` | Generates `num_toy` pseudo-experiments via eigendecomposition of covariance + Poisson draws |
| `Minimization_Lee_strength_FullCov(Lee_init, flag_fixed)` | `TLee.cxx:219` | Minuit2 MIGRAD minimisation of CNP chi² over Lee signal strength; fills `minimization_*` members |
| `Exe_Feldman_Cousins(low, hgh, step, num_toy, ifile)` | `TLee.cxx:151` | Full FC scan writing `file_FC_NNNNNN.root` |
| `Exe_Fledman_Cousins_Asimov(...)` | `TLee.cxx:98` | Asimov FC scan writing `file_Asimov.root` |
| `Exe_Fiedman_Cousins_Data(...)` | `TLee.cxx:44` | Data FC scan writing `file_data.root` |

### Inputs Consumed / Outputs Produced

- **Inputs read by `Set_Spectra_MatrixCov`:**
  - `spectra_file` ROOT file: histograms `histo_N` (prediction channels), `hdata_obsch_N` (data), `mat_collapse` (TMatrixD collapse/transform), `cov_mat_add` (additional covariance).
  - `flux_Xs_directory/cov_N.root`: fractional covariance matrices `frac_cov_xf_mat_N`.
  - `detector_directory/cov_*.root`: fractional detector covariance matrices `frac_cov_det_mat_N`.
  - `mc_directory/N.log`: text files with per-bin MC statistical uncertainties as a function of LEE strength.
- **Outputs written:**
  - `file_FC_NNNNNN.root` (TTree `tree`): Feldman-Cousins toys.
  - `file_Asimov.root` (TTree `tree_Asimov`): Asimov FC.
  - `file_data.root` (TTree `tree_data`): data FC scan.
  - PNG figures from `Plotting_singlecase` when `saveFIG=true`.

### Algorithm Sketch

After loading (`Set_Spectra_MatrixCov`) and POT-scaling (`Set_POT_implement`), `Set_Collapse` applies the sparse collapse matrix to project the "old-world" binning (all MC channels combined) into the "new-world" observation binning, simultaneously scaling the LEE-signal channels by `scaleF_Lee`. The absolute covariance is formed as Cₙₑₓ = Tᵀ·Cₒₗₐ·T where T is the collapse transform. `GetChi2` uses the CNP (Combined Neyman-Pearson) statistic: for bins with prediction below a DocDB-32520 threshold table, the diagonal statistical term is replaced by `(pred−meas)²/[2(pred−meas+meas·log(meas/pred))]`; otherwise `pred` is used as the Poisson variance. `Exe_Goodness_of_fit` computes both a free chi² (no constraint) and a constrained chi² where support channels condition the target prediction via the Schur complement. `Minimization_Lee_strength_FullCov` wraps the chi² in a Minuit2 lambda functor with LEE strength as the single free parameter. Feldman-Cousins intervals are obtained by scanning `scaleF_Lee` on a grid and generating MC toys at each true LEE value. → see 07_algorithms.md#cnp-chi2, 07_algorithms.md#sideband-constraint, 07_algorithms.md#feldman-cousins

### Notes and Red Flags

- **`TLee.cxx:2505`** — `TFile *file_spectra = new TFile(roostr, "read")` is opened inside `Set_Spectra_MatrixCov` but is never closed. All `TFile*` objects opened for flux-XS covariance matrices in the subsequent loop (`map_file_flux_Xs_frac[idx] = new TFile(...)`, lines 2613–2614) and detector files (`map_file_detector_frac[idx] = new TFile(...)`, lines 2666) are similarly never closed. This constitutes a TFile handle leak on every call to `Set_Spectra_MatrixCov`. → see 05_bugs.md#B-09
- **`TLee.cxx:1130–1230`** — There are seven consecutive `if (index==0)` blocks. Each block unconditionally overwrites the `flag_axis_userAA`, `flag_axis_userAB`, and `userAA_*`/`userAB_*` variables set by all previous blocks. Only the last such block (lines 1205–1216) has any effect at runtime. The preceding six blocks are dead code that misleads readers about which axis parameters are active. → see 05_bugs.md#B-10
- **`TLee.cxx:253–261`** — Inside the CNP chi² computation in `Minimization_Lee_strength_FullCov`, the statistical covariance diagonal is computed as `3/(1/val_meas + 2/val_pred)` when `val_meas != 0`. When `val_meas == 0`, the code sets `val_stat_cov = val_pred/2`. However, the denominator `1/val_meas` is evaluated only after the `val_meas==0` branch check; the `val_meas==0` case correctly falls through to `val_pred/2`. The CNP formula for `val_meas==0` is `val_pred/2` (correct per the CNP prescription), but this differs from the formula used in `GetChi2` at line 560 (`matrix_stat_cov(idx,idx) = matrix_pred_temp(0,idx)` as the default). The two chi² paths are inconsistent in their treatment of the `meas=0` case. → see 05_bugs.md#B-11
- **`TLee.cxx:2331–2337`** — `Set_TransformMatrix()` has an empty body with only a print statement. It is referenced in comments as corresponding to `Set_Spectra_MatrixCov`, suggesting the transform-matrix setup was intended to be refactored out but was never completed.
- **`TLee.cxx:44–95` and `TLee.cxx:98–149`** — The Feldman-Cousins method names contain typographic inconsistencies: `Exe_Fiedman_Cousins_Data` (line 44, "Fiedman") and `Exe_Fledman_Cousins_Asimov` (line 98, "Fledman"). The correct method `Exe_Feldman_Cousins` at line 151 uses the proper spelling. These typos propagate to the public header (`TLee.h:202–203`).

---

*End of document. Cross-references to 05_bugs.md use placeholder IDs B-01 through B-11 assigned above; the bugs document should adopt these IDs. Cross-references to 07_algorithms.md use anchor names that the algorithms document should provide.*
