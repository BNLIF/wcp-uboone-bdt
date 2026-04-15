# Algorithm Reference: MicroBooNE LEE Analysis

This document is a technical reference for the statistical and numerical algorithms
used in the wcp-uboone-bdt analysis framework. Every claim is tied to a specific
file and line number so that a new group member can read the code alongside this text.

---

## Table of Contents

1. [CNP chi-squared with prediction protection](#1-cnp-chi-squared-with-prediction-protection)
2. [Conditional-constraint goodness-of-fit (Schwartz decomposition)](#2-conditional-constraint-goodness-of-fit-schwartz-decomposition)
3. [Feldman-Cousins delta-chi-squared sensitivity scans](#3-feldman-cousins-delta-chi-squared-sensitivity-scans)
4. [LEE signal-strength fit](#4-lee-signal-strength-fit)
5. [Systematic covariance pipeline](#5-systematic-covariance-pipeline)
6. [Wiener-SVD unfolding](#6-wiener-svd-unfolding)
7. [Bayesian MC-stat convolution (Bayes class)](#7-bayesian-mc-stat-convolution-bayes-class)
8. [Gaussian-process regression (GP kernel and regressor)](#8-gaussian-process-regression-gp-kernel-and-regressor)

---

## 1. CNP chi-squared with prediction protection

### Formula

The Combined-Neyman-Pearson (CNP) statistic is designed to handle low-count bins
where the standard Pearson chi-squared breaks down. The diagonal statistical
covariance for bin i is:

```
sigma_stat^2(i) = 3 / (1/obs + 2/pred)      [obs > 0]
sigma_stat^2(i) = pred/2                      [obs = 0]
sigma_stat^2(i) = 1e-6                        [obs = 0 and pred = 0]
```

This is the harmonic mean of obs and pred with weights 1 and 2 respectively.
When obs >> pred it approaches pred; when pred >> obs it approaches obs/2.

The full chi-squared is then the matrix form:

```
chi2 = delta^T * (Cov_stat + Cov_syst)^{-1} * delta
```

where `delta = pred - obs` is a row vector.

### Implementation mapping

The formula is computed in `src/TLee.cxx:537-571` inside
`TLee::GetChi2(matrix_pred_temp, matrix_meas_temp, matrix_syst_abscov_temp)`.

| Formula term | Variable | Location |
|---|---|---|
| delta | `matrix_delta` | TLee.cxx:541 |
| sigma_stat^2(i) diagonal | `matrix_stat_cov(idx,idx)` | TLee.cxx:551-563 |
| Cov_syst | `matrix_syst_abscov_temp` | argument passed in |
| Cov_total = Cov_stat + Cov_syst | `matrix_total_cov` | TLee.cxx:567 |
| (Cov_total)^{-1} | `matrix_total_cov_inv` | TLee.cxx:568 |
| chi2 | `chi2 = (delta * inv * delta^T)(0,0)` | TLee.cxx:569 |

The same diagonal formula is also computed inline at TLee.cxx:254-262 (inside the
Minuit2 FCN lambda) and at TLee.cxx:1290-1320 (goodness-of-fit path).

### Prediction protection table

When integer observation counts are small (1-10 events) and the prediction is
very low, the CNP denominator can become numerically unstable.
The code replaces `sigma_stat^2` with the Neyman-Poisson estimator
`(pred-obs)^2 / [2(pred - obs + obs*ln(obs/pred))]` only when
`pred < array_pred_protect[obs]`. The threshold table, referenced to docDB-32520:

```c++
double array_pred_protect[11] = {0, 0.461, 0.916, 1.382, 1.833, 2.298,
                                  2.767, 3.225, 3.669, 4.141, 4.599};
```
(TLee.cxx:547)

So for obs=1 the protection kicks in only if pred < 0.461; for obs=5 only if pred < 2.298.
Index 0 is unused (the condition requires obs >= 1).

**Where in the code:**
- Primary definition: `src/TLee.cxx:537-572` (`TLee::GetChi2`)
- Inline copies (FCN lambda and GoF path): `src/TLee.cxx:254-263`, `1289-1319`

---

## 2. Conditional-constraint goodness-of-fit (Schwartz decomposition)

### Purpose

The analysis observes many sidebands (control regions) in addition to the
signal-region bins. After fitting all bins simultaneously, a goodness-of-fit
test on the signal region alone would be misleading because the prediction there
has been updated ("pulled") by the observed control-region counts.
The Schwartz (multivariate-normal conditional) update gives the correct
posterior mean and covariance for the Y bins given what was observed in the X bins.

### The constraint formulas

Let the full prediction vector and full covariance matrix be partitioned as:

```
mu = [ mu_Y ]     Sigma = [ Sigma_YY   Sigma_YX ]
     [ mu_X ]             [ Sigma_XY   Sigma_XX ]
```

Given observations x in the X (sideband) bins, the updated Y prediction is:

```
mu_Y' = mu_Y + Sigma_YX * Sigma_XX^{-1} * (x - mu_X)
```

and the updated covariance is:

```
Sigma_YY' = Sigma_YY - Sigma_YX * Sigma_XX^{-1} * Sigma_XY
```

These are the standard multivariate-normal conditional update equations.

### Bin assignment in TLee

The primary goodness-of-fit call assigns bins as follows (TLee.cxx:327-333):

- **Y (signal region, target):** channels 0-7 and 26-33
  (`vc_target_detailed_chs` contains `idx` for `idx in [0,8)` and `[26,34)`)
- **X (sidebands, support):** channels 8-25, 34-51, and 52-136
  (`vc_support_detailed_chs` contains `idx` for `idx in [8,26)`, `[34,52)`, `[52,137)`)

The index 26 boundary separates the two neutrino-energy reconstruction categories
(FC and PC). Within each category the first 8 bins (low reconstructed energy)
form the signal window.

### Matrix construction

Inside `TLee::Exe_Goodness_of_fit(int num_Y, int num_X, ...)` at TLee.cxx:1084:

| Step | Code | Line |
|---|---|---|
| Extract Sigma_XX | `matrix_XX = matrix_cov_total.GetSub(num_Y, num_Y+num_X-1, num_Y, num_Y+num_X-1)` | 1522 |
| Invert Sigma_XX | `matrix_XX_inv.Invert()` | (after GetSub) |
| Extract Sigma_YX | `matrix_YX = matrix_cov_total.GetSub(0, num_Y-1, num_Y, num_Y+num_X-1)` | 1522-area |
| Compute mu_Y' | `matrix_Y_under_X = matrix_pred_Y + matrix_YX*matrix_XX_inv*(matrix_meas_X-matrix_pred_X)` | TLee.cxx:386 |
| Compute Sigma_YY' | `matrix_YY_under_XX = matrix_YY - matrix_YX*matrix_XX_inv*matrix_XY` | TLee.cxx:387 |
| Chi-squared on Y | `matrix_wicons_chi2 = delta * Sigma_YY'^{-1} * delta^T` | TLee.cxx:401 |

The no-constraint chi-squared (using only Sigma_YY, not the conditioned version) is
computed just before at TLee.cxx:1322-1330 for comparison.

**Where in the code:**
- Driver: `src/TLee.cxx:1084` (`TLee::Exe_Goodness_of_fit`)
- Constraint matrices inside FCN with the active code path: `src/TLee.cxx:377-402`
- Bin assignment: `src/TLee.cxx:327-333`

---

## 3. Feldman-Cousins delta-chi-squared sensitivity scans

### The test statistic

For each hypothesized signal strength mu_true, the FC statistic is:

```
delta_chi2(mu_true) = chi2(mu_true, fixed) - chi2(mu_best_fit, free)
```

where `chi2(mu, fixed)` evaluates the chi-squared with Lee_strength fixed to mu,
and `chi2(mu_best_fit, free)` minimizes over Lee_strength. The confidence interval
for a given data set is the set of mu_true values for which delta_chi2 is below
the critical value determined from the toy distribution at that mu_true.

### The three variants

**`Exe_Fiedman_Cousins_Data`** (TLee.cxx:44, note the typo in the method name):
Operates on real (or provided fake) data. Computes the global best-fit Lee_strength
once (`Minimization_Lee_strength_FullCov(1, 0)`), then scans mu_true values and
for each evaluates `chi2(mu_true, fixed)`. Output: a TTree named `tree_data` in
`file_data.root` with branches `Lee_bestFit_data`, `chi2_gmin_data`, and the scan
vector `chi2_null_scan_data`.

**`Exe_Fledman_Cousins_Asimov`** (TLee.cxx:98, note the typo):
For each mu_true value, sets the Asimov data set (`Set_toy_Asimov`, which sets
fake data equal to the prediction with no Poisson fluctuation) and evaluates the
chi-squared at each scanned hypothesis. There are no random toys; this gives the
median expected sensitivity.

**`Exe_Feldman_Cousins`** (TLee.cxx:151):
Full toy-based FC scan. For each mu_true:
1. `Set_Collapse()` builds the prediction at Lee_strength = mu_true (TLee.cxx:179).
2. `Set_Variations(num_toy)` generates `num_toy` fluctuated spectra (TLee.cxx:181).
3. The toy loop (TLee.cxx:189-203) calls `Set_toy_Variation(itoy)` then runs
   `Minimization_Lee_strength_FullCov` twice: once with fixed mu_true to get
   `chi2_null`, once free to get `chi2_gmin`. The difference is stored per toy.

### Toy generation: eigendecomposition + Poisson

`Set_Variations` (TLee.cxx:499-532):
1. Constructs the symmetric covariance matrix `DSmatrix_cov` from
   `matrix_absolute_cov_newworld` (TLee.cxx:505-510).
2. Eigendecomposes it: `TMatrixDSymEigen DSmatrix_eigen(DSmatrix_cov)` (TLee.cxx:511).
3. For each toy, draws independent Gaussian random numbers scaled by
   `sqrt(eigenvalue[i])` (TLee.cxx:519), transforms back to the original bin
   basis via `matrix_eigenvector * matrix_element` (TLee.cxx:525).
4. Adds the systematic variation to the prediction, clamps to zero, then
   Poisson-fluctuates: `rand->PoissonD(val_with_syst)` (TLee.cxx:529).

**Where in the code:**
- `Exe_Fiedman_Cousins_Data`: `src/TLee.cxx:44-95`
- `Exe_Fledman_Cousins_Asimov`: `src/TLee.cxx:98-149`
- `Exe_Feldman_Cousins` (full toys): `src/TLee.cxx:151-214`
- Toy loop: `src/TLee.cxx:189-203`
- Eigendecomposition toy throw: `src/TLee.cxx:499-532`

---

## 4. LEE signal-strength fit

### Parameterization

The prediction is:

```
pred_i(Lee_strength) = sum_j  T_ij(Lee_strength) * input_j
```

where `T_ij` is the transformation matrix that maps old-world bins to new-world
(observable) bins. For rows corresponding to LEE signal bins, each entry of
`matrix_transform_Lee` is multiplied by `scaleF_Lee` (TLee.cxx:2225):

```c++
if( map_Lee_oldworld.find(ibin)!=map_Lee_oldworld.end() )
    matrix_transform_Lee(ibin, jbin) *= scaleF_Lee;
```

So varying `Lee_strength` from 0 to 1 interpolates between pure background and
the full LEE signal template; values greater than 1 over-predict the LEE excess.

### What Set_Collapse does

`TLee::Set_Collapse()` (TLee.cxx:2218) is called every time `scaleF_Lee` changes:
1. Builds `matrix_transform_Lee` by scaling the Lee-signal rows of `matrix_transform`
   by `scaleF_Lee`.
2. Projects the old-world prediction: `matrix_pred_newworld = matrix_pred_oldworld * matrix_transform_Lee`.
3. Projects all systematic covariance matrices similarly:
   `matrix_absolute_cov_newworld = T^T * matrix_absolute_cov_oldworld * T`.
4. Adds the MC-stat contribution bin-by-bin from a precomputed TGraph:
   `gh_mc_stat_bin[ibin]->Eval(scaleF_Lee)` (TLee.cxx:2257). This TGraph was
   built at analysis startup by scanning over Lee_strength values and computing
   the Bayesian MC-stat variance at each one (TLee.cxx:2829-2840).

### Minuit2 minimization

`TLee::Minimization_Lee_strength_FullCov` (TLee.cxx:219) sets up a Minuit2 MIGRAD
minimization over a single parameter `Lee_strength`.

- The FCN lambda is defined at TLee.cxx:232 and closed at line 413.
- Inside the FCN: `scaleF_Lee = Lee_strength`, `Set_Collapse()`, then the CNP
  chi-squared is computed (TLee.cxx:246-276).
- The lower bound of Lee_strength is 0 (TLee.cxx:421).
- If `flag_fixed=true` the parameter is frozen (TLee.cxx:423): used when evaluating
  chi-squared at a specific hypothesis point rather than minimizing.
- `min_Lee.Minimize()` calls MIGRAD (TLee.cxx:427).
- Results stored in `minimization_Lee_strength_val` and `minimization_chi2`
  (TLee.cxx:438-439).

**Where in the code:**
- `Set_Collapse`: `src/TLee.cxx:2218-2260`
- FCN lambda: `src/TLee.cxx:232-413`
- MIGRAD call: `src/TLee.cxx:427`
- MC-stat TGraph eval: `src/TLee.cxx:2257`

---

## 5. Systematic covariance pipeline

The full systematic covariance is assembled from four independent sources.
Each is computed in a separate include file, all pulled into
`src/master_cov_matrix.cxx` (lines 41-44).

### 5a. Flux and cross-section covariance (`gen_xf_cov_matrix`, mcm_2.h)

Located in `src/mcm_2.h:4`, this function iterates over neutrino interaction
universes (one per flux or cross-section systematic parameter variation).
For each universe `i` within a given systematic group `j`, the vector of
bin counts is collected into `x[]` (mcm_2.h:141-144) and an outer product
`temp_mat(n,m) += x[n] * x[m]` accumulates the second-moment matrix (mcm_2.h:147-151).
After dividing by `nsize` (the number of universes), this gives the unnormalized
covariance for that group. All groups are summed into `cov_xf_mat` (mcm_2.h:163).
The resulting matrix is an absolute (not fractional) covariance in units of counts^2.

The outer loop over systematic groups is at `src/mcm_2.h:73`.

### 5b. Detector systematic covariance (`gen_det_cov_matrix`, mcm_1.h)

Located in `src/mcm_1.h:10`. This is a two-stage procedure:

**Stage 1 — 1000-sample bootstrap via TPrincipal (mcm_1.h:77-143):**
The detector variation samples are drawn by re-sampling events according to their
per-event weights (Poisson bootstrap). For each of 1000 iterations the full
POT-normalized histogram is filled, and each resulting count vector is
added to a `TPrincipal` object (mcm_1.h:139). After the loop,
`*cov_mat_bootstrapping = *(TMatrixD*)prin.GetCovarianceMatrix()` (mcm_1.h:143),
and `*vec_mean_diff = *(prin.GetMeanValues())` captures the mean fractional shift
of the detector variation with respect to the central value (mcm_1.h:149).
This mean shift vector is then passed to GP smoothing (mcm_1.h:152).

**Stage 2 — 16000 Gaussian throws with global amplitude modulation (mcm_1.h:169-188):**
The bootstrapped covariance is eigendecomposed (mcm_1.h:161-165).
For each of 16000 throws:
- A correlated shift is drawn: `matrix_element(j,0) = Gaus(0, sqrt(eigenvalue(j)))` (mcm_1.h:173).
- The shift is transformed back to bin space: `matrix_variation = eigenvector * matrix_element` (mcm_1.h:177).
- The GP-smoothed detector bias `vec_mean_diff` is added: `matrix_variation(j,0) += (*vec_mean_diff)(j)` (mcm_1.h:181).
- A single global relative-error multiplier is drawn: `rel_err = random3.Gaus(0,1)` (mcm_1.h:178).
- The stored vector is `x[j] = rel_err * matrix_variation(j,0)` (mcm_1.h:183).
- All 16000 rows are added to a second `TPrincipal`; its covariance becomes the
  final detector covariance `cov_det_mat` (mcm_1.h:189-194).

The purpose of the `rel_err` amplification is to allow the overall detector
systematic scale to float as if it were an unconstrained nuisance parameter,
producing the correct marginal uncertainty on the bin counts.

### 5c. MC-stat covariance (`gen_pred_stat_cov_matrix`, mcm_pred_stat.h)

Located in `src/mcm_pred_stat.h:1`. The MC statistical uncertainty comes from the
finite number of simulated events. For each of 5000 bootstrap iterations
(`nround = 5000`, mcm_pred_stat.h:136), each MC event is independently reweighted
by drawing a new per-event weight from a Bayesian posterior (via the `Bayes` class,
see Section 7). The varied histogram gives an alternate count vector `x[]`.
The covariance is accumulated as:
`cov_mat(n,m) += (x[n] - mean[n]) * (x[m] - mean[m]) / nround` (mcm_pred_stat.h:185).

### 5d. Data-stat covariance (`gen_data_stat_cov_matrix`, mcm_data_stat.h)

Located in `src/mcm_data_stat.h:1`. Uses a 5000-iteration Poisson bootstrap
(mcm_data_stat.h:69). In each iteration, each data event is included 0 or more
times according to a Poisson(1) draw. The resulting count histogram is submitted
to a `TPrincipal` accumulator. The final covariance is extracted from the
principal component analysis.

### 5e. GP smoothing of detector bias (GPSmoothing.C, mcm_1.h)

The mean-shift vector `vec_mean_diff` (the average detector-variation minus CV
fractional shift) is noisy because only 1000 bootstrap samples are available.
Before using this vector to bias the 16000-throw Stage 2 samples, it is smoothed
by a GP posterior mean.

`GPSmoothing(vec_mean_diff, cov_mat_bootstrapping, ...)` is called at mcm_1.h:152.
Inside `src/GPSmoothing.C:50`, an `RBFKernel` is constructed from parameters read
from a configuration file (`./configurations/gp_input.txt`), then a `GPRegressor`
is fitted to the noisy bias vector (GPSmoothing.C:113-116). The GP posterior mean
`reg.PosteriorMean()` replaces `vec_mean_temp`, and the GP posterior covariance
`reg.PosteriorCov()` replaces `cov_mat_bootstrapping_temp` (GPSmoothing.C:118-125).
These smoothed values are what enter Stage 2.

**Where in the code:**
- Flux/XS outer loop: `src/mcm_2.h:73`
- Detector Stage 1 bootstrap: `src/mcm_1.h:77-143`
- Detector Stage 2 Gaussian throws: `src/mcm_1.h:169-188`
- MC-stat bootstrap: `src/mcm_pred_stat.h:136-188`
- Data-stat bootstrap: `src/mcm_data_stat.h:69`
- GP smoothing call: `src/mcm_1.h:152`

---

## 6. Wiener-SVD unfolding

### The unfolding model

The measurement vector M is related to the true signal S by:

```
M = R * S + background
```

where R is the response matrix (maps signal space bins to measurement space bins).
The goal is to infer S from M after subtracting background.

### Step 1: whitening

The measurement covariance Cov is decomposed by SVD:
`TDecompSVD decV(Covariance)` (WienerSVD.cxx:125).
Let `Cov = V * diag(sigma^2) * V^T`. The whitening matrix is:
`Q = diag(1/sigma) * V^T` (WienerSVD.cxx:126-143).

Applying Q rotates and scales the measurement space so that all bins have
unit variance: `M' = Q*M`, `R' = Q*R` (WienerSVD.cxx:145-146).

### Step 2: regularization matrix C

`Matrix_C(n, type)` (WienerSVD.cxx:29) builds an n×n regularization matrix:

| Type | Meaning | Matrix form |
|---|---|---|
| 0 | identity (no regularization) | C = I |
| 1 | first-derivative | C(i,i)=-1, C(i,i+1)=+1 |
| 2 | second-derivative | C(i,i)=-2+eps, C(i,j)=+1 for |i-j|=1 |
| 23 | second-derivative, stitched 3D | type-2 but with zero off-diagonals at slice boundaries |
| 33 | third-derivative, stitched 3D | higher-order; breaks at slice boundaries |
| 332/333 | full 3D via C3_3D helper | calls WienerSVD_3D.C |

For types 22/32 (2D stitched) the slice boundary indices are hardcoded as:
`dim_edges = { 0, 3, 7, 11, 14, 18, 22, 26, 31, 36 }` (WienerSVD.cxx:40).
These represent 9 angular slices of a muon-neutrino cross-section measurement
with varying numbers of momentum bins per slice.

For types 23/33 (3D stitched) the active boundary table is:
```
dim_edges = { 0, 3, 6, 10, 13, 16, 19, 22, 25,
             28, 31, 35, 39, 42, 45, 50, 55, 60,
             63, 66, 70, 74, 77, 82, 88, 94, 100,
            105,108,112,115,118,122,126,129,133, 138 }
```
(WienerSVD.cxx:47-50). This defines 35 slices across 138 total bins (4 Enu × 9
costheta × variable Pmu bins), matching `nbins_tot_3D = 138` in WienerSVD_3D.C:20.
At each edge index the second-derivative stencil is truncated to prevent coupling
across physically distinct kinematic slices.

The signal vector is also pre-conditioned by a normalization matrix:
`C0 = C0 * normsig` where `normsig(i,i) = 1/Signal(i)^Norm_type` (WienerSVD.cxx:157-160).

### Step 3: SVD of the conditioned response

`TDecompSVD udv(R')` decomposes the whitened, C-conditioned response matrix into
`R' = U * D * V^T` (WienerSVD.cxx:169-174).

### Step 4: Wiener filter

In the rotated V^T basis, the signal power spectrum is `S_i = (V^T * C * Signal)(i)`.
The Wiener filter diagonal matrix W is (WienerSVD.cxx:211):

```
W(i,i) = flag_WienerFilter * S(i)^2 / (D(i)^2 * S(i)^2 + 1)
```

where:
- `S(i)` is the i-th component of the signal in the V basis (signal power)
- `D(i)` is the i-th singular value of R' (response strength)
- `flag_WienerFilter` is a normalization constant (nominally 1)
- The denominator `D^2 * S^2 + 1` balances signal-to-noise against regularization

When `flag_WienerFilter = 0` the filter degenerates to `W(i,i) = 1/D(i)^2`
(pure pseudo-inverse, no regularization) (WienerSVD.cxx:217).

### Unfolded result and covariance

The unfolded signal is (WienerSVD.cxx:225):
```
unfold = C^{-1} * V * W * D^T * U^T * M'
```

The additional smearing matrix that maps the true signal to the unfolded space is
(WienerSVD.cxx:226):
```
AddSmear = C^{-1} * V * W0 * V^T * C
```
where `W0(i,i) = D(i)^2 * W(i,i)`.

The unfolded covariance is:
`UnfoldCov = covRotation * Covariance * covRotation^T` (WienerSVD.cxx:232)
where `covRotation = C^{-1} * V * W * D^T * U^T * Q`.

**Where in the code:**
- Main function: `src/WienerSVD.cxx:119-236`
- Regularization matrix: `src/WienerSVD.cxx:29-110`
- Wiener filter: `src/WienerSVD.cxx:194-222`
- 3D bin geometry: `src/WienerSVD_3D.C:17-20`

---

## 7. Bayesian MC-stat convolution (Bayes class)

### The propagated-Poisson PDF

The core likelihood function is `LEEana::Prop_Poisson_Pdf` (`src/bayes.cxx:377`).
For a component with measured count `meas`, variance `sigma2`, and POT scaling
weight `weight`, the function computes the likelihood of true expectation `mu` as:

```
eff_weight = sigma2 / meas       (effective weight, i.e. mean MC weight)
eff_meas   = meas / eff_weight   (effective number of unweighted events)

P(mu | meas, sigma2, weight) =
  exp(eff_meas - mu/(eff_weight*weight)
      + eff_meas * log(mu / (eff_weight * eff_meas * weight)))
```

This is proportional to a Poisson likelihood with effective count `eff_meas`
evaluated at a rescaled mean `mu/(eff_weight*weight)`. It reduces to a standard
Poisson when all MC events have equal weight (sigma2 = meas, so eff_weight = 1).
When meas = 0 it returns `exp(-mu/(eff_weight*weight))` (bayes.cxx:408).

### FFT convolution for multi-component bins

Each reconstructed bin receives contributions from multiple MC samples
(track-like, shower-like, dirt, etc.), each with its own effective weight.
The total count distribution is the convolution of per-component distributions.

In `do_convolution()` (bayes.cxx:104):
1. The first component's PDF is stored as `f_conv` (bayes.cxx:138).
2. Each subsequent component is convolved with the running product using
   `TF1Convolution(f_conv, f1, llimit, hlimit, true)` (bayes.cxx:143).
   The `true` flag requests FFT-based convolution.
3. `conv->SetNofPointsFFT(10000)` sets the FFT grid (bayes.cxx:144).
4. Each intermediate result is stored as a new TF1 with NPX=60000 (bayes.cxx:134).
5. After all convolutions, the result is sampled onto a TGraph `g1` with 5000 points
   for fast evaluation (bayes.cxx:162-166), then wrapped in a TF1 `f_conv_num`
   with NPX=5000 (bayes.cxx:172-177).

### The reweight correction in `get_covariance`

`get_covariance()` (bayes.cxx:294) integrates `sum += P(x) * (x - mean)^2 * corr`
over the total PDF, where the correction factor is:

```
corr = exp(-(num_component - 1) * log(x / mean))
     = (x / mean)^{-(num_component - 1)}
```

When `mean = 0` the log is taken of x directly (bayes.cxx:312).
This factor downweights large-x tails. Mathematically it undoes the convolution
broadening: summing `num_component` independent variables, each with mean ~mean,
produces an overall distribution with mean = num_component * mean_per_component.
Dividing by this factor corrects for the fact that we want the uncertainty on the
total, not on each component separately.

### Bisection method for percentile extraction

`calculate_lower_upper(nSigma)` (bayes.cxx:183) finds the confidence interval
using bisection. Precision is `analytic_result_precision = 1e-4` (bayes.cxx:206).
The maximum number of iterations is 30 (bayes.cxx:229, 251, 278). The loop
condition is `fabs(val_typeA - val_typeB) > 1e-4` (bayes.cxx:219), which for a
support range of O(100) units requires about log2(100/1e-4) ≈ 20 iterations,
well within the cap.

**Where in the code:**
- `Prop_Poisson_Pdf`: `src/bayes.cxx:377-415`
- `do_convolution` with FFT setup: `src/bayes.cxx:104-179`
- TF1Convolution + SetNofPointsFFT: `src/bayes.cxx:143-144`
- `get_covariance` with reweight: `src/bayes.cxx:294-331`
- Bisection loop: `src/bayes.cxx:219-233` and `241-256`

---

## 8. Gaussian-process regression (GP kernel and regressor)

### RBF kernel: squared Mahalanobis distance

`RBFKernel::Element(pt1, pt2)` (GPKernel.cxx:66) returns:

```
K(pt1, pt2) = coeff * exp(-mag / 2)
```

where `mag` is the squared Mahalanobis distance with 5 separate length scales:

```
mag = sum_{i=0}^{4}  ((x1[i] - x2[i]) / L[i])^2
```

and `coeff = fPars[5]` is the overall amplitude (GPKernel.cxx:69).
The 5 dimensions correspond to the 5 kinematic axes (e.g., Enu, costheta, Pmu,
plus two additional axes from the `GPPoint` 5-element coordinate array).

### Log-scale trick and mutation side-effect

In `RBFKernel::Mag` (GPKernel.cxx:43), if `doLogScales[i]` is true for dimension i:

```c++
p1x[i] = log(p1x[i]);
p2x[i] = log(p2x[i]);
parameters.push_back(log(fPars[i]));
```

(GPKernel.cxx:51-53). This makes the distance metric fractional: a change from
x=1 to x=2 has the same kernel distance as x=2 to x=4. The length scale is also
log-transformed, so `fPars[i]` stores the ratio, not the absolute scale.

**Important bug:** `p1.X()` returns a raw pointer to the internal data of the
`GPPoint` object (which is an array). Writing `p1x[i] = log(p1x[i])` modifies
the point coordinates in place. Since the same `GPPoint` objects are used for
both training and prediction, subsequent calls to `Mag` with the same points will
see already-logarithm-transformed coordinates. This means GP smoothing is only
correct when the log-scale transformation is applied consistently and not mixed
with linear-scale passes. In practice, `GPSmoothing.C` calls `Fit` then `Predict`
with the same point set, so within a single call the mutation is self-consistent.

### Hyperparameter optimization

`GPRegressor::SolveHyperParameters` (GPRegressor.cxx:65) uses ROOT's GSL-BFGS2
minimizer on the negative log-marginal-likelihood:

```
LML = 0.5 * y^T * K^{-1} * y  +  0.5 * log|K|  +  const
```

(MarginalLikelihood::DoEval, GPRegressor.cxx:178-188). The log-determinant uses
the Cholesky upper triangle: `fKDet += 2 * log(choU[i][i])` (GPRegressor.cxx:171).

Hyperparameters are log-transformed before optimization:
`par[i] = log(par_start[i])` (GPRegressor.cxx:74). The optimizer explores the
space of `theta = log(hyperparameter)` on the interval `[-5, 5]` (GPRegressor.cxx:79).
`RBFKernel::SetThetas` (GPKernel.cxx:96-101) undoes the transform: `fPars[i] = exp(thetas[i])`.
If any `par_start[i] <= 0`, `log(par_start[i])` is undefined, which would silently
produce NaN starting values and cause optimizer failure. In `GPSmoothing.C` the
initial parameters are read from the configuration file (GPSmoothing.C:84), so
this must never include zeros or negatives.

In the MicroBooNE analysis, `solveHyperParams = false` is passed to `reg.Fit`
(GPSmoothing.C:115), meaning the hyperparameters are used as given from the
configuration file without further optimization. The `SolveHyperParameters` path
exists but is not exercised in the detector-systematic smoothing workflow.

### Posterior prediction and variance clamping

After fitting, `GPRegressor::Predict` (GPRegressor.cxx:94) computes:

1. Cross-covariance `k(X*, X_train)` = `PosteriorK_T` (GPRegressor.cxx:96).
2. Posterior mean: `fPosteriorMean = k(X*, X_train) * alpha` scaled back from
   normalized space by `fY_Ts` and shifted by `fY_Tm` (GPRegressor.cxx:102-104).
   Here `alpha = K_train^{-1} * y_train` was computed during `Fit` (GPRegressor.cxx:58).
3. Posterior covariance: `PosteriorCov = sigma_11 - k(X*, X_train) * K_train^{-1} * k(X_train, X*)` (GPRegressor.cxx:109-118).
4. Any negative diagonal values are clamped to zero (GPRegressor.cxx:120):
   ```c++
   if(mPosteriorCov(i,i) < 0.) { mPosteriorCov(i,i) = 0.; }
   ```
   Negative posterior variance can arise from floating-point cancellation when the
   posterior uncertainty is very small relative to the prior. Clamping prevents
   downstream NaN propagation.

### Where GP is used in the analysis

- Called from: `src/mcm_1.h:152` as
  `GPSmoothing(vec_mean_diff, cov_mat_bootstrapping, "./configurations/gp_input.txt", flag_gp)`
- Implementation: `src/GPSmoothing.C:50-146`
- Kernel: `src/GPKernel.cxx:43-89`
- Regressor and log-marginal-likelihood: `src/GPRegressor.cxx:9-201`
