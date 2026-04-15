# 08 — Open Questions

This document collects items that require input from Xin (the framework author)
to resolve. Each question arises from a point where the code is internally
consistent but its intent cannot be determined by reading alone — the answer
may classify the item as an intentional physics design choice, a known
approximation, or a defect.

Questions are grouped by component. Each entry gives the file and line range
where the relevant code lives, a one-sentence context, the specific question,
and what the answer changes for bug triage or documentation.

---

### OQ-01  [detector covariance]  src/mcm_1.h:169-188

**Context:** `gen_det_cov_matrix` uses a two-stage Monte Carlo. Stage 1 runs
1000 Poisson-bootstrap throws to estimate the mean detector-variation shift
(`vec_mean_diff`) and the bootstrapped covariance. Stage 2 draws 16000
Gaussian throws from the Stage-1 eigenvectors; each throw `x[j]` is formed as
`x[j] = rel_err * matrix_variation(j,0)`, where `rel_err = random3.Gaus(0,1)`
is a single scalar drawn fresh per throw and applied globally to all bins.

**Question:** Is the global amplitude modulation by `rel_err` intentional? The
current construction folds a unit-normal scalar into every correlated throw,
which inflates the marginal variance of each bin by an additional factor of
`Var(rel_err) = 1` on top of the correlated shape variation. What is the
physical motivation for this term — is it meant to represent an overall
normalization uncertainty for the detector variation, or is it a residual from
an earlier formulation that should have been removed?

**Why it matters:** If intentional, the algorithm documentation (07_algorithms.md
§5b) should explain the physical meaning of the global-scale nuisance. If it is
a leftover term, removing it halves the per-bin marginal variance of the detector
covariance and constitutes a physics-level correction.

---

### OQ-02  [Bayesian MC-stat covariance]  src/bayes.cxx:294-320

**Context:** `Bayes::get_covariance` integrates the convolved posterior PDF
weighted by `corr = (x/mean)^{-(num_component-1)}` (equivalently
`exp(-(num_component-1)*log(x/mean))`). The document 07_algorithms.md §7
describes this as "undoing convolution broadening" by downweighting large-x
tails.

**Question:** Is this correction derived from the posterior of a product of
independent Poisson PDFs? If so, what is the precise derivation, and what does
`num_component` count — the number of distinct MC samples contributing to the
bin (beam, dirt, EXT, ...) or the number of convolution steps? Is this formula
documented in an internal MicroBooNE note or a conference proceeding that can
be cited?

**Why it matters:** If the correction is correct and documented externally,
add the reference to 07_algorithms.md. If it is an approximation whose validity
range is not established, it should be flagged as a potential source of
mis-estimated MC-stat uncertainties.

---

### OQ-03  [channel layout]  src/TLee.cxx:253-333

**Context:** `Minimization_Lee_strength_FullCov` partitions the covariance
matrix into Y (signal) and X (sideband) blocks using the literal boundaries
`0-7`, `26-33` (signal), `8-25`, `34-51`, `52-136` (sidebands). These same
boundaries appear in the GoF bin-assignment (TLee.cxx:327-333) and are the
basis for the `num_Y = 26+26 = 52` and total `= 137` magic numbers at
TLee.cxx:283, 332, 333 (see also B-15 in 05_bugs.md).

**Question:** Do the channel ranges 0-7 / 26-33 (Y) and 8-25 / 34-51 / 52-136
(X) exactly match the current live analysis layout for all supported run
configurations? The boundary `26` separates FC from PC categories, and `52`
marks the start of the extended sideband region. Have these boundaries changed
since the code was written, or are they still the authoritative layout?

**Why it matters:** If the layout has changed and the literals are stale, the
wrong bins are included in the signal region and the chi-squared is computed
over an incorrect partition. Confirming the layout allows the literals to be
replaced with named constants, closing B-15.

---

### OQ-04  [EM charge scale]  inc/WCPLEEANA/cuts.h:24, 3719-3783

**Context:** The global constant `em_charge_scale = 0.95` (cuts.h:24) is
applied in `get_reco_Enu_corr`, `get_reco_showerKE_corr`, and several branches
of `get_kine_var` for the DATA branch only. It is consistently applied inside
`is_cc_pi0` and `is_pi0` via `get_kine_var`. However, `is_NCpio_sel` at
cuts.h:3766 also uses pi0 energy variables; whether those variables go through
the same `em_charge_scale` path depends on how the caller invokes `get_kine_var`.

**Question:** Is the asymmetric treatment of `em_charge_scale` between data and
MC in these predicates intentional? Specifically, is `is_NCpio_sel` expected to
apply the same 0.95 factor as `is_cc_pi0` and `is_pi0`, or is the NC pi0
channel calibrated separately? And what is the origin of the 0.95 value — is it
from a dedicated EM charge-scale calibration note?

**Why it matters:** If the factor is inconsistently applied across pi0 selection
predicates, events near the pi0 mass window boundary are cut differently between
data and MC in some channels but not others. If the 0.95 is a measured
calibration constant, a citation should be added to 03_selection_layer.md.

---

### OQ-05  [truth variable dispatcher]  inc/WCPLEEANA/cuts.h:270-275

**Context:** `get_truth_var(var_name, eval)` is a one-branch dispatcher: the
only recognized string is `"truth_energyInside"`, which returns
`eval.truth_energyInside`. All other values of `var_name` fall through and
return an uninitialized stack variable (the function has no `else` or default
return path).

**Question:** Is `get_truth_var` intentionally restricted to `"truth_energyInside"`
— i.e., it is only ever called with that argument in the current analysis — or
was it meant to be a general truth-variable dispatcher (analogous to
`get_kine_var`) that was never extended? Are there analysis configurations that
call `get_truth_var` with other variable names?

**Why it matters:** If the function is correctly narrow, a comment should say so
and the missing `else { return 0.; }` default should be added to prevent
undefined behaviour. If it is an incomplete dispatcher, the missing variables are
a silent correctness gap that should be filled.

---

### OQ-06  [Wiener-SVD binning]  src/WienerSVD.cxx:40-50

**Context:** For smoothness-matrix types 22/32 the slice boundaries are the
hardcoded vector `{ 0, 3, 7, 11, 14, 18, 22, 26, 31, 36 }`, and for types
23/33 a 36-element boundary vector is used (WienerSVD.cxx:47-50). A commented-out
alternative block at lines 41-45 has different boundary values, suggesting the
layout changed at least once. See also B-16 in 05_bugs.md.

**Question:** Do the currently-active `dim_edges` arrays exactly match the
binning used in the cross-section analysis that calls these types? The 2D table
appears to encode 9 angular slices with 3-5 momentum bins each (total 36 bins);
the 3D table encodes 35 slices over 138 bins. Are these consistent with the
current analysis configuration files?

**Why it matters:** A mismatch between the hardcoded boundary array and the
actual histogram binning silently applies the wrong smoothness regularisation
matrix, biasing the unfolded cross-section. If the values are confirmed correct
they should be documented as constants with physical labels.

---

### OQ-07  [GP log-scale mutation]  src/GPKernel.cxx:50-52

**Context:** `RBFKernel::Mag(p1, p2)` receives `GPPoint` by value, obtains raw
pointers `p1x = p1.X()` and `p2x = p2.X()`, and writes `p1x[i] = log(p1x[i])`
directly. Because `GPPoint::X()` returns a pointer to the object's internal
array (not a copy), this mutates the local-copy point's coordinates in place
for every call to `Mag`. B-02 in 05_bugs.md documents the incorrect early-exit
guard consequence; the question here is about intent.

**Question:** Was the log-scaling of the GP input space intentional, meaning
that the kernel is designed to operate in log-coordinate space for certain
dimensions? If so, does it matter that the mutation is performed on `p1x/p2x`
(which are pointers into the by-value copy's array) rather than a separate
buffer? Under `GPSmoothing.C`'s usage pattern — calling `Fit` then `Predict`
with the same point set — is the mutation self-consistent, or are there call
sequences where the same `GPPoint` is passed to `Mag` more than once and the
second call sees already-logarithm-transformed coordinates?

**Why it matters:** If intentional and provably self-consistent under the current
call pattern, document that assumption explicitly. If the mutation can produce
doubly-transformed coordinates in any reachable path, the fix (local copy before
log-transform) is required.

---

### OQ-08  [event-matching key]  apps/merge_det.cxx:925,947 vs apps/convert_cv_spec.cxx:519-520

**Context:** `merge_det.cxx` matches events between the CV tree and the
detector-variation tree using a 2-tuple `(run, event)` as the map key
(B-20 in 05_bugs.md). The closely related `convert_cv_spec.cxx` explicitly
switched to a 3-tuple `(run, subrun, event)` for this reason. Both perform the
same logical pairing operation for detector-systematic inputs.

**Question:** Which event-matching key is correct for the current data format?
In the MicroBooNE data files used for detector systematics, are event numbers
guaranteed to be unique within a run (making 2-tuple matching safe), or can
two subruns within the same run legally share event numbers? If 3-tuple matching
is required, is `merge_det.cxx` known to produce incorrect pairings on any
existing input dataset?

**Why it matters:** If 2-tuple matching is incorrect, the detector covariance
matrix is computed from mis-paired CV/detvar events, which is a silent physics
error. If it is safe for all current detector-variation samples, a code comment
should document the guarantee.

---

### OQ-09  [zero-prediction fallback]  apps/det_cov_matrix.cxx:133

**Context:** When the prediction in a bin is zero but the off-diagonal covariance
entry is non-zero, `det_cov_matrix.cxx` sets the fractional covariance diagonal
to `1./16.` (i.e., `(0.25)^2`, corresponding to a 25% relative uncertainty).
The in-line comment reads `// 25% uncertainties ...` with no further explanation.
B-21 in 05_bugs.md flags this as a magic number without physical motivation.

**Question:** What is the physical meaning of assigning 25% uncertainty to a
zero-prediction bin? Is this value derived from a known detector-variation
amplitude at low statistics, taken from a MicroBooNE technical note, or is it
an ad-hoc placeholder from the initial code draft? Is it still considered
appropriate for the current analysis configuration?

**Why it matters:** If there is a documented physical basis, a comment and
citation should be added and B-21 can be closed as intentional. If it is a
placeholder, it should be replaced by a configurable parameter or a
physics-motivated estimate.

---

### OQ-10  [MC-stat TGraph interpolation]  src/TLee.cxx:2257, 2829-2840

**Context:** Inside `Set_Collapse`, the MC-stat contribution to the bin
covariance is added as `gh_mc_stat_bin[ibin]->Eval(scaleF_Lee)`, where
`gh_mc_stat_bin[ibin]` is a `TGraph`. According to 07_algorithms.md §4, this
TGraph was built at startup by scanning over `Lee_strength` values and
computing the Bayesian MC-stat variance at each one (TLee.cxx:2829-2840).

**Question:** Is this TGraph precomputation the mechanism by which the MC-stat
covariance varies as a function of LEE signal strength, and was it generated
once and stored in the MC-stat log directory (`mc_directory/N.log` read by
`Set_Spectra_MatrixCov`) or recomputed fresh on each program run? If recomputed,
at what grid of `scaleF_Lee` values, and is a coarse grid a known source of
interpolation error?

**Why it matters:** If the TGraph was generated externally and stored, then
changes to the LEE signal template require regenerating those log files before
the fit will use the correct MC-stat covariance. If it is recomputed at runtime,
the scan grid spacing should be documented. Either way, this mechanism is not
described in any existing document.

---

### OQ-11  [reweight universe RNG seeding]  src/mcm_2.h:1990,2033

**Context:** B-01 in 05_bugs.md (HIGH) documents that `gRandom->SetSeed(j*reweight*77777)`
produces a near-deterministic lattice of seeds rather than independent draws,
collapsing the reweight covariance toward a rank-1 structure. The bug report
suggests replacing per-universe re-seeding with a single seed before the loop.

**Question:** Was the per-universe seeding introduced to make the covariance
calculation reproducible across different run configurations (so that universe j
always uses the same random numbers regardless of how many total universes are
computed)? If reproducibility per-universe was the intent, a better mechanism
would be a universe-indexed deterministic RNG (e.g., `TRandom3` seeded with a
hash of `systematic_index * N_universes + j`). Please confirm whether
reproducibility per universe was the design goal, as this determines which fix
is appropriate.

**Why it matters:** The fix strategy differs depending on intent: if
reproducibility matters, a deterministic per-universe seed is needed; if not, a
single seed before the loop is sufficient. The current code silently produces
wrong covariances either way.

---

### OQ-12  [rw_type variable]  src/master_cov_matrix.cxx:75, src/mcm_2.h:2026

**Context:** In the `CovMatrix` constructor, `flag_reweight` is set to `true`
whenever any configured reweight entry has `flag_reweight_i==1`, but `rw_type`
is never updated from its default value of 0. Code in `mcm_2.h` that branches
on `rw_type==3` (the `"UBGenieFluxSmallUni"` path) therefore never executes
from a standard configuration.

**Question:** Is the `rw_type==3` branch dead by design — i.e., is
`UBGenieFluxSmallUni` handling superseded by the standard reweight path — or is
`rw_type` intended to be set from configuration but the setter was never wired
in? If it should be settable, what value in `rw_cv_input.txt` would trigger it?

**Why it matters:** If the branch is dead by design, it should be removed or
guarded with `#if 0` to prevent future confusion. If it is reachable by
configuration, the missing setter is a defect that silently disables a
systematic variation.

---

### OQ-13  [TPrincipal for covariance estimation]  src/mcm_1.h:77,139; src/mcm_data_stat.h:69

**Context:** Both `gen_det_cov_matrix` (Stage 1) and `gen_data_stat_cov_matrix`
use ROOT's `TPrincipal` to accumulate bootstrap samples and then extract the
covariance via `prin.GetCovarianceMatrix()`. `TPrincipal` is a PCA class; while
it does internally track the sample covariance, this usage is documented in
02_core_framework.md §5 as "unusual and fragile".

**Question:** Was `TPrincipal` chosen deliberately because its covariance
accumulation is numerically more stable than a naive outer-product sum, or was
it used pragmatically because it was a convenient available class? Are there
known precision issues with retrieving the covariance from `TPrincipal` when the
number of bootstrap rows (1000 or 5000) is comparable to the number of bins
(O(50-100))? Would a direct `(x - mean)(x - mean)^T` accumulator be equivalent?

**Why it matters:** If `TPrincipal` introduces numerical artefacts (e.g.,
centering errors or internal `double` accumulation overflow), the bootstrap
covariance is wrong. If it is simply a convenient equivalent, the code can be
replaced with a simpler and more transparent accumulator as part of refactoring.

---

### OQ-14  [data POT hardcode]  src/master_cov_matrix.cxx:720; src/mcm_1.h:57

**Context:** `double data_pot = 5e19` is the fallback POT used in covariance
normalisation if the run-period lookup fails (B-14 in 05_bugs.md, LOW). The
value 5×10¹⁹ POT is the nominal MicroBooNE Run 1 open-data exposure.

**Question:** Under what circumstances does the run-period lookup fail in
practice? Is the `5e19` fallback ever reached during a normal production run,
or is it only reachable when configuration files are malformed? If it is
reachable in normal use for certain run periods, the silent fallback is a
live normalisation error.

**Why it matters:** If the fallback is never reached in production, a `throw`
or hard abort on lookup failure would be safer. If it is reached for some run
periods, those periods are silently using the wrong POT normalisation for
their covariance, which is a physics error.

---

*End of document. Questions OQ-01 through OQ-14 are the primary set. The
recommended triage order is: OQ-01 (detector cov amplitude), OQ-03 (channel
layout), OQ-08 (event matching), OQ-11 (RNG seeding intent) — these four have
the highest potential physics impact if the answers indicate defects.*
