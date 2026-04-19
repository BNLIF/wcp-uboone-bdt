# 06 — Efficiency and Performance Notes

This document catalogues observed performance and resource-use inefficiencies in the wcp-uboone-bdt framework. The focus is on patterns whose cost multiplies with data size, number of toys, or number of minimiser iterations rather than one-time startup overhead. Each entry was verified against the source before inclusion. "Fix sketch" sections describe approach only; no actual patches are given.

---

### E-01  inc/WCPLEEANA/cuts.h:1635–1720 (`get_cut_pass`)

**Deferred:** The lazy-eval rewrite needs a benchmark confirming the saving outweighs the new per-event struct allocation; per-event hot path must not regress. Requires a timing fixture.

**Impact:** per-event — rebuilds an `O(30)` map on every call; with ~10 M events in a full run the map construction dominates the per-event cost of the selection loop.

**Observation:** `get_cut_pass` constructs a fresh `std::map<std::string, bool> map_cuts_flag` on every invocation (line 1652 onward). Each entry requires a heap-allocated tree node and a string comparison for insertion. The map is then queried by the `ch_name` and `add_cut` strings with further map lookups. The ~30 entries correspond to truth-level flags whose values could be computed lazily only when the relevant channel is queried.

**Fix sketch:** Compute the full flag set once per event entry and cache it in a small stack-allocated struct or a flat `bool` array indexed by an enum, avoiding heap allocation entirely. Alternatively, evaluate only the flags required by the requested `ch_name` string using a direct `if`/`else if` chain rather than a map insert-then-lookup.

**Cross-reference:** → 03_selection_layer.md §get_cut_pass

---

### E-02  inc/WCPLEEANA/cuts.h:238–263 (`get_weight` string compare chain)

**Fixed:** commit f027a3d — replaced the 14-branch TString if/else chain with a static `unordered_map<string,int>` dispatch + switch. O(14) string compares → O(1) hash lookup per call. Test: `test/test_get_weight.cxx` (16 doctest cases).

**Impact:** per-event — called at least once per event per systematic variation; a linear scan through ~14 string comparisons for every weight evaluation.

**Observation:** `get_weight` dispatches on the `weight_name` argument using a chain of `else if (weight_name == "cv_spline")` ... `else if (weight_name == "add_weight")` comparisons (lines 238–263). `TString` equality is a character-by-character comparison. With 14 branches, most events pay for 7–14 comparisons. When called inside a 1 000-universe loop the total string comparison work per event is O(14 000) character operations.

**Fix sketch:** Map `weight_name` to an integer enum once at job initialisation (using an `std::unordered_map<std::string,int>`) and switch on the integer inside the hot loop. Alternatively, accept a typed enum directly as the parameter instead of a `TString`.

**Cross-reference:** → 03_selection_layer.md §get_weight

---

### E-03  src/bayes.cxx:132–154 (`Bayes::add_meas_component` / `do_convolution`)

**Deferred:** Lowering `NPX`/`NofPointsFFT` requires an accuracy budget for credible-interval error < 0.1%; needs physics-owner sign-off before changing.

**Impact:** per-bin in the posterior credible-interval computation — each new measurement component adds a TF1 evaluated at `NPX=60000` points and convolved via FFT with `SetNofPointsFFT(10000)`.

**Observation:** Every call to `add_meas_component` creates a TF1 with `SetNpx(60000)` (line 134, 140, 148) and wraps it in a `TF1Convolution` with `SetNofPointsFFT(10000)` (line 144). A convolution with 10 000 FFT points per TF1 call means that for `N` measurement components the total FFT work is `O(N * 10000 * log(10000))` per bin per credible-interval evaluation. For bins with many Poisson components this is the dominant computational cost.

**Fix sketch:** Profile to find the minimum `NPX` and `NofPointsFFT` that preserve < 0.1% error on the credible interval for the relevant signal rates. Reducing `NPX` to 5 000–10 000 and the FFT points to 2 000–4 000 typically suffices for the Gaussian-like tails of the Poisson-convolution and gives a 4–10x speedup.

**Cross-reference:** → 02_core_framework.md §Bayes posterior calculation

---

### E-04  src/TLee.cxx:499–529 (`TLee::Set_Variations`)

**Deferred:** Parallelising the toy loop requires a thread-safety audit of `TLee` shared state and ROOT IMT compatibility. `std::async`/OpenMP approach is sound but needs validation with actual multi-channel samples.

**Impact:** per-toy setup — performs a full eigendecomposition of the `bins_newworld × bins_newworld` covariance matrix once before toy generation; this is acceptable, but then generates `num_toy` throws inside the same function call with no parallelism.

**Observation:** `Set_Variations` performs a `TMatrixDSymEigen` decomposition (line 511) of an `O(100×100)` symmetric matrix, which is an `O(n^3)` operation executed once per call. The subsequent toy loop (lines 515–531) generates `num_toy` (typically 1 000–10 000) correlated Gaussian throws as matrix-vector products. The matrix product `matrix_eigenvector * matrix_element` (line 525) is `O(n^2)` per toy, giving total throw cost `O(n^2 * num_toy)`. For `n=137` and `num_toy=10000` this is ~188 M floating-point operations, all single-threaded.

**Fix sketch:** The throw loop is embarrassingly parallel; parallelise with `std::async` or OpenMP. If the framework is run many times with the same covariance matrix (e.g. in a scan over LEE strength), cache the eigenvectors and skip the decomposition on repeated calls.

**Cross-reference:** → 02_core_framework.md §TLee toy generation

---

### E-05  src/TLee.cxx:247 (`TLee::Set_Collapse` inside Minuit2 FCN)

**Deferred:** Caching the systematic covariance + rank-1 update requires restructuring `Minimization_Lee_strength_FullCov`; large refactor touching the fit core. Needs careful validation against a reference fit result.

**Impact:** per-minimiser-call — `Set_Collapse` rebuilds the collapsed prediction vector and covariance matrix on every function evaluation inside the Minuit2 minimisation loop.

**Observation:** The Minuit2 lambda FCN at line 232 calls `Set_Collapse()` at line 247 on every function evaluation. `Set_Collapse` (defined at line 2218) applies the collapse matrix to the full prediction and covariance, an `O(n_full^2)` operation. With Minuit2/MIGRAD typically requiring hundreds to thousands of function evaluations, and the full covariance being `O(137×137)`, this rebuilds ~18 000-element matrix products per evaluation.

**Fix sketch:** Factor `Set_Collapse` to separate the parts that depend on `scaleF_Lee` (the only parameter being varied) from the parts that are constant. Only the LEE-channel prediction changes with the scale factor; the collapse matrix and sideband portions can be cached and updated with a rank-1 correction proportional to `Lee_strength`.

**Cross-reference:** → 02_core_framework.md §TLee minimisation

---

### E-06  apps/merge_hist.cxx:80–113 and 128–171 (double pass over input files)

**Deferred (cosmetic only):** The second loop (lines 128–171) does not reopen TFiles — it only iterates the already-loaded `map_name_histogram`. The FD-exhaustion root cause was fixed by B-10 (Wave 2). Merging the loops would be a cosmetic restructuring with <1% speedup.

**Impact:** per-file — the map of input files is iterated twice in sequence; the first pass (lines 80–113) reads POT and histogram metadata, and the second pass (lines 128–171) re-reads histograms already retrieved in the first pass.

**Observation:** Both loops iterate over `map_inputfile_info` (the same set of files) and both read histogram objects from `temp_file`. The first loop stores histogram pointers in `map_name_histogram`; the second loop retrieves them again from the same map for the `filetype==5` data histograms and clones them. The TFile opened in the first pass is never closed before the second pass begins, so all files remain open simultaneously throughout both passes.

**Fix sketch:** Merge both loops into a single pass over `map_inputfile_info`, performing all histogram retrieval and clone operations in one traversal. Close each TFile immediately after extracting all needed data from it.

**Cross-reference:** → 04_apps_pipeline.md §merge_hist

---

### E-07  inc/WCPLEEANA/tagger.h; inc/WCPLEEANA/eval.h; inc/WCPLEEANA/pfeval.h; inc/WCPLEEANA/kine.h; inc/WCPLEEANA/weights.h (`SetBranchAddress` setup)

**Deferred:** Selective branch-disable requires enumerating every address registered by `set_tree_address` across all five headers and adding a `SetBranchStatus("X",1)` call for each. The branch sets are broad (eval.h alone touches flash_*, match_*, truth_*, pl_*, gl_* sub-groups); a partial selective-disable risks silently zeroing un-listed branches at GetEntry time. Needs a per-analysis enable-list audit before changing.

**Impact:** per-file (tree open) — 859 `SetBranchAddress` calls are performed at tree setup time across the five branch-address headers, with `tagger.h` alone contributing 605.

**Observation:** The branch-address setup headers define inline functions that call `tree->SetBranchAddress(name, &field)` for every single field in the struct, including many fields that are never read by the analysis (e.g. the full complement of lol/cosmict/numu sub-scores). ROOT activates all branches for reading when `SetBranchAddress` is called, increasing the per-entry I/O cost proportionally to the number of active branches.

**Fix sketch:** Call `tree->SetBranchStatus("*", 0)` followed by `tree->SetBranchStatus(name, 1)` only for branches actually read by the analysis path being executed. Alternatively, restructure the branch-address setup into per-group enable functions so that, e.g., a cross-section analysis enables only the kinematic branches and not the full cosmic-tagger suite.

**Cross-reference:** → 03_selection_layer.md §branch address setup; → 04_apps_pipeline.md §bdt_convert tree reading

---

### E-08  src/mcm_1.h:77–188 (`gen_det_cov_matrix` two-stage bootstrap + amplification)

**Deferred:** Pre-caching per-file histograms before the bootstrap loop is the right fix; requires profiling to confirm stage-1 event-loop dominates and to determine whether the 16 000 amplification throws can be safely reduced.

**Impact:** per-file (covariance matrix computation) — a two-stage Monte Carlo is used: 1 000 bootstrap throws to estimate the mean detector variation, followed by 16 000 amplification throws from the eigenvectors of that estimate. The full bootstrap stage re-fills histograms on every throw.

**Observation:** Stage 1 (lines 77–141) performs 1 000 bootstrap resamples. Each throw calls `fill_det_histograms` (line 85) and then loops over all covariance channels to accumulate the prediction vector `x[i]`. The `fill_det_histograms` call itself iterates over all events for each detector-variation file on every throw, making the total stage-1 cost `O(1000 * N_events * N_files)`. Stage 2 (lines 167–188) draws 16 000 Gaussian throws from the stage-1 covariance eigenvectors; this is purely algebraic (`O(16000 * n^2)` where `n` is the number of bins) and is much cheaper.

**Fix sketch:** For stage 1, cache the per-file histogram contents before entering the bootstrap loop so that each throw only needs to rescale and add pre-filled histograms rather than re-reading event data. For stage 2, the 16 000 draws are larger than typically needed for a covariance estimate with `n ~ 50–100` bins; profile to determine the minimum number of throws that stabilises the matrix elements to < 1% and reduce accordingly.

**Cross-reference:** → 02_core_framework.md §CovMatrix detector variation

---

### E-09  apps/merge_hist.cxx:86; apps/xs_cov_matrix.cxx:190; apps/det_cov_matrix.cxx:147 (`new TFile` inside loops without reuse)

**Deferred (already mitigated):** The primary `merge_hist.cxx` FD accumulation was fixed by B-10 (Wave 2) — `temp_file` is now closed and deleted after each iteration. The `xs_cov_matrix`/`det_cov_matrix` single-file patterns are low priority. RAII wrapping remains a clean-up for a later pass.

**Impact:** per-file — each iteration of the outer loop opens a new ROOT TFile but does not close the previous one, holding all files open simultaneously.

**Observation:** In `merge_hist.cxx` line 86, `temp_file = new TFile(out_filename)` is reassigned each iteration without `Close()`/`delete` on the old pointer (see also B-10 in 05_bugs.md). In `xs_cov_matrix.cxx` line 190 and `det_cov_matrix.cxx` line 147 single-file output is fine, but the pattern of creating `TFile` objects on the heap with raw `new` and no corresponding `delete` is repeated throughout the application layer, relying on ROOT's global file ownership to avoid leaks at process exit.

**Fix sketch:** Use RAII — wrap each `TFile*` in an `std::unique_ptr<TFile>` with a custom deleter that calls `Close()` before `delete`. This guarantees timely file closure and eliminates the file-descriptor accumulation without requiring manual `Close()` tracking.

**Cross-reference:** → 04_apps_pipeline.md §output file management
