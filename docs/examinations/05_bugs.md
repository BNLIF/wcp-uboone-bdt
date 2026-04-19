# 05 — Bug Catalogue

This document catalogues confirmed code defects in the wcp-uboone-bdt analysis framework, ordered by potential impact on physics results. Each finding was verified by reading the actual source line(s) before inclusion; line numbers are correct as of the current HEAD. Severity is assigned on a three-level scale: **High** — produces silently wrong physics outputs; **Medium** — latent or conditional, can corrupt results in identifiable circumstances; **Low** — brittle assumptions or hardcoded constants that could silently break on dataset changes; **Style** — naming or dead-code issues with no runtime consequence.

---

## High Severity

### B-01  [HIGH]  src/master_cov_matrix.cxx:1990,1995 (also 2033,2038)

**Symptom:** The 1 000 Gaussian universe throws in the `"reweight"` systematic branch are driven by a nearly deterministic quasi-lattice of seeds rather than independent random draws.

**Explanation:** The call `gRandom->SetSeed(j*reweight*77777)` casts the floating-point product `j * reweight * 77777` to an unsigned integer. Because `reweight` is typically very close to 1.0, the seeds for consecutive universe indices `j` and `j+1` differ by only `~77777`, which is a tiny fraction of the full 32-bit seed space. ROOT's `TRandom3` is seeded to a strongly correlated sub-lattice, so the 1 000 draws are not statistically independent. The resulting covariance matrix for the reweighting systematic is effectively collapsed to a rank-1 (or very-low-rank) structure instead of the intended full-rank Gaussian scatter. The same pattern is repeated inside the `"UBGenieFluxSmallUni"` branch at lines 2033 and 2038.

**Fix sketch:** Replace the per-universe seed with a single seed set once before the loop (e.g. `gRandom->SetSeed(run_number * 131071 + systematic_index)`) and let `gRandom->Gaus()` advance the state naturally across universes. Alternatively, construct a `TRandom3` with seed 0 (clock-based) once and draw from it inside the loop without re-seeding.

**Fixed:** commit `d99fa05` — seed moved before the universe loop using `(unsigned int)(weight.run * 131071u + weight.event)`; both the `"reweight"` and `"UBGenieFluxSmallUni"` branches updated. Covered by `test/test_master_cov_seeds.sh`.

**Cross-reference:** → 02_core_framework.md §CovMatrix / systematic universe construction

---

### B-02  [HIGH]  src/GPKernel.cxx:50–59

**Symptom:** Within a single `RBFKernel::Mag` evaluation, the log-transform of coordinates is applied *before* the zero-length-scale guard check, causing the guard to test log-transformed values instead of the original coordinates.

**Explanation:** Lines 50–52 overwrite `p1x[i]` and `p2x[i]` with their logarithms in place; the guard `if (GPKernel::fPars[i]==0 && p1x[i]!=p2x[i])` at line 59 is then comparing log-transformed coordinates. For a dimension whose scale is zero (disabled), the equality `p1x[i] != p2x[i]` is evaluated on the log values, not the originals, so two points that originally differed only in a log-scale dimension will give the wrong early-exit result of `1e6`. Because `GPPoint` is passed by value, the mutation is confined to local copies and does not corrupt training-set coordinates across calls; the damage is limited to incorrect early-exit logic within the affected call.

**Fix sketch:** Apply the log transformation into a separate local array rather than overwriting `p1x` and `p2x`, so the dimension-zero guard at line 59 still operates on the original coordinate values.

**Fixed:** commit `b24aaa7` — guard moved before the log-transform loop; zero-length-scale dimensions are skipped entirely in the transform (set to 0.0) so `log()` is never called on them. Covered by `test/test_gpkernel.cxx`.

**Cross-reference:** → 02_core_framework.md §GPRegressor / kernel evaluation

---

### B-03  [HIGH]  inc/WCPLEEANA/cuts.h:209,223–228,231

**Symptom:** In `get_weight`, when a reweighting variable falls between no pair of custom bin edges, `wbin` is used uninitialised as a vector index, corrupting the event weight with a stack-garbage value.

**Explanation:** `int wbin;` is declared at line 209 with no initialiser. In the `equal_binning == false` branch the loop at lines 223–228 sets `wbin = b` only when `var <= bins[b+1] && var > bins[b]`; if `var` exactly equals `bins[0]`, or falls below `bins[0]` with `underflow == false` (the code falls through to the `else if(var>min_var)` block only when `var > min_var`), the loop exits without a match. After the loop, line 231 unconditionally executes `addtl_weight *= reweight[wbin]`, indexing `reweight` with an indeterminate stack integer.

**Fix sketch:** Initialise `wbin = 0` (or the appropriate underflow bin) at declaration, and add an explicit `found` flag or `break`-with-sentinel after the loop so that an out-of-range `var` falls back to a safe default weight rather than using garbage as an index.

**Fixed:** commit `7b91d96` — `wbin` initialised to `-1`; reweight access guarded by `if(wbin >= 0 && wbin < (int)reweight.size())`. Covered by `test/test_cuts_wbin.cxx`.

**Cross-reference:** → 03_selection_layer.md §get_weight / custom-binning reweighting

---

### B-04  [HIGH]  src/TLee.cxx:1130–1229

**Symptom:** In `Exe_Goodness_of_fit`, all seven consecutive axis-configuration blocks test `index == 0`, so only the last block's settings (axis range 0–2600 MeV, 26 bins) are ever applied; the intended channel-specific axis ranges for indices 1–6 are dead code.

**Explanation:** The seven `if(index==0){...}` blocks at lines 1130, 1148, 1166, 1179, 1192, 1205, and 1218 are plainly intended to correspond to `index == 0` through `index == 6` (one per goodness-of-fit channel). Because all conditions are identical, each block overwrites the same variables, and the final block — with `userAA_index_hgh = 26` and range 0–2600 MeV — is the one that actually takes effect for every channel. Plots for channels with different kinematic ranges (e.g. the pi-zero channel at lines 1166–1177) are drawn with the wrong axis.

**Fix sketch:** Change the condition of each successive block from `index==0` to `index==1`, `index==2`, ..., `index==6`. Alternatively, use a `switch(index)` statement or a lookup table indexed by channel number to select axis parameters.

**Fixed:** commit `06d1f69` — blocks 2–7 changed to `if(index==1)` … `if(index==6)`. Covered by `test/test_tlee_index.sh`.

**Cross-reference:** → 02_core_framework.md §TLee / goodness-of-fit output

---

## Medium Severity

### B-05  [MEDIUM]  inc/WCPLEEANA/Configure_Lee.h (entire file)

**Symptom:** Every translation unit that includes this header acquires its own copy of all configuration variables, with no guard against multiple inclusion; if two `.cxx` files include it, the linker sees duplicate definitions.

**Explanation:** The file has no `#ifndef` include guard and no `#pragma once`. It defines non-`const` variables such as `TString spectra_file`, `bool flag_display_graphics`, `double Lee_strength_for_outputfile_covariance_matrix`, and many others inside `namespace config_Lee`. In C++, a non-`const` variable definition in a header produces multiple-definition link errors if the header is included in more than one translation unit. In the current build the header appears to be included from a single `.cxx`, but any future refactoring risks an ODR violation.

**Fix sketch:** Add `#pragma once` (or an `#ifndef` guard) at the top of the file. Convert the mutable variables to `extern` declarations in the header and move the definitions to a corresponding `.cxx` file, or make the constants `inline constexpr`.

**Cross-reference:** → 02_core_framework.md §TLee configuration

---

### B-06  [MEDIUM]  src/bayes.cxx:30–32,38

**Symptom:** The `f_conv_vec` cleanup loop is commented out while `f_conv` (which points to the last element of `f_conv_vec`) is separately deleted, risking a double-delete or use-after-free when convolution builds a chain of more than one component.

**Explanation:** In `add_meas_component` (line 151), every new convolution TF1 is appended to `f_conv_vec` and `f_conv` is updated to point to the most recently created element. The destructor (lines 30–32) has the entire `f_conv_vec` deletion loop commented out, but line 38 still calls `delete f_conv`. The pointer `f_conv` is therefore the only remaining reference deleted at destruction; the intermediate entries in `f_conv_vec` (all except the last) leak memory. If the comment were re-enabled, the final element would be deleted twice.

**Fix sketch:** Either uncomment and rely solely on the `f_conv_vec` deletion loop (and set `f_conv` to `nullptr` before the destructor finishes so the separate guard at line 38 becomes a no-op), or stop storing intermediate TF1 objects in `f_conv_vec` at all and rely exclusively on the `conv_vec` + single `f_conv` cleanup path.

**Cross-reference:** → 02_core_framework.md §Bayes posterior calculation

---

### B-07  [MEDIUM]  src/GPRegressor.cxx:191–201

**Symptom:** `MarginalLikelihood::DoDerivative` uses `fAlpha` and `fKInv` without calling `Solve(par)` first, so if the optimiser evaluates the gradient at a parameter point not preceded by an `DoEval` call, the derivative is computed with a stale kernel matrix.

**Explanation:** `DoDerivative` at line 191 builds `dK` from the current kernel parameters but then immediately uses `fAlpha` and `fKInv`, which are member fields set only by `Solve`. ROOT's `ROOT::Fit::Fitter` (used in `SolveHyperParameters`) typically calls `DoEval` before `DoDerivative` in the same parameter point, but this is a convention, not an API contract. If the optimiser (GSLMultiMin/BFGS2) evaluates gradient-only steps at a new parameter point without a preceding function evaluation, the result is wrong.

**Fix sketch:** Add an explicit `Solve(par)` call at the top of `DoDerivative`, or assert that the current parameter vector matches the one used in the last `Solve` call. Since `Solve` is cheap relative to kernel matrix inversion the simplest fix is unconditional re-solve.

**Cross-reference:** → 02_core_framework.md §GPRegressor hyper-parameter optimisation

---

### B-08  [MEDIUM]  src/GPRegressor.cxx:71

**Symptom:** `new double[n]` allocated for the initial parameter array `par` in `SolveHyperParameters` is never freed, leaking memory every time hyper-parameter optimisation is run.

**Explanation:** At line 71, `double* par = new double[n]` is allocated, filled, passed to `fitter.SetFCN`, and then the function returns without a `delete[] par`. The fitter copies the values internally, so the pointer is genuinely orphaned. For a typical workflow running hyper-parameter optimisation once this is a small fixed-size leak, but it is still a resource error.

**Fix sketch:** Replace `new double[n]` with `std::vector<double> par(n)` and pass `par.data()` to `SetFCN`, which eliminates the manual allocation entirely.

**Cross-reference:** → 02_core_framework.md §GPRegressor

---

### B-09  [MEDIUM]  src/WienerSVD.cxx:125–127 and 169

**Symptom:** Neither `TDecompSVD` call (for the covariance matrix at line 125 and for the response matrix at line 169) checks its return value, so a failed decomposition silently propagates zeroed or garbage singular vectors into the unfolding result.

**Explanation:** `TDecompSVD` can fail (rank-deficient or ill-conditioned input) and signals this through its return value from `Decompose()`. The code calls `decV.GetSig()`, `udv.GetU()`, etc., unconditionally. ROOT's lazy decomposition means these getters implicitly trigger decomposition; if it fails the returned matrices contain undefined values which then silently propagate to the unfolded spectrum.

**Fix sketch:** Call `decV.Decompose()` (and `udv.Decompose()`) explicitly, check the `bool` return, and throw an exception or print a meaningful error and return a sentinel if decomposition fails.

**Cross-reference:** → 04_apps_pipeline.md §WienerSVD unfolding

---

### B-10  [MEDIUM]  apps/merge_hist.cxx:86

**Symptom:** Each iteration of the outer file loop opens a new `TFile` via `new TFile(out_filename)` without closing or deleting the previous one, so all input files remain open simultaneously for the entire duration of the loop.

**Explanation:** The variable `temp_file` is overwritten at every iteration (line 86) with a new `TFile*`. The previous pointer is orphaned: no `Close()` or `delete` is called. ROOT's global file-management list keeps all files open until the process ends. For an analysis run with many input files this exhausts file descriptors and inflates ROOT's internal memory.

**Fix sketch:** Add `if (temp_file) { temp_file->Close(); delete temp_file; }` before the `new TFile(...)` assignment, or use a RAII wrapper (`std::unique_ptr<TFile>` with a custom deleter).

**Cross-reference:** → 04_apps_pipeline.md §merge_hist

---

### B-11  [MEDIUM]  inc/WCPLEEANA/kine.h:46–49

**Symptom:** `clear_kine_info` calls `->clear()` on four raw pointer members (`kine_energy_particle`, `kine_energy_info`, `kine_particle_type`, `kine_energy_included`) without null-checking them first.

**Explanation:** The `KineInfo` struct stores these as raw `std::vector<...>*` pointers. `clear_kine_info` is called before each tree entry is read; if a `KineInfo` instance is ever used before its branch addresses are wired (or if a branch is absent in a particular file), the pointers are uninitialised, and `->clear()` is undefined behaviour.

**Fix sketch:** Add null-pointer guards (`if (ptr) ptr->clear();`) around each pointer dereference, or initialise all pointer members to `nullptr` in a constructor or `init_pointers` helper.

**Cross-reference:** → 03_selection_layer.md §KineInfo branch setup

---

### B-12  [MEDIUM]  src/bayes.cxx:229–282

**Symptom:** When the bisection loop for the credible-interval bound fails to converge within 30 iterations, the code emits a `std::cerr` message and `break`s, then silently uses whatever half-interval boundary it had reached as the result.

**Explanation:** The loop at lines 229 and 241 exits via `break` when `line > 30`, but the caller receives the value of `(val_typeA + val_typeB)/2` computed just before the break, which may not satisfy the required precision `analytic_result_precision`. No flag, exception, or return-code signals the failure to the caller.

**Fix sketch:** Set a `bool converged = false` flag, set it to `true` on normal exit, and check it after the loop to either throw a `std::runtime_error` or return `NaN`/sentinel so the caller can handle non-convergence explicitly.

**Cross-reference:** → 02_core_framework.md §Bayes credible interval

---

### B-13  [MEDIUM]  apps/bdt_convert.cxx:91

**Symptom:** The training-list file is read with `while(!infile.eof())`, which causes the last record to be processed twice when the final `>>` reads exactly to EOF.

**Explanation:** `std::istream::eof()` is set only *after* a read that hits the end-of-file, not before it. The standard anti-pattern `while(!infile.eof())` therefore enters the loop one extra time after the last successful read; the failed `infile >> tmp_type >> run >> subrun` leaves the variables unchanged from the previous iteration, so that record is inserted again, potentially inflating the training/test set with a duplicate entry.

**Fix sketch:** Replace `while(!infile.eof())` with `while(infile >> tmp_type >> run >> subrun)`, which checks the stream state after each read and exits cleanly at EOF.

**Cross-reference:** → 04_apps_pipeline.md §bdt_convert

---

## Low Severity

### B-14  [LOW]  src/master_cov_matrix.cxx:720; src/mcm_1.h:57 (comment)

**Symptom:** The data POT is hardcoded as `5e19` in the covariance matrix code, so any run with a different exposure silently uses the wrong normalisation unless the constant is manually updated.

**Explanation:** `double data_pot = 5e19;` at line 720 of `master_cov_matrix.cxx` and the analogous comment in `mcm_1.h:57` pin the data exposure to the MicroBooNE Run 1 open-data value. If a different run's data is passed in, the POT used to normalise MC predictions is wrong, causing a systematic offset in all covariance estimates.

**Fix sketch:** Read the POT from the data input file (which already stores it in the `T` tree's `pot` branch, as seen in `merge_hist.cxx:88–90`) and propagate it through the covariance-matrix functions rather than using a literal.

**Cross-reference:** → 02_core_framework.md §CovMatrix POT normalisation

---

### B-15  [LOW]  src/TLee.cxx:283,332,333

**Symptom:** The channel-count literals `26+26` and `137` for the nueCC FC+PC and total-channel boundaries are hardcoded in `Minimization_Lee_strength_FullCov`, so adding or removing analysis channels silently gives wrong covariance submatrix slices.

**Explanation:** The code partitions the full covariance matrix using `int num_Y = 26+26` (52 nueCC bins) and `137` (total bins including side-bands) at lines 283, 332, and 333. These are arithmetic on the channel layout documented nowhere nearby. Any change to the number of analysis bins (e.g. adding a new sideband channel) requires manual updates in at least three places.

**Fix sketch:** Replace the literals with named constants or compute the values from the actual matrix dimension and the configuration in `Configure_Lee.h`. At minimum, `static_assert` or a runtime check should verify that the hardcoded sum matches the matrix dimensionality.

**Cross-reference:** → 02_core_framework.md §TLee channel layout

---

### B-16  [LOW]  src/WienerSVD.cxx:40–50

**Symptom:** The `dim_edges` bin-boundary arrays for the 2D and 3D regularisation matrix types are hardcoded integers that encode the specific analysis channel layout; a change to histogram binning silently misaligns the smoothness matrix.

**Explanation:** The vectors `{ 0, 3, 7, 11, 14, 18, 22, 26, 31, 36 }` (type 22/32) and the 36-element vector for type 23/33 (lines 40–50) encode the exact bin-count structure of the multi-dimensional unfolding problem. These must match the `TMatrixD` dimensions passed at call time but there is no runtime check. A commented-out alternative block (lines 41–45) with different values suggests the layout has changed at least once without a corresponding refactoring.

**Fix sketch:** Pass the bin-edge structure as a parameter to `Matrix_C` or derive it from the input matrix dimensions, and add an assertion that `dim_edges.back() == n`.

**Cross-reference:** → 04_apps_pipeline.md §WienerSVD unfolding

---

### B-17  [LOW]  apps/bdt_convert.cxx; apps/filter_goodruns.cxx; apps/numi_filter.cxx

**Symptom:** The good-run list is copy-pasted as a large `good_run_list_vec` initialiser in all three application source files; updating the list requires editing three independent locations.

**Explanation:** All three files contain an identical (or near-identical) multi-hundred-element `std::vector<int>` initialiser for `good_run_list_vec`. There is no shared constant or configuration file. A run added to one list but missed in another produces inconsistent event selections between analysis steps.

**Fix sketch:** Extract the good-run list into a shared header or a plain-text configuration file read at runtime (a format already used elsewhere in the framework for other configuration), and have all three applications read from the single source.

**Cross-reference:** → 04_apps_pipeline.md §good-run filtering

---

### B-18  [LOW]  inc/WCPLEEANA/cuts.h:3910; inc/WCPLEEANA/tagger.h:573

**Symptom:** `is_numuCC_cutbased` compares `tagger_info.cosmict_flag == 0` where `cosmict_flag` is declared as `float`, which is an exact floating-point equality test against an integer literal.

**Explanation:** `TaggerInfo::cosmict_flag` is of type `float` (tagger.h line 573). The cut at cuts.h line 3910 tests `tagger_info.cosmict_flag == 0`. Although the BDT output is stored as 0.0 or 1.0 and can be represented exactly in IEEE 754 single precision, the pattern is fragile: any arithmetic along the data chain that produces a value slightly different from exact 0.0 (e.g. type conversion via double followed by narrowing) would silently pass or fail the cut. The same pattern occurs wherever `numu_cc_flag == 1` is tested on a field declared `float`.

**Fix sketch:** Convert the relevant flag fields to `int` or `bool`, or replace floating-point equality tests with a small-epsilon comparison or an explicit cast, e.g. `static_cast<int>(tagger_info.cosmict_flag) == 0`.

**Cross-reference:** → 03_selection_layer.md §numuCC cut-based selection

---

### B-19  [LOW]  apps/xs_cov_matrix.cxx:192–193

**Symptom:** The output ROOT key for the cross-section covariance matrix is `"cov_xf_mat_N"` (with `xf` standing for "flux") rather than `"cov_xs_mat_N"`, so any downstream code looking for `"cov_xs_mat_N"` will not find this object.

**Explanation:** At line 192, `cov_xs_mat->Write(Form("cov_xf_mat_%d",run))` and `frac_cov_xs_mat->Write(Form("frac_cov_xf_mat_%d",run))` at line 193 use the key suffix `xf` (matching the flux-covariance naming convention) rather than `xs`. Variable and comment context surrounding the block make clear this file is the cross-section covariance output. If downstream analysis reads `"cov_xs_mat_%d"` it will silently fail to find the matrix.

**Fix sketch:** Change the format string from `"cov_xf_mat_%d"` to `"cov_xs_mat_%d"` (and similarly for `frac_cov_xf_mat`) to match the cross-section context, and verify all downstream readers use the same key.

**Cross-reference:** → 04_apps_pipeline.md §xs_cov_matrix

---

### B-20  [LOW]  apps/merge_det.cxx:925,947

**Symptom:** Event matching between the CV and detector-variation trees uses a 2-tuple `(run, event)` as the map key, ignoring `subrun`, which causes incorrect matches when the same run contains multiple subruns with the same event number.

**Explanation:** At lines 925 and 947, `map_re_entry_cv` and `map_re_entry_det` are keyed on `std::pair<int,int>{run, event}`. In MicroBooNE data, event numbers are unique only within a (run, subrun) pair; two subruns in the same run can legally share an event number. The adjacent `map_rs_re_cv` (line 913) correctly uses `(run, subrun)` → `{(run, event)}` nesting, but the flat lookup map that drives the actual matching drops the subrun dimension.

**Fix sketch:** Change the key type of `map_re_entry_cv` and `map_re_entry_det` to `std::tuple<int,int,int>` keyed on `(run, subrun, event)` and update the fill and lookup sites accordingly.

**Cross-reference:** → 04_apps_pipeline.md §merge_det event matching

---

### B-21  [LOW]  apps/det_cov_matrix.cxx:133

**Symptom:** When a prediction bin is zero while the off-diagonal covariance element is non-zero, the fractional covariance diagonal is set to the hardcoded value `1./16.` (25% relative uncertainty squared), regardless of the actual uncertainty size.

**Explanation:** At line 133, the comment reads `// 25% uncertainties ...` and sets `(*frac_cov_det_mat)(i,j) = 1./16.` for the `i==j` case. This is a placeholder that assigns a fixed 25% fractional uncertainty to any zero-prediction bin with a non-zero covariance entry, without physical motivation or traceability to the actual detector variation amplitude.

**Fix sketch:** Replace the literal with either a configurable parameter or a physics-motivated estimate derived from the non-zero neighbouring bins. At minimum, document in a comment why 25% is appropriate for these bins and under which run conditions.

**Cross-reference:** → 04_apps_pipeline.md §det_cov_matrix

---

## Style

### B-22  [STYLE]  src/Util.cxx:29 vs inc/WCPLEEANA/Util.h:17

**Symptom:** The function is declared as `MatrixMatrix` in the header but defined as `MatrixMatirx` (transposed `ri`) in the implementation; the linker resolves both because only one object file is involved, but the mismatch prevents any external caller of the header declaration from linking.

**Explanation:** `Util.h` line 17 declares `void MatrixMatrix(TMatrixD M1, TMatrixD M2)`. `Util.cxx` line 29 defines `void MatrixMatirx(...)`. These are two different symbols: the declared but never-defined `MatrixMatrix` is an unresolved external for any user, and the defined `MatrixMatirx` is unreachable via the header API.

**Fix sketch:** Correct the typo in `Util.cxx` to `MatrixMatrix`.

**Cross-reference:** → 02_core_framework.md §Util helpers

---

### B-23  [STYLE]  src/TLee.cxx:44 and 98

**Symptom:** Two of the three Feldman-Cousins method implementations are named `Exe_Fiedman_Cousins_Data` and `Exe_Fledman_Cousins_Asimov` (misspellings of "Feldman"), while the third is correctly named `Exe_Feldman_Cousins`.

**Explanation:** `Exe_Fiedman_Cousins_Data` (line 44) transposes `el` → `ie`, and `Exe_Fledman_Cousins_Asimov` (line 98) drops the `a`. These are callable entry points whose names appear in external scripts or Makefile targets; any caller using the correct spelling "Feldman" will fail to link.

**Fix sketch:** Rename both methods to `Exe_Feldman_Cousins_Data` and `Exe_Feldman_Cousins_Asimov` respectively, updating all call sites.

**Cross-reference:** → 02_core_framework.md §TLee Feldman-Cousins

---

### B-24  [STYLE]  src/TLee.cxx:2331–2337

**Symptom:** `TLee::Set_TransformMatrix()` is defined with an empty body (only a `cout` debug line), so the feature it was meant to provide is absent.

**Explanation:** The function body at lines 2331–2337 contains only a status print and a comment referencing `Set_Spectra_MatrixCov`. The implementation was never completed. Any caller of `Set_TransformMatrix` obtains no transformation, potentially causing silent no-ops in the analysis chain if the function is ever wired in.

**Fix sketch:** Either implement the function or mark it `[[deprecated]]` / remove it and its declaration, replacing any call sites with a direct call to the equivalent logic in `Set_Spectra_MatrixCov`.

**Cross-reference:** → 02_core_framework.md §TLee matrix setup

---

### B-25  [STYLE]  src/WienerSVD.cxx:190

**Symptom:** `bool test_wiener = false;` is permanently hardcoded, silently disabling a debugging/validation code path.

**Explanation:** The variable `test_wiener` (line 190) controls an alternate Wiener-filter computation path but is hard-set to `false` with no compile-time or runtime switch. The dead code branch is never exercised, making it effectively unmaintained.

**Fix sketch:** Either remove the `test_wiener` branch entirely (if no longer needed) or expose it as a function parameter or preprocessor flag so it can be activated without source edits.

**Cross-reference:** → 04_apps_pipeline.md §WienerSVD

---

### B-26  [STYLE]  inc/WCPLEEANA/bdt.h

**Symptom:** `bdt.h` has no include guard and no `#pragma once`, though its content is declarations only (function prototypes in a namespace), so the ODR risk is currently low.

**Explanation:** Multiple inclusion of `bdt.h` would produce redeclaration errors for the function prototypes inside `namespace LEEana`. The file is currently included in a controlled way, but the absence of a guard is a maintenance hazard.

**Fix sketch:** Add `#pragma once` at the top of the file.

**Cross-reference:** → 03_selection_layer.md §BDT score evaluation
