# Apps Pipeline Examination

## 1. Build System

The top-level `wscript` is 31 lines and uses the `wcb` (Wire-Cell Build) framework, which wraps Waf. The `build()` function calls `bld.smplpkg('WCPuBooNE_BDT_APP', use=use)`, which by convention globs `apps/*.cxx` for source files. The glob pattern is flat: it matches only files directly in `apps/`, so `apps/old/*.cxx` is never compiled. There is no explicit exclusion rule; the nesting of `old/` under a subdirectory is what keeps those files out of the build.

The `HAVE_VLNEVAL` conditional (lines 13-14 and 27-28) gates two things: a Boost `program_options` library check in `configure()`, and the addition of `VLNEVAL` and `BOOST` to the dependency list in `build()`. This allows the `eval_vlne` app (which uses `vlneval/zoo/VLNEnergyModel.h` and `boost::program_options`) to be compiled only when the vlneval deep-learning energy estimator library is present. All other apps build unconditionally against ROOT.

---

## 2. Pipeline Overview

```
Raw WCP output ROOT files (wcpselection/ directory with 5 TTrees:
  T_BDTvars, T_eval, T_pot, T_PFeval, T_KINEvars)
         |
         | [per-file input .root]
         v
    bdt_convert
    (applies TMVA BDT scores, filters good runs via hardcoded
     good_run_list_vec and low_lifetime_runs, applies cuts)
         |
         | [condensed checkout .root, same 5 TTrees under wcpselection/]
         v
 +-----------------+---------------------+
 |                 |                     |
 convert_checkout  convert_checkout      merge_det
 _hist             _hist_xs              (pairs CV + det-sys checkout trees
 (standard LEE     (xs-measurement       by run+event number matching)
  channels)         channels)                  |
 |                 |                           | [merged det .root with
 | [per-file       | [per-file                 |  T_eval_cv + T_eval_det etc.]
 |  histograms]    |  histograms]              v
 v                 v                     det_cov_matrix
 merge_hist        merge_hist_xs         (detector covariance ROOT file
 |                 |                      vec_mean, cov_det_mat,
 | [merge.root]    | [merge_xs.root]      frac_cov_det_mat)
 v                 v
 Combined          Combined xs
 histogram         histogram
 file              file
         |
         | [merge.root + cov ROOT files from all cov-matrix apps]
         |
   xf_cov_matrix ----> cov_xf_mat ROOT file (flux+xsec systs)
   xs_cov_matrix ----> cov_xs.root (cross-section systs)
   stat_cov_matrix --> run_data_stat.root (data stat covariance)
   stat_pred_cov_matrix --> MC stat covariance ROOT file
         |
         v
    read_TLee_v20
    (reads spectra + all covariance files via Configure_Lee.h / TLee.h,
     collapses to observation channels, runs chi2 minimization,
     performs GoF tests, optional FC scan)
         |
         v
    file_collapsed_covariance_matrix.root
    (TMatrices for all covariance components, pred/data spectra,
     tree with GoF values; optional file_out_NNN.root for FC toys)
```

---

## 3. Complete App Triage Table

| App name | Category | Pipeline step | Core class/header used | Lines |
|---|---|---|---|---|
| applyNuMIGeomtryWeights | converter | pre-checkout: apply NuMI geometry weights to MC | WCPLEEANA/eval.h | 592 |
| bdt_convert | filter/converter | raw input -> condensed checkout tree | WCPLEEANA/bdt.h, eval.h, tagger.h, TMVA::Reader | 3242 |
| check_failures | util | diagnose failed subrun jobs | WCPLEEANA/eval.h | 790 |
| check_xf_weight_xs | util | inspect flux/xs weight branches | WCPLEEANA/eval.h | 61 |
| CNN_event_info_filter | filter | filter events by CNN score | WCPLEEANA/eval.h, tagger.h | 437 |
| convert_checkout_hist | converter | checkout tree -> per-file histograms (standard) | WCPLEEANA/master_cov_matrix.h, cuts.h | 366 |
| convert_checkout_hist_xs | converter | checkout tree -> per-file histograms (xs mode) | WCPLEEANA/master_cov_matrix.h, cuts.h | 390 |
| convert_cv | converter | CV checkout tree -> merged paired tree (run+event key) | WCPLEEANA/eval.h, tagger.h | 643 |
| convert_cv_spec | converter | CV checkout tree -> merged paired tree (run+subrun+event key) | WCPLEEANA/eval.h, tagger.h | 642 |
| det_cov_matrix | cov-matrix | generate detector systematic covariance | WCPLEEANA/master_cov_matrix.h | 166 |
| eval_vlne | converter | apply deep-learning energy estimator (VLNEVAL) | vlneval/zoo/VLNEnergyModel.h, boost | 530 |
| event_info_filter | filter | general-purpose event selection filter | WCPLEEANA/eval.h, tagger.h | 431 |
| filter_goodruns | filter | standalone good-run filter for checkout trees | WCPLEEANA/eval.h (implicit) | 686 |
| gen_training_list | util | generate BDT training file lists | (stdlib only) | 64 |
| merge_det | converter | pair CV + detector-systematic checkout trees by event ID | WCPLEEANA/eval.h, tagger.h | 1063 |
| merge_glee | converter | merge GLEE-format checkout trees | WCPLEEANA/eval.h, tagger.h | 485 |
| merge_hist | converter/plotter | per-file histograms -> combined merge.root | WCPLEEANA/master_cov_matrix.h, bayes.h | 1072 |
| merge_hist_xs | converter/plotter | per-file histograms -> combined xs merge (xs mode) | WCPLEEANA/master_cov_matrix.h, bayes.h | 310 |
| merge_pelee_filter | filter | filter PeLEE-format ntuples | WCPLEEANA/eval.h, tagger.h | 524 |
| merge_pelee_nuwro_truth | converter | merge PeLEE NuWro truth branches | WCPLEEANA/eval.h | 505 |
| merge_xf | converter | merge flux-systematic checkout trees | WCPLEEANA/eval.h, tagger.h | 705 |
| numi_filter | filter | NuMI beamline good-run filter | WCPLEEANA/eval.h | 769 |
| plot_hist | plotter | plot data/MC comparison histograms | WCPLEEANA/master_cov_matrix.h | 1733 |
| plot_hist2 | plotter | extended plot_hist with LEE type-2 overlay style | WCPLEEANA/master_cov_matrix.h | 2039 |
| plot_hist_xspaper | plotter | xs-paper specific plots including response matrix | WCPLEEANA/master_cov_matrix.h, eval.h | 941 |
| pot_counting | util | sum POT from data checkout trees | WCPLEEANA/pot.h | 112 |
| pot_counting_mc | util | sum POT from MC checkout trees | WCPLEEANA/pot.h | 52 |
| print_event | util | print event-level info from merge.root | WCPLEEANA/master_cov_matrix.h | 187 |
| prune_checkout_trees | converter | strip most branches from checkout trees | WCPLEEANA/eval.h | 64 |
| prune_weightmar18_trees | converter | prune weight trees (Mar 2018 knob set) | (stdlib + ROOT) | 216 |
| prune_weightsep24_partial_trees | converter | prune weight trees (Sep 2024, partial GENIE knobs) | (stdlib + ROOT) | 247 |
| prune_weightsep24_trees | converter | prune weight trees (Sep 2024 full knob set) | (stdlib + ROOT) | 220 |
| prune_weightsep24_trees_numi | converter | prune weight trees (Sep 2024, NuMI beam) | (stdlib + ROOT) | 235 |
| read_TLee_v20 | fitter | LEE chi2 fit, GoF, FC scan | WCPLEEANA/TLee.h, Configure_Lee.h | 962 |
| sideband_filter | filter | filter events into sideband regions | WCPLEEANA/eval.h, tagger.h | 465 |
| stat_cov_matrix | cov-matrix | data statistical covariance | WCPLEEANA/master_cov_matrix.h | 146 |
| stat_pred_cov_matrix | cov-matrix | MC prediction statistical covariance | WCPLEEANA/master_cov_matrix.h | 176 |
| test_apply_cuts | util | test harness for cut logic | WCPLEEANA/eval.h, tagger.h | 493 |
| test_gp | util | test Gaussian Process regressor | WCPLEEANA/GPKernel.h, GPRegressor.h | 74 |
| test_readout | util | test TTree readout | WCPLEEANA/eval.h | 141 |
| time_dep | util | time-dependence studies on BNB data | WCPLEEANA/eval.h | 220 |
| time_dep_ext | util | time-dependence studies on EXT data | WCPLEEANA/eval.h | 136 |
| wiener_example | util | example Wiener-SVD unfolding | (standalone) | 208 |
| xf_cov_matrix | cov-matrix | flux + cross-section systematic covariance | WCPLEEANA/master_cov_matrix.h | 159 |
| xs_cov_matrix | cov-matrix | cross-section measurement systematic covariance | WCPLEEANA/master_cov_matrix.h | 223 |

---

## 4. Algorithmic Core Deep-Dives

### 4.1 bdt_convert.cxx (3242 lines)

**Purpose.** Transforms a raw WCP output ROOT file (containing five TTrees under `wcpselection/`) into a condensed checkout tree. It evaluates multiple TMVA BDT scores per event, applies good-run quality cuts, and writes filtered copies of all five trees to the output file.

**Input/output files and TTree names.**
- Input: single ROOT file, argv[1]; reads `wcpselection/T_BDTvars`, `T_eval`, `T_pot`, `T_PFeval`, `T_KINEvars`.
- Output: single ROOT file, argv[2]; writes same five trees under `wcpselection/`.

**Core logic summary.** The program first parses optional flags for BDT weight-cut value (`-c`), a training-list file (`-l`) to exclude known run/subruns, a global file-type string (`-g`), a NuMI mode flag (`-n`), and a skip-cuts flag (`-s`). It then constructs an in-memory set from a hardcoded `good_run_list_vec` spanning runs 1 through 5, plus two auxiliary exclusion sets for low-electron-lifetime runs and low-neutrino-count NuMI RHC runs (lines 126-703). A preliminary pass over `T_eval` populates a `remove_set` of (run, subrun) pairs to exclude (e.g. bad subruns from a training list). In the main event loop (line 3090), BDT scores for roughly 20 sub-classifiers are evaluated via `TMVA::Reader` objects (e.g. `cal_br3_3_bdt`, `cal_nue_score`), the composite `nue_score` XGBoost model is applied, combined numu and nue flags are derived, and events failing the good-run or lifetime cuts (lines 3178-3180) are skipped. A second loop over `T_BDTvars` fills the output TTrees. A final loop over `T_pot` writes filtered POT records.

**Red flags.**
- `bdt_convert.cxx:126`: `good_run_list_vec` is a ~565-line static initializer hardcoding run numbers for all five run periods directly in the source file. Any future run must be added here and the binary recompiled.
- `bdt_convert.cxx:693-699`: `low_lifetime_runs` and `low_neutrino_count_numi_run2RHC` are further hardcoded lists. Three separate static-initializer lists performing the same logical function (run-level quality filter) with no shared configuration file.
- `bdt_convert.cxx:109-115`: Input `TFile` `file1` is never explicitly closed; only `file2` (the output) is written/closed.

---

### 4.2 convert_checkout_hist.cxx (366 lines) and convert_checkout_hist_xs.cxx (390 lines)

**Purpose.** Reads a single condensed checkout tree (from `bdt_convert`) and fills analysis histograms for all defined channels and weighting modes, producing a per-input-file ROOT histogram file used by `merge_hist`.

**Input/output files and TTree names.**
- Input: argv[1]; reads `wcpselection/T_BDTvars`, `T_eval`, `T_pot`, `T_PFeval`, `T_KINEvars`.
- Output: argv[2]; named `TH1F` histograms plus a `T` tree with the POT value.

**Core logic summary.** A `CovMatrix` object provides histogram specifications from configuration files. POT is summed from `T_pot`, overridable via `cov.get_ext_pot()`. Three histogram sets are built from `cov.get_histograms()`: CV (mode 0), error-squared (mode 1), and cross-term (mode 2). Branch status flags selectively enable needed branches. For each event, cuts in `WCPLEEANA/cuts.h` select the observation channel and the histogram is filled with the event weight. The `_xs` variant adds `cov.add_xs_config()` and allocates additional `TH1F` signal and `TH2F` response-matrix histograms for cross-section channels.

**Red flags.**
- Both files use `#include "init.txt"` at line 58/51: a non-standard preprocessor include of a code fragment, making the translation unit non-self-contained.
- Neither file closes the input `TFile` before program exit.

---

### 4.3 merge_hist.cxx (1072 lines)

**Purpose.** Aggregates the per-file histograms produced by `convert_checkout_hist` across all input files listed in the `CovMatrix` configuration, normalizing by POT and optionally adding LEE signal, to produce the final combined `merge.root` used as input to `read_TLee_v20`.

**Input/output files and TTree names.**
- Input: per-file histogram ROOT files as registered in `cov.get_map_inputfile_info()`; each contains a `T` TTree with a `pot` branch and named `TH1F` histograms.
- Output: `merge.root` (hardcoded name, line 998); contains `hdata_obsch_N` and `hmc_obsch_N` histograms and a `cov_mat_add` TMatrixD.

**Core logic summary.** `cov.get_map_inputfile_info()` returns a map from filename to a 7-tuple (file type, run period, output filename, external POT, file number, etc.). Each registered file is opened and POT is read from its `T` tree (lines 86-89); data POT per period is recorded separately. Three histogram variants (CV, err2, cross) are loaded into `map_name_histogram`. A second pass identifies data files (filetype 5 or 15) and creates output histogram skeletons. `cov.fill_data_histograms()` and `cov.fill_pred_histograms()` populate them with POT-normalized sums. An optional Bayesian-error branch (flag `flag_err==2`) computes posterior uncertainties bin by bin. Output is written to `merge.root` at line 998.

**Red flags.**
- `merge_hist.cxx:86`: `temp_file = new TFile(out_filename)` overwrites the pointer each iteration with no intervening `Close()`. All input files are left open simultaneously, relying on ROOT's global file list for cleanup at process exit. For large file sets this can exhaust file descriptors.
- `merge_hist.cxx:998`: The output filename `merge.root` is hardcoded as a string literal with no command-line override mechanism.

---

### 4.4 merge_det.cxx (1063 lines)

**Purpose.** Combines a central-value (CV) checkout tree and a detector-systematic variant checkout tree into a single output file, retaining only events present in both samples (matched by run+event number). This paired output is required by `det_cov_matrix` to compute detector systematic uncertainties.

**Input/output files and TTree names.**
- Input: argv[1] = CV checkout ROOT; argv[2] = detector-systematic checkout ROOT; both under `wcpselection/`.
- Output: argv[3]; writes `T_eval_cv`, `T_eval_det`, `T_BDTvars_cv`, `T_BDTvars_det`, `T_KINEvars_cv`, `T_KINEvars_det`, `T_PFeval_cv`, `T_PFeval_det`, `T_pot_cv`, `T_pot_det` under `wcpselection/`.

**Core logic summary.** Both input files are opened and all five TTree pointers retrieved (lines 41-54). A preselection pass builds `map_re_entry_cv` keyed by `(run, event)` (lines 909-927) and `map_re_entry_det` similarly (lines 931-949), admitting only events passing basic quality cuts (match_found, stm flags). The intersection is computed by iterating `map_re_entry_cv` and looking up each key in `map_re_entry_det` (lines 981-986), yielding `map_cv_det_index` of paired entry indices. The main fill loop (lines 987-1010) reads both trees in lockstep and fills all output trees. POT trees are written in a separate pass intersecting run+subrun keys.

**Red flags.**
- `merge_det.cxx:913,935`: Event matching uses `(run, event)` as the key, without subrun. For the case where event numbers are not globally unique within a run (possible in multi-subrun production), this could silently mismatch events across subruns. Compare with `convert_cv_spec.cxx` which switches to a `(run, subrun, event)` 3-tuple for this reason.

---

### 4.5 det_cov_matrix.cxx / xf_cov_matrix.cxx / xs_cov_matrix.cxx / stat_cov_matrix.cxx

**Purpose.** Each of these four short drivers (146-223 lines) constructs one component of the systematic covariance matrix and writes it to a ROOT file. `det_cov_matrix` handles detector systematics using the merged CV+det tree, `xf_cov_matrix` handles flux and cross-section universe reweighting, `xs_cov_matrix` handles the cross-section measurement systematic (adding signal and response-matrix histograms), and `stat_cov_matrix` handles data statistical covariance. `stat_pred_cov_matrix` (176 lines) handles MC prediction statistical covariance.

**Input/output files and TTree names.**
- All four read per-file histogram ROOT files listed in the `CovMatrix` configuration (same files as `merge_hist`).
- Outputs: `det_cov_matrix` writes `vec_mean_N`, `cov_det_mat_N`, `frac_cov_det_mat_N` to a file named by the configuration; `xf_cov_matrix` writes `cov_xf_mat_N`, `frac_cov_xf_mat_N`; `xs_cov_matrix` writes to `./hist_rootfiles/XsFlux/cov_xs.root` (hardcoded, line 147); `stat_cov_matrix` writes to `./hist_rootfiles/run_data_stat.root` (hardcoded, line 45).

**Core logic summary.** All four share the same skeleton: construct a `CovMatrix` with `./configurations/cov_input.txt` plus a variant-specific config file, call the relevant `gen_*_cov_matrix()` method, compute fractional versions by dividing element `(i,j)` by `mean_i * mean_j`, and write all matrices plus per-channel histograms to the output file.

**Red flags.**
- `det_cov_matrix.cxx:133`: diagonal fallback for zero-prediction bins is hardcoded as `1./16.` (i.e., 25% uncertainty squared), with a comment `// 25% uncertainties ...`. This magic number is not configurable.
- `xs_cov_matrix.cxx:147-149`: Output filename is hardcoded as `"./hist_rootfiles/XsFlux/cov_xs.root"`, not derived from the configuration, despite all other cov-matrix apps using `outfile_name` from configuration.
- `stat_cov_matrix.cxx:45`: Output hardcoded as `"./hist_rootfiles/run_data_stat.root"`.
- `xs_cov_matrix.cxx:192`: Writes the xs covariance under the key `"cov_xf_mat_N"` (same name as the flux covariance) rather than `"cov_xs_mat_N"`, which could cause confusion when loading in `read_TLee_v20`.

---

### 4.6 read_TLee_v20.cxx (962 lines)

**Purpose.** The top-level statistical analysis driver. It loads all covariance matrices and predicted/observed spectra via `TLee` and `config_Lee`, performs chi2 minimization over the LEE signal strength, computes goodness-of-fit for multiple channel combinations, and optionally runs a Feldman-Cousins scan.

**Input/output files and TTree names.**
- Input: covariance ROOT files and spectrum files specified via `Configure_Lee.h` (paths set at compile/config time). Also reads `file_collapsed_covariance_matrix.root` for some GoF sub-paths.
- Output: `file_collapsed_covariance_matrix.root` (line 203); contains all `matrix_absolute_*_cov_newworld` TMatrices, `matrix_pred_newworld`, `matrix_data_newworld`, and a `tree` TTree with GoF chi2 values and configuration flags.

**Core logic summary.** Accepts `-p scaleF_POT` and `-f ifile` flags. A `TLee` object is configured from the `config_Lee` namespace; `Set_Spectra_MatrixCov()`, `Set_POT_implement()`, `Set_TransformMatrix()`, and `Set_Collapse()` load all spectra, covariance inputs, and collapse to observation channels. All covariance sub-matrices are written to `file_collapsed_covariance_matrix.root` (lines 229-250). `Exe_Goodness_of_fit()` is called for each combination of target and support channels (nueCC FC/PC, numuCC, CCpi0, NCpi0; lines 262-484). A `TGraph`-based delta-chi2 scan over LEE strength runs lines 651-697. Dead-code blocks in `if(0){}` (lines 821-938) contain FC toy-throwing and Asimov sensitivity studies. `TApplication::Run()` at line 958 optionally holds the display open.

**Red flags.**
- `read_TLee_v20.cxx:28` (comment): References `ROOTSYS=/home/xji/data0/software/root_build` in a usage comment in the Makefile note, indicating the code was originally developed against a personal ROOT installation.
- Extensive `if(0){}` dead code blocks (lines 118-198, 336-445, 519-530, 821-938) containing shape-only covariance decompositions, FC toy generation, and Asimov sensitivity calculations that are unreachable without editing the source. These should be behind a command-line flag.

---

## 5. Fork Families

### convert_cv vs convert_cv_spec (643 vs 642 lines)

The only substantive delta is the event-matching key. `convert_cv` uses a `std::pair<int,int>` of `(run, event)` to populate `map_re_entry_cv` (line 519 of `convert_cv`). `convert_cv_spec` uses a `std::tuple<int,int,int>` of `(run, subrun, event)` (line 519-520 of `convert_cv_spec`), eliminating the ambiguity when event numbers repeat across subruns. There are also minor whitespace differences and two additional branch-status enables in `convert_cv`. A single `bool flag_use_subrun` command-line flag could collapse these into one binary; the key data structure is the only logic difference.

### merge_hist vs merge_hist_xs (1072 vs 310 lines)

`merge_hist_xs` adds `cov.add_xs_config()` at line 64 and allocates additional `map_name_xs_hists` containers for `_signal` (TH1F) and `_R` (TH2F response matrix) histogram variants alongside the standard CV histograms. The standard `merge_hist` has no cross-section measurement channels and therefore no response-matrix handling. The logic divergence is non-trivial: a runtime flag would need to gate the allocation and filling of two additional histogram types per channel. Feasible but would require moderate refactoring.

### convert_checkout_hist vs convert_checkout_hist_xs (366 vs 390 lines)

The `_xs` variant adds `cov.add_xs_config()` (line 31) and allocates `TH1F* htemp1` (signal histogram) and `TH2F* htemp2` (response matrix) for channels identified as xs channels via `cov.is_xs_chname()`. This mirrors the `merge_hist`/`merge_hist_xs` split. Same collapse-to-one-flag argument applies.

### prune_weight* family (4 apps: prune_weightmar18_trees, prune_weightsep24_trees, prune_weightsep24_partial_trees, prune_weightsep24_trees_numi)

The `mar18` vs `sep24` delta is primarily the list of GENIE/Geant4 knob names handled: `mar18` includes `reinteractions_piminus/piplus/proton_Geant4` knobs and fewer GENIE knobs; `sep24` adds `RPA_CCQE_Reduced_UBGenie` and drops the Geant4 reinteraction knobs (diff lines 32, 56-58, 99, 126-128). The `partial` variant adds logic to save entries with fewer than the full 1000 universes by padding the weight vector to 1000 via copies (lines 196-217 of `_partial`). The `numi` variant differs only in the knob list applied to NuMI geometry weights. A `--knob-set {mar18,sep24,numi}` flag plus a `--partial` flag could unify all four.

### plot_hist family (3 apps: plot_hist, plot_hist2, plot_hist_xspaper)

`plot_hist2` adds a `lee_strength_type2` variable (line 45) that allows drawing a LEE spectrum at strength >100 as a visual overlay without including it in chi2 or error calculations. `plot_hist_xspaper` adds XS paper-specific inset plots, a response matrix pad, and opens two additional hardcoded input files (`./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root` and `run3.root`, lines 97-98) to compute efficiency. These diverge substantially enough that a simple flag would not capture the `_xspaper` specialization; the other two could be collapsed with a `--lee-type2` flag.

### bdt_convert + filter_goodruns + numi_filter (shared good_run_list_vec)

All three apps embed the identical `good_run_list_vec` static initializer (or a close variant). `filter_goodruns` begins with the comment "replicated code from bdt_convert but standalone good runs filter app" (line 2). `numi_filter` has its own copy starting at line 65. Any change to the good-run list must be applied in all three files. This is a textbook violation of DRY; the list should live in a shared header or a configuration file read at runtime.

---

## 6. Hardcoded Paths

| File | Line | Path |
|---|---|---|
| `apps/applyNuMIGeomtryWeights.cxx` | 34 | `/data1/xqian/MicroBooNE/processed_checkout_rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root` (default, overridable via `-i`) |
| `apps/applyNuMIGeomtryWeights.cxx` | 35 | `/home/xqian/wire-cell/wcp-uboone-bdt/scripts/NuMI_Geometry_Weights_Histograms.root` (default, overridable via `-w`) |
| `apps/applyNuMIGeomtryWeights.cxx` | 37 | `/data1/xqian/MicroBooNE/processed_checkout_rootfiles/prodgenie_numi_intrinsic_nue_overlay_run1/nucleoninexsec_FluxUnisim.root` (default output, overridable via `-o`) |
| `apps/xs_cov_matrix.cxx` | 147 | `./hist_rootfiles/XsFlux/cov_xs.root` (hardcoded, not overridable) |
| `apps/stat_cov_matrix.cxx` | 45 | `./hist_rootfiles/run_data_stat.root` (hardcoded, not overridable) |
| `apps/merge_hist.cxx` | 998 | `merge.root` (hardcoded output name, no override) |
| `apps/plot_hist_xspaper.cxx` | 97 | `./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root` (hardcoded input) |
| `apps/plot_hist_xspaper.cxx` | 98 | `./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run3.root` (hardcoded input) |
| `apps/plot_hist_xspaper.cxx` | 282 | `merge_xs.root` (hardcoded input) |
| `apps/read_TLee_v20.cxx` | 28 | `/home/xji/data0/software/root_build` (in comment only, not executed) |

---

## 7. Dead Code: apps/old/

| File | First-line purpose |
|---|---|
| `nueCC_convert.cxx` | BDT conversion for nueCC channel; comment states code was modified from ROOT TMVA tutorial. Predecessor to `bdt_convert.cxx`. |
| `numuCC_convert.cxx` | BDT conversion for numuCC channel; same TMVA tutorial origin comment. Predecessor to the numuCC branch of `bdt_convert.cxx`. |

Both files are excluded from the build by the `smplpkg` glob, which does not descend into subdirectories.

---

## 8. scripts/ Summary

The `scripts/` directory contains standalone ROOT macros that are analysis utilities rather than compiled pipeline components. `plot_FC_new.cc` is a compiled C++ ROOT macro (requires explicit `root -l` or standalone compilation) that reads `file_collapsed_covariance_matrix.root` and produces Feldman-Cousins confidence-level plots; it includes `<iostream>` and ROOT headers directly, indicating it is intended for interactive or batch ROOT sessions. `plot_systematics.cc` is a ROOT interpreted macro (`void plot_systematics()`) that reads `file_collapsed_covariance_matrix.root` and renders systematic covariance sub-matrices; it depends on a local `./src/draw.icc` include. `extractNuMIGeometryWeights.C` is a ROOT macro defining `void extractNuMIGeometryWeights(string, string, string)` that ingests a CV checkout file and the binary `NuMI_Geometry_Weights_Histograms.root` (also in `scripts/`) to produce per-event NuMI geometry-correction weight histograms; this function was later wrapped as the compiled `applyNuMIGeomtryWeights` app. The `read_ratio_v01/v02.C` and `print_knob.C` files are additional interactive analysis macros. None of these scripts is part of the automated pipeline; they are post-processing visualization and cross-check utilities intended for interactive use after the main pipeline has run.
