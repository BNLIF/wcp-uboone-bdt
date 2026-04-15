# wcp-uboone-bdt Code Examination

**What this is:** A read-only static audit of the
[wcp-uboone-bdt](https://github.com/BNLIF/wcp-uboone-bdt) repository,
the post-reconstruction analysis framework for the MicroBooNE Low-Energy
Excess search. The examination covers architecture, bug identification,
algorithm documentation, efficiency, and open design questions.
No source files were modified. Completed 2026-04-14.

---

## Reading Order

| # | Document | Contents |
|---|----------|----------|
| 1 | [01_architecture.md](01_architecture.md) | Start here. Five-stage pipeline overview, component dependency diagram, configuration files, hardcoded paths, build system, source tree layout. |
| 2 | [07_algorithms.md](07_algorithms.md) | Read second to understand what the code is computing: CNP chi-squared, sideband constraint (Schwartz decomposition), Feldman-Cousins scan, LEE signal-strength fit, systematic covariance pipeline, Wiener-SVD unfolding, Bayesian MC-stat convolution, Gaussian-process regression. Every claim is tied to a file and line number. |
| 3 | [02_core_framework.md](02_core_framework.md) | Technical examination of the core src/ classes: Util, WienerSVD, GPKernel/GPRegressor/GPSmoothing, Bayes, CovMatrix (master_cov_matrix), TLee. |
| 4 | [03_selection_layer.md](03_selection_layer.md) | The selection and weighting layer in inc/WCPLEEANA/: cuts.h (get_weight, get_cut_pass, signal predicates, em_charge_scale), Configure_Lee.h, tagger.h, eval/kine/pfeval/weights/pot headers, bdt.h. |
| 5 | [04_apps_pipeline.md](04_apps_pipeline.md) | All ~35 compiled executables: complete triage table, deep-dives on bdt_convert, merge_hist, merge_det, convert_checkout_hist, det/xf/xs/stat cov matrix drivers, read_TLee_v20. Fork families and hardcoded paths. |
| 6 | [05_bugs.md](05_bugs.md) | Prioritized bug catalogue. Read the High items first. Each entry has a symptom, explanation, fix sketch, and cross-reference. |
| 7 | [06_efficiency.md](06_efficiency.md) | Performance and resource-use findings: per-event map construction, string-comparison dispatch, FFT parameter sizing, TFile leaks, double-pass file iteration, branch-address overhead. |
| 8 | [08_open_questions.md](08_open_questions.md) | Items that cannot be resolved by code reading alone and require author input. |

---

## Bug Counts (05_bugs.md)

| Severity | Count | Examples |
|----------|-------|---------|
| High | 4 | B-01 correlated RNG seeds, B-02 GP early-exit on log-transformed coords, B-03 uninitialised wbin index, B-04 dead axis-config blocks |
| Medium | 9 | B-05 Configure_Lee.h ODR risk, B-06 Bayes TF1 leak/double-delete, B-07 stale GP derivative, B-08 GP parameter leak, B-09 unchecked SVD decomposition, B-10 TFile handle accumulation, B-11 null KineInfo pointer deref, B-12 silent bisection non-convergence, B-13 eof() double-read |
| Low | 8 | B-14 hardcoded POT fallback, B-15 channel-boundary literals, B-16 hardcoded dim_edges, B-17 triplicated good-run list, B-18 float equality on flag, B-19 wrong ROOT key for xs cov, B-20 2-tuple event matching, B-21 magic 1/16 fallback |
| Style | 5 | B-22 MatrixMatirx typo, B-23 Feldman spelling, B-24 empty Set_TransformMatrix, B-25 hardcoded test_wiener=false, B-26 bdt.h missing include guard |
| **Total** | **26** | |

---

## Scope Limitations

The following areas were **not** examined in depth:

- **scripts/*.C** — standalone ROOT macros (plot_FC_new.cc, plot_systematics.cc, extractNuMIGeometryWeights.C, read_ratio_v01/v02.C, print_knob.C) were surveyed for pipeline role only; no line-by-line review.
- **git history** — no examination of commit history, authorship, or prior versions of any file.
- **Runtime behavior** — no code was executed. All findings are static, based on source reading. I/O performance numbers and memory estimates are analytical.
- **apps/old/** — nueCC_convert.cxx and numuCC_convert.cxx are excluded from the build and were not read.
- **tagger.h line-by-line** — the 2994-line BDT tagger variable struct was examined for structure and branch-address setup only; individual variable definitions were not audited.
- **TMVA / XGBoost weight files** — external BDT model files consumed by bdt_convert are not part of this repository and were not examined.
- **configurations/ files** — these do not exist in the repository; their expected schema was inferred from the parsing code in master_cov_matrix.cxx.

---

## How to Use 05_bugs.md

**Severity rubric:**

- **High** — the defect produces silently wrong physics outputs. These should be fixed before any result is used for a publication or unblinding.
- **Medium** — latent or conditional defects that corrupt results in identifiable circumstances (specific input patterns, optimizer call sequences, or large file counts). Fix before production runs on new datasets.
- **Low** — brittle assumptions or hardcoded constants that could silently break on dataset changes (new run periods, changed binning, different POT). Fix as part of any analysis extension.
- **Style** — naming or dead-code issues with no runtime consequence. Fix opportunistically.

Each entry in 05_bugs.md follows the pattern: **Symptom** (what goes wrong at runtime), **Explanation** (the code path that causes it with file:line references), **Fix sketch** (approach only, no patch), and **Cross-reference** (the examination document that provides additional context for the surrounding code). The bug IDs (B-01 through B-26) are used as cross-references throughout the other documents.

For the four High items, the fix sketches are self-contained and can be implemented independently. B-01 (RNG seeding) and B-04 (dead axis blocks) are the most straightforward. B-02 (GP early-exit) and B-03 (uninitialised wbin) each require understanding the surrounding algorithm context in 07_algorithms.md before applying the fix.
