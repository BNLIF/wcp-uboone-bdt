# Selection Layer Examination

**Repository:** `/nfs/data/1/xqian/wcp-uboone-bdt/`  
**Scope:** `inc/WCPLEEANA/` header files  
**Date:** 2026-04-14

---

## Table of Contents

1. [cuts.h — the selection engine](#1-cutsh--the-selection-engine)
   - 1.1 [File overview and function groups](#11-file-overview-and-function-groups)
   - 1.2 [Kinematic extractors (get_reco_*)](#12-kinematic-extractors-get_reco_)
   - 1.3 [Variable dispatcher (get_kine_var)](#13-variable-dispatcher-get_kine_var)
   - 1.4 [Weight calculator (get_weight) — deep dive](#14-weight-calculator-get_weight--deep-dive)
   - 1.5 [Truth variable dispatcher (get_truth_var)](#15-truth-variable-dispatcher-get_truth_var)
   - 1.6 [Channel cut dispatcher (get_cut_pass)](#16-channel-cut-dispatcher-get_cut_pass)
   - 1.7 [Reweight cut dispatcher (get_rw_cut_pass)](#17-reweight-cut-dispatcher-get_rw_cut_pass)
   - 1.8 [Signal predicates (is_*)](#18-signal-predicates-is_)
   - 1.9 [Operator precedence in is_pi0 and is_cc_pi0](#19-operator-precedence-in-is_pi0-and-is_cc_pi0)
   - 1.10 [em_charge_scale asymmetry](#110-em_charge_scale-asymmetry)
   - 1.11 [Copy-paste topology: is_0p / is_1p / is_0pi / is_numuCC_1mu0p](#111-copy-paste-topology-is_0p--is_1p--is_0pi--is_numucc_1mu0p)
   - 1.12 [Legacy MCC8 binning: mcc8_pmuon_costheta_bin](#112-legacy-mcc8-binning-mcc8_pmuon_costheta_bin)
   - 1.13 [Alternative variable index: alt_var_index](#113-alternative-variable-index-alt_var_index)
2. [Configure_Lee.h](#2-configure_leeh)
3. [tagger.h](#3-taggerh)
4. [Event-record headers: eval.h, kine.h, pfeval.h, weights.h, pot.h](#4-event-record-headers)
   - 4.1 [eval.h](#41-evalh)
   - 4.2 [kine.h](#42-kineh)
   - 4.3 [pfeval.h](#43-pfevalh)
   - 4.4 [weights.h](#44-weightsh)
   - 4.5 [pot.h](#45-poth)
5. [bdt.h](#5-bdth)

---

## 1. cuts.h — the selection engine

**File:** `inc/WCPLEEANA/cuts.h`  
**Length:** 4067 lines  
**Header guard:** `#ifndef UBOONE_LEE_CUTS` / `#define UBOONE_LEE_CUTS` (lines 1–2, `#endif` at line 4067)

This is the central analysis logic file. It includes `tagger.h`, `kine.h`, `eval.h`, and `pfeval.h` and provides all cut predicates, kinematic computations, weight evaluation, and signal-selection channel dispatch used across the analysis.

### 1.1 File overview and function groups

| Group | Approximate line range | Purpose |
|-------|----------------------|---------|
| Global constant | 24 | `em_charge_scale = 0.95` — EM energy rescaling for data |
| Forward declarations | 28–105 | All function prototypes within `namespace LEEana` |
| Kinematic extractors | 108–159 | `get_reco_Enu_corr`, `get_reco_showerKE_corr`, `get_reco_Eproton`, `get_reco_Epion` |
| Truth proton predicate | 162–170 | `is_true_0p` |
| Weight calculator | 172–268 | `get_weight` |
| Truth variable dispatcher | 270–275 | `get_truth_var` |
| Reco variable dispatcher | 278–~1100 | `get_kine_var` (large if-else chain over ~50 named variables) |
| XS signal number | ~1100–1633 | `get_xs_signal_no` |
| Channel cut dispatcher | 1635–3633 | `get_cut_pass` (main selection routing) |
| Reweight cut dispatcher | 3636–3655 | `get_rw_cut_pass` |
| Sideband predicates | 3658–3685 | `is_far_sideband`, `is_near_sideband`, `is_LEE_signal` |
| Truth-match predicates | 3691–3706 | `is_truth_nueCC_inside`, `is_truth_numuCC_inside` |
| FC predicate | 3711–3716 | `is_FC` |
| Pi0 predicates | 3719–3762 | `is_cc_pi0`, `is_pi0` |
| NC predicates | 3766–3783 | `is_NCpio_sel`, `is_NCdelta_sel`, `is_NC` |
| numuCC predicates | 3787–3914 | `is_numuCC`, `is_numuCC_tight`, `is_0p`, `is_1p`, `is_0pi`, `is_numuCC_1mu0p`, `is_numuCC_lowEhad`, `is_numuCC_cutbased` |
| nueCC predicates | 3917–3933 | `is_nueCC`, `is_loosenueCC` |
| Preselection | 3935–3955 | `is_generic`, `is_preselection` |
| MCC8 2D binning | 3958–4030 | `mcc8_pmuon_costheta_bin` |
| Alt-var index | 4034–4065 | `alt_var_index` |

---

### 1.2 Kinematic extractors (get_reco_*)

**`get_reco_Enu_corr`** (lines 108–125): Reconstructed neutrino energy with EM charge scale correction. For data (`flag_data=true`) it loops over `kine_energy_particle` and applies `em_charge_scale` to electrons (`pdgcode==11`, `kine_energy_info==2`), then adds `kine_reco_add_energy`. For MC it returns the raw `kine.kine_reco_Enu`. This is the canonical neutrino energy used throughout all channel selections.

**`get_reco_showerKE_corr`** (lines 127–133): Applies `em_charge_scale` to `pfeval.reco_showerKE` for data only.

**`get_reco_Eproton`** (lines 136–147): Sums kinetic energies for protons (PDG 2212) above a 35 MeV threshold.

**`get_reco_Epion`** (lines 149–159): Sums kinetic energies for charged pions (PDG ±211) with no threshold (a commented threshold at 10 keV was removed).

---

### 1.3 Variable dispatcher (get_kine_var)

Lines 278 to approximately 1100. A large if-else chain over a `TString var_name` argument. It maps roughly 50 named variables to their computed or extracted values. Selected entries:

- `"kine_reco_Enu"` → `get_reco_Enu_corr(kine, flag_data)` (line 283)
- `"reco_showerKE"` → `get_reco_showerKE_corr * 1000.` (line 285, note unit conversion to MeV)
- `"kine_pio_energy_max"` / `"kine_pio_energy_min"` (lines 294–303): applies `em_charge_scale` for data
- `"pi0_mass"` (lines 315–327): applies `em_charge_scale` to `kine_pio_mass` for data
- `"nue_score"` (line 331): capped at 15.99
- `"shower_energy"` (lines 338–342): applies `em_charge_scale` for data
- `"electron_energy"` (lines 343–348): applies `em_charge_scale` to `pfeval.reco_showerMomentum[3]` for data, converts GeV to MeV
- `"median_dQdx"` (lines 358–369): hand-computes median from 7 stored `mip_vec_dQ_dx` values
- `"muon_costheta"` (lines 425–430): returns -2 as sentinel when muon momentum not found

The function has no `else` fallthrough error handler at the end — unrecognized `var_name` values silently fall through and return an uninitialized or zero value. This is a latent correctness risk if a mis-spelled variable name is passed.

---

### 1.4 Weight calculator (get_weight) — deep dive

**Lines 172–268.** The function computes an `addtl_weight` via the reweighting loop, then combines it with named weight products.

**The exact wbin section (lines 209–231) is reproduced in full below for reader judgment:**

```cpp
        int wbin;
        bool flag_pass = get_rw_cut_pass(cut_str, eval, pfeval, tagger, kine);
        if (flag_pass){
          if (var>max_var && overflow) addtl_weight = reweight.back();
          else if(var>max_var) addtl_weight = 1;
          else if (var<min_var && underflow) addtl_weight = reweight[0];
          else if (var>min_var){
            if(equal_binning){
              double bin_len = (max_var-min_var)/reweight.size();
              if(underflow && overflow) bin_len = (max_var-min_var)/(reweight.size()-2);
              else if (underflow || overflow) bin_len = (max_var-min_var)/(reweight.size()-1);
              wbin = floor((var-min_var)/bin_len);
            }else{
              std::vector<double> bins = std::get<9>(rw_info_i);
              for(int b=0; b<bins.size()-1; b++){
                if(var<=bins[b+1] && var>bins[b]){
                  wbin = b;
                  break;
                }
              }
            }
            if(underflow) wbin++;
            addtl_weight *= reweight[wbin];
          }
        }
```

**Issue confirmed:** `int wbin;` at line 209 is declared without initialization. It is only assigned inside the two inner branches — the `equal_binning` path and the custom-bins loop. In the custom-bins path (`!equal_binning`), the loop searches for a matching bin and sets `wbin = b` then `break`s. If `var` falls exactly on or between no pair of consecutive bin edges (e.g., if the `bins` vector is malformed, or if `var` exactly equals `min_var` — note the loop uses strict `var>bins[b]`), the loop completes without setting `wbin`. Execution then falls to `if(underflow) wbin++;` and `addtl_weight *= reweight[wbin];`, which reads `wbin` while it is uninitialized, producing undefined behavior. A compiler with optimization enabled may produce any integer value for `wbin`, causing an out-of-bounds vector access.

Additionally, note the outer `else if (var>min_var)` condition: events with `var == min_var` exactly fall through all four branches and receive no reweighting. In the underflow case without that flag, the event also receives no reweighting, which may be correct by design but is worth verifying.

---

### 1.5 Truth variable dispatcher (get_truth_var)

**Lines 270–275:**

```cpp
double LEEana::get_truth_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval,
                              TaggerInfo& tagger, TString var_name){
  if(var_name == "truth_energyInside"){
    return eval.truth_energyInside;
  }else {std::cout<<"Unknown truth var, check configurations"<<std::endl;}
  return 0;
}
```

**Issue confirmed:** The function recognizes exactly one variable name (`"truth_energyInside"`). For every other name it prints a message and silently returns 0. The `get_weight` function calls `get_truth_var` to obtain the variable for reweighting; if a reweighting configuration specifies any truth variable other than `"truth_energyInside"`, the variable value will always be 0 and reweighting bins will always index from the underflow. This is not a compile-time error but is a runtime logic trap.

---

### 1.6 Channel cut dispatcher (get_cut_pass)

**Lines 1635–3633.** This is the largest function in the file. Every call rebuilds a `std::map<std::string, bool> map_cuts_flag` from scratch (line 1652). The map contains approximately 35 truth-level classification entries (numuCCinFV, nueCCinFV, NCinFV, CCMEC, CCQE, CCRES, CCDIS, NCDeltainFV, etc.) all recomputed from the raw event branches on each invocation. Following map construction, the function applies the `add_cut` string logic (a `_`-delimited AND of flag names from the map) and then branches over roughly 100 named channels (`ch_name`) to apply the physical selection criteria.

**Performance implication:** Because `map_cuts_flag` is rebuilt on every call, callers that invoke `get_cut_pass` for multiple channels on the same event pay the cost of the full map construction each time. The map is local to the function and cannot be cached between calls.

**FV definition** is hardcoded inline at line 1649:
```
eval.truth_vtxX > -1 && eval.truth_vtxX <= 254.3
eval.truth_vtxY > -115.0 && eval.truth_vtxY <= 117.0
eval.truth_vtxZ > 0.6 && eval.truth_vtxZ <= 1036.4
```
This duplicates FV coordinates that also appear elsewhere and creates a maintenance risk.

---

### 1.7 Reweight cut dispatcher (get_rw_cut_pass)

**Lines 3636–3655.** Recognizes five named cut strings: `"NCPi0"`, `"NCPi0_Np"`, `"NCPi0_0p"`, `"NCDeltaNp_scale"`, `"NCDelta0p_scale"`. Any other string prints a message and returns `false`. This is called from `get_weight` during the reweighting loop.

---

### 1.8 Signal predicates (is_*)

**`is_nueCC`** (lines 3917–3924): `numu_cc_flag >= 0 && nue_score > 7.0`

**`is_loosenueCC`** (lines 3927–3932): `numu_cc_flag >= 0 && nue_score > 4.0`

**`is_numuCC`** (lines 3787–3793): `numu_cc_flag >= 0 && numu_score > 0.9`

**`is_numuCC_tight`** (lines 3796–3802): adds `pfeval.reco_muonMomentum[3] > 0`

**`is_NC`** (lines 3778–3783): `(!cosmict_flag) && numu_score < 0`

**`is_preselection`** (lines 3943–3955): requires `match_found == 1`, `stm_eventtype != 0`, all STM veto flags zero, `stm_clusterlength > 0`. Handles two era variants: when `is_match_found_int` is set, it reads `match_found_asInt` instead of the `Bool_t match_found` (accommodating file format differences between data-taking eras).

**`is_generic`** (lines 3935–3940): extends preselection with `stm_clusterlength > 15`.

---

### 1.9 Operator precedence in is_pi0 and is_cc_pi0

The two functions differ in how they handle the `kine_pio_flag` condition.

**`is_cc_pi0`** (lines 3719–3739, data branch at line 3731, MC branch at line 3735):

```cpp
if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1* em_charge_scale > 40 && ...)
```

Both data and MC paths require `(kine_pio_flag==1 AND vtx_dis<9)`. The parentheses around the flag-1 condition are explicit; no flag-2 alternative is allowed. This represents a "with vertex only" selection consistent with the comment at line 76.

**`is_pi0`** (lines 3743–3762, data branch at line 3754, MC branch at line 3758):

```cpp
if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && ...)
```

Here the inner parentheses contain `flag==1 && vtx_dis<9 || flag==2`. In C++, `&&` has higher precedence than `||`, so the compiler parses this as:

```
( (flag==1 && vtx_dis<9) || (flag==2) )
```

This is the intended logic: accept either a flag-1 pi0 with vertex attachment, or a flag-2 pi0 without vertex requirement. The parenthesization is therefore correct by C++ operator precedence rules and matches the comment at line 72 ("with and without vertex"). **There is no precedence bug here**, but the pattern looks suspicious at a glance because the outer parentheses suggest the author intended to group the entire flag condition before applying `&&` to the energy cuts — which is exactly what happens.

By contrast, `is_cc_pi0` intentionally drops the `|| flag==2` alternative, making it strictly vertex-attached only.

---

### 1.10 em_charge_scale asymmetry

`em_charge_scale = 0.95` is declared at line 24 as a non-const `double` in `namespace LEEana`. It is applied:

- In `get_reco_Enu_corr` (line 114): to electron-type energy deposits for **data only**
- In `get_reco_showerKE_corr` (line 129): to shower KE for **data only**
- In `get_kine_var` for several variables (lines 296, 300, 321, 340, 345): applies to pi0 energies, pi0 mass, shower energy, electron energy — **data only** in all cases
- In `is_cc_pi0` (lines 3729–3731): `em_charge_scale` applied to `kine_pio_energy_1`, `kine_pio_energy_2`, and `kine_pio_mass` **only in the data branch** (lines 3723–3733). The MC branch (lines 3734–3736) uses raw values.
- In `is_pi0` (lines 3746–3755): same pattern — scale applied in data branch only (lines 3746–3755), not in MC branch (lines 3757–3759).

The asymmetry is consistent and intentional: MC simulation is assumed not to need the charge-scale correction (it is an empirical data/MC energy-scale adjustment). However, the result is that the pi0 mass and energy thresholds in the selection are effectively tighter for data than for MC by a factor of 0.95. This is a deliberate design choice, not a bug, but it means that MC and data pass different effective threshold values.

---

### 1.11 Copy-paste topology: is_0p / is_1p / is_0pi / is_numuCC_1mu0p

All four functions (lines 3805–3887) share an identical body structure:

1. Guard: `if (tagger_info.numu_cc_flag >= 0)`
2. Declare `int Nproton = 0; int Npion = 0;`
3. Loop over `kine_energy_particle` with identical thresholds: proton PDG 2212 at 35 MeV, pion PDG ±211 at 10 MeV
4. Differ only in the final condition:
   - `is_0p` (line 3820): `if(Nproton==0) flag = true`
   - `is_1p` (line 3841): `if(Nproton==1) flag = true`
   - `is_0pi` (line 3863): `if(Npion==0) flag = true`
   - `is_numuCC_1mu0p` (line 3884): adds `tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0` to the outer guard, then `if(Nproton==0) flag = true`

The loop body is duplicated four times with no shared helper. A change to particle thresholds (e.g., the proton 35 MeV cutoff) must be made in four places. The comment in each function says "KE threshold: 50 MeV, 1.5 cm?" but the code uses 35 MeV — the comment is stale.

---

### 1.12 Legacy MCC8 binning: mcc8_pmuon_costheta_bin

**Lines 3958–4030.** A two-dimensional binning function returning an integer bin number (1–42) or -10000 for out-of-range events. It encodes the MCC8 double-differential muon momentum × cos(theta) analysis binning:

| cos(theta) slice | p_muon bins | Bin numbers |
|-----------------|-------------|-------------|
| [-1, -0.5) | 5 bins: [0,0.18), [0.18,0.30), [0.30,0.45), [0.45,0.77), [0.77,2.5) | 1–5 |
| [-0.5, 0) | same 5 bins | 6–10 |
| [0, 0.27) | same 5 bins | 11–15 |
| [0.27, 0.45) | 4 bins (no [0,0.18) split) | 16–19 |
| [0.45, 0.62) | 4 bins | 20–23 |
| [0.62, 0.76) | 4 bins | 24–27 |
| [0.76, 0.86) | 5 bins (adds [1.28,2.5) split) | 28–32 |
| [0.86, 0.94) | 5 bins | 33–37 |
| [0.94, 1.00) | 5 bins | 38–42 |

Nine cos(theta) slices with variable numbers of momentum bins. The forward-going slices (cos(theta) > 0.76) have finer momentum granularity. Exactly 42 bins total. The function uses `and`/`or` (C++ alternative tokens) rather than `&&`/`||`.

---

### 1.13 Alternative variable index: alt_var_index

**Lines 4034–4065.** Lazy-loads a 2D binning configuration from a text file (default `./configurations/alt_var_xbins.txt`) into a static-like `std::map<std::string, TH1F> map_var_hist` (declared at line 104 in the namespace). On first call it parses the file; on subsequent calls the map is already populated. Returns a linearized 2D bin index `bin1 + (bin2-1)*nBins1`, or -1 if either variable is out of range (underflow/overflow excluded).

Note: `map_var_hist` is a global variable in the `LEEana` namespace defined at line 104 as `std::map<std::string, TH1F> map_var_hist;`. This means it is shared across translation units that include `cuts.h` (see also the Configure_Lee.h discussion below about ODR issues).

---

## 2. Configure_Lee.h

**File:** `inc/WCPLEEANA/Configure_Lee.h`  
**Length:** 137 lines

**Include guard:** None. The file has no `#ifndef` / `#define` guard and no `#pragma once`. Every translation unit that includes it gets a fresh copy of all its variable definitions.

**Namespace:** `namespace config_Lee { ... }` wrapping lines 1–137.

**Variable declarations:** All are non-const mutable variables defined (not merely declared) at namespace scope:

```cpp
TString spectra_file = "./new_TLee_input_opendata5e19/merge.root";  // line 63
TString flux_Xs_directory = "./new_TLee_input_opendata5e19/flux_Xs/";  // line 64
TString detector_directory = "./new_TLee_input_opendata5e19/det/";  // line 65
TString mc_directory = "./new_TLee_input_opendata5e19/mc_stat/";  // line 66
int channels_observation = 7;   // line 70
int syst_cov_flux_Xs_begin = 1; // line 73
int syst_cov_flux_Xs_end   = 19;// line 74
int syst_cov_mc_stat_begin = 0; // line 76
int syst_cov_mc_stat_end   = 99;// line 77
bool flag_display_graphics = 1; // line 89
bool flag_syst_flux_Xs    = 1;  // line 93
// ... many more bool flags ...
double Lee_strength_for_GoF = 0; // line 107
```

**Implications:**

1. **One Definition Rule (ODR) violation risk:** If more than one `.cxx` file includes `Configure_Lee.h`, each translation unit defines its own copy of these variables. This causes linker errors for non-inline, non-template variables under the ODR unless only one compilation unit ever includes this header (which would be an unusual pattern). In a header-only codebase driven from a single `main()` TU, this may never trigger, but it is fragile.

2. **Mutability:** Every flag is freely mutable. Nothing is `const` or `constexpr`. A function receiving `config_Lee::flag_syst_detector` by value would silently use a stale value if the flag were changed elsewhere.

3. **Large commented-out blocks:** Lines 5–61 contain approximately a dozen commented-out alternative path sets from previous analysis iterations (various fake datasets, NuMI configurations, etc.), creating maintenance noise. The currently active configuration (lines 63–66) corresponds to the `new_TLee_input_opendata5e19` dataset.

---

## 3. tagger.h

**File:** `inc/WCPLEEANA/tagger.h`  
**Length:** 2994 lines  
**Header guard:** `#ifndef UBOONE_LEE_TAGGER` / `#define UBOONE_LEE_TAGGER`

### Structure

The file declares `struct TaggerInfo` (lines 6–772) and three free functions in `namespace LEEana`:
- `clear_tagger_info` (lines 780–~1100): zeroes all scalar fields and calls `->clear()` on all vector-pointer fields
- `set_tree_address` (lines ~1100–~2200): maps ROOT TBranch addresses for reading
- `put_tree_address` (lines ~2200–2991): creates output TBranch bindings for writing

The three blocks are parallel in structure — they cover the same fields in the same order.

### Field count

`TaggerInfo` contains approximately **500 named data members**: roughly 350 scalar `float` fields and approximately 75 `std::vector<float>*` pointer fields. The scalar fields represent individual tagger sub-classifier scores and input variables (cosmic tagger, gap tagger, MIP ID, pi0 ID, shower-to-wall, bad-reconstruction sub-classifiers, track-overclustering, etc.). The vector fields hold per-object arrays used by taggers that evaluate multiple candidates (e.g., `br3_3_v_energy`, `cosmict_10_length`, `numu_cc_1_particle_type`).

### Float-holding-integer pattern

Several logically boolean or integer fields are declared as `float` instead of the appropriate type. Examples confirmed at lines 739–763:

```cpp
float match_isFC;       // boolean stored as float
float truth_isCC;       // boolean stored as float
float truth_vtxInside;  // boolean stored as float
float truth_nuPdg;      // PDG code (integer) stored as float
float kine_pio_flag;    // integer flag stored as float
float event_type;       // integer category stored as float
```

These fields are used for ROOT TTree training output (BDT input preparation). The ROOT branch strings confirm this — e.g., `T_tagger->Branch("truth_isCC",&tagger_info.truth_isCC,"data/F")` at line 3961. Storing PDG codes and categorical integers as float introduces rounding risk for large PDG codes, though within the values used in this analysis (PDG codes ≤ 2212, flag values 0/1) there is no actual precision loss.

### Era-version toggle mechanism

`clear_tagger_info` at lines 780–782 reads two boolean flags:
```cpp
tagger_info.flag_nc_gamma_bdt = false;
tagger_info.flag_nc_gamma_0track_bdt = false;
```
The `set_tree_address` and `put_tree_address` functions check these flags before binding `nc_delta_score`, `nc_pio_score`, `nc_delta_0track_score`, and `nc_delta_ntrack_score` branches. This allows processing files from eras before these BDTs were added without branch-not-found errors.

---

## 4. Event-record headers

### 4.1 eval.h

**File:** `inc/WCPLEEANA/eval.h` (310 lines)  
**Header guard:** `#ifndef UBOONE_LEE_EVAL`

**Purpose:** Defines `struct EvalInfo`, the primary event-level record holding flash matching, STM (Space-Time Manager) veto results, truth-level neutrino properties, and event weights.

**Key fields:** `run/subrun/event`, `match_found` / `match_found_asInt` (era toggle), `match_isFC`, `stm_eventtype`, `stm_lowenergy`, `stm_LM`, `stm_TGM`, `stm_STM`, `stm_FullDead`, `stm_clusterlength`, `truth_nuEnergy`, `truth_energyInside`, `truth_nuPdg`, `truth_isCC`, `truth_vtxInside`, `truth_vtxX/Y/Z`, `weight_spline`, `weight_cv`, `weight_lee`, plus era-gated PeLEE (`flag_pl`) and gLEE (`flag_gl`) extended fields.

**Pointer initialization:** No raw pointer members. No `init_pointers` / `del_pointers` needed.

**Era-version toggle:** `is_match_found_int` (set at line 139 when branch `flash_found_asInt` exists) and `is_file_type` (line 143) control whether legacy Bool_t branches or newer Int_t branches are read. `flag_pl` and `flag_gl` gate optional PeLEE and gLEE branches. The `set_tree_address` function (line 97) probes branch existence with `tree0->GetBranch(...)` before binding.

**MC-only flag:** The `flag` parameter (default 1) on `set_tree_address` and `put_tree_address` gates binding of all truth and weight branches, allowing data files to be processed without those branches.

### 4.2 kine.h

**File:** `inc/WCPLEEANA/kine.h` (139 lines)  
**Header guard:** `#ifndef UBOONE_LEE_KINE`

**Purpose:** Defines `struct KineInfo` holding reconstructed energy quantities: neutrino energy, particle energy list, pi0 kinematics, and optional deep-learning energy estimates.

**Key pointer fields:** `kine_energy_particle`, `kine_energy_info`, `kine_particle_type`, `kine_energy_included` — all `std::vector<float/int>*`.

**Pointer initialization:** There is NO `init_pointers` or `del_pointers` function. The pointers are declared but never initialized in the struct definition. They receive valid addresses only when `set_tree_address` calls `tree0->SetBranchAddress(...)` which sets them to point at ROOT-managed memory. If `clear_kine_info` is called before `set_tree_address` has been called, the pointer dereferences at lines 46–49 are undefined behavior:

```cpp
void LEEana::clear_kine_info(KineInfo& tagger_info){
  tagger_info.kine_reco_Enu=0;
  tagger_info.kine_reco_add_energy=0;
  tagger_info.kine_energy_particle->clear();   // line 46 — no null check
  tagger_info.kine_energy_info->clear();        // line 47 — no null check
  tagger_info.kine_particle_type->clear();      // line 48 — no null check
  tagger_info.kine_energy_included->clear();    // line 49 — no null check
```

**Confirmed:** there is no null-check before any of these four dereferences. If the struct is default-constructed (all pointer members are indeterminate), calling `clear_kine_info` crashes.

**Era toggle:** `has_dl_ee` (set to `false` at line 92, set to `true` at line 102 if branch `vlne_v4_numu_full_primaryE` exists) gates the deep-learning energy branch bindings.

### 4.3 pfeval.h

**File:** `inc/WCPLEEANA/pfeval.h` (604 lines)  
**Header guard:** `#ifndef UBOONE_LEE_PFEVAL`

**Purpose:** The most complete event record, holding full PFParticle-level reconstruction and truth-matching information: vertex coordinates, shower/muon/proton momenta, MCFlux kinematics, per-track truth arrays (up to 10000 tracks), and all truth interaction type fields.

**Key pointer fields:** `truth_process` (`vector<string>*`), `truth_daughters` and `reco_daughters` (`vector<vector<Int_t>>*`), `fMC_trackPosition` (`TObjArray*`), `reco_process` (`vector<string>*`).

**Pointer initialization:** `init_pointers` (lines 146–154) allocates all five pointers with `new`. `del_pointers` (lines 156–162) deletes them. `clear_pfeval_info` (lines 165–~270) checks `flag_init_pointers` before calling `init_pointers`, preventing double-initialization.

**Era toggle:** `flag_NCDelta`, `flag_showerMomentum`, `flag_recoprotonMomentum`, `flag_pf_truth`, `flag_pf_reco` gate branch binding in `set_tree_address`.

### 4.4 weights.h

**File:** `inc/WCPLEEANA/weights.h` (247 lines)  
**Header guard:** `#ifndef UBOONE_LEE_WEIGHTS`

**Purpose:** Defines `struct WeightInfo` for storing GENIE universe weights and flux systematics. Holds `weight_cv`, `weight_spline`, `weight_lee` scalars and approximately 25 `std::vector<float>*` members for named systematic universes (flux unisim vectors, GENIE knob vectors).

**Era toggle:** `flag_sep_28` (line 8) controls whether September-28 additions (`reinteractions_piminus_Geant4`, etc.) replace the older `RPA_CCQE_Reduced_UBGenie` field.

**Pointer initialization:** No explicit `init_pointers` function. Like `KineInfo`, the vectors receive addresses through ROOT branch binding.

### 4.5 pot.h

**File:** `inc/WCPLEEANA/pot.h` (39 lines)  
**Header guard:** `#ifndef UBOONE_LEE_POT`

**Purpose:** Defines `struct POTInfo` with six fields: `runNo`, `subRunNo`, `pot_tor875`, `pot_tor875good`, `spill_tor875`, `spill_tor875good`. Provides `set_tree_address` and `put_tree_address` for the sub-run POT counting tree.

No pointer members, no era toggle, no `init_pointers`. Simplest header in the set.

---

## 5. bdt.h

**File:** `inc/WCPLEEANA/bdt.h` (650 lines)  
**Header guard:** None. The file opens with `namespace LEEana{` directly (line 1) and has no `#ifndef` guard.

### Declarations

The file begins (lines 1–127) with forward declarations of approximately 20 `float cal_*_bdt(...)` and `float cal_*_bdts_xgboost(...)` functions inside `namespace LEEana`. These take a `TaggerInfo&` and a `TMVA::Reader&` as their first two arguments, plus additional `float&` arguments for the specific input variables needed by each BDT.

### Implementations

**The implementations ARE present in this file.** Despite the context hint suggesting they might be absent, bdt.h contains the full function bodies starting at line 131. All `cal_*` functions are implemented inline within the header:

- `cal_nc_delta_bdts_xgboost` (lines 131–138): calls `reader.EvaluateMVA("MyBDT")` and applies a logit transform `log10((1+val1)/(1-val1))`
- `cal_nc_delta_0track_bdts_xgboost`, `cal_nc_delta_ntrack_bdts_xgboost`, `cal_nc_pio_bdts_xgboost`, `cal_numu_bdts_xgboost` (lines 140–179): all apply the same logit transform pattern
- The per-sub-tagger functions (`cal_cosmict_10_bdt`, `cal_numu_1_bdt`, `cal_numu_2_bdt`, `cal_br3_3_bdt`, etc., lines 182–~650): iterate over the per-object vector fields in `TaggerInfo`, copy values into the `float&` reference arguments (which must be bound to TMVA reader input variables), call `reader.EvaluateMVA("MyBDT")`, and accumulate the minimum or maximum score over all candidates

No separate `.cxx` file defines these functions. The implementation is header-only. The only BDT-related `.cxx` in the repo is `apps/bdt_convert.cxx`, which is an application that uses these functions, not a definition source.

**No include guard** means that any translation unit including `bdt.h` more than once (or any project that includes it from two headers both pulled into the same TU) will get duplicate definitions of all these functions, causing linker errors. In practice the framework appears to be structured so that `bdt.h` is included through a single chain.

---

## Summary of Key Findings

| Finding | Location | Severity |
|---------|----------|----------|
| `wbin` used uninitialized in custom-bin reweighting path when no bin matches | cuts.h:209–231 | High — undefined behavior / potential out-of-bounds vector access |
| `get_truth_var` returns 0 silently for all var names except `"truth_energyInside"` | cuts.h:270–275 | High — silent logic error for any other reweighting variable |
| `get_cut_pass` rebuilds full `map_cuts_flag` on every call | cuts.h:1652 | Medium — performance cost proportional to call frequency |
| `clear_kine_info` dereferences four pointers without null check | kine.h:46–49 | High — crash if called before `set_tree_address` |
| `Configure_Lee.h` has no include guard, all globals are mutable, non-const | Configure_Lee.h:1–137 | Medium — ODR violation risk; globals silently mutable |
| `bdt.h` has no include guard | bdt.h:1 | Medium — multiple-inclusion linker error risk |
| `is_0p`/`is_1p`/`is_0pi`/`is_numuCC_1mu0p` body duplicated four times; stale threshold comment | cuts.h:3805–3887 | Low — maintenance hazard |
| `get_kine_var` has no fallthrough error handler for unknown variable names | cuts.h:278–~1100 | Low — silent incorrect values for mis-spelled variables |
| `em_charge_scale` is mutable non-const at namespace scope | cuts.h:24 | Low — could be inadvertently modified |
