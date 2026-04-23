# Shuyi's TLee Refactor — File Comparison

**Branch:** `xin_claude_improvement`  
**Source:** `Xiangpan_files/` (five files placed there by Xin from Shuyi/LSY's validated version)  
**Adopted on:** 2026-04-23

Shuyi validated the `wcp-uboone-wcp` and `LEEana` packages and found no issues.
Her version of the TLee code is more convenient than the repository version because
the procedure for reading the input spectra is driven by the ROOT file contents
rather than by hardcoded channel maps in `TLee.cxx`.

---

## Files changed

| File | Location | Lines before | Lines after | Net |
|---|---|---|---|---|
| `TLee.h`           | `inc/WCPLEEANA/` | 206 | 205 | −1 |
| `Configure_Lee.h`  | `inc/WCPLEEANA/` | 139 | 76  | −63 |
| `TLee.cxx`         | `src/`           | 2848 | 2558 | −290 |
| `draw.icc`         | `src/`           | 101 | 101 | 0 (byte-identical) |
| `read_TLee_v20.cxx`| `apps/`          | 962 | 729 | −233 |

---

## Per-file summary

### `TLee.h` — minor cleanup

**Removed members** (reweight systematics pathway, 6 members):
- `bool flag_syst_reweight;`
- `bool flag_syst_reweight_cor;`
- `TMatrixD matrix_input_cov_reweight;`
- `TMatrixD matrix_input_cov_reweight_cor;`
- `TMatrixD matrix_absolute_reweight_cov_newworld;`
- `TMatrixD matrix_absolute_reweight_cor_cov_newworld;`

**Changed in constructor:** four existing flags that were previously uninitialized
are now explicitly set to `false` at construction:
`flag_syst_flux_Xs`, `flag_syst_detector`, `flag_syst_additional`,
`flag_syst_mc_stat`.

**`flag_Lee_minimization_after_constraint`** member and ctor initialization
removed (the branch it guarded is gone from `TLee.cxx`).

No methods were added, removed, or renamed. The public API is unchanged aside
from the six removed data members.

---

### `Configure_Lee.h` — removed flags, new LEE-channel array

**Removed variables:**
- `bool flag_syst_reweight`
- `bool flag_syst_reweight_cor`
- `bool flag_chi2_data_H0`
- `bool flag_dchi2_H0toH1`
- `bool flag_Lee_minimization_after_constraint`
- ~50 commented-out alternative path blocks (all the `fakeset*`,
  `NumiReinteraction*`, `opendata*`, `numi_summary*`, etc.)

**New variables:**
- `int array_LEE_ch[4] = {8,9,0,0};` — declarative list of LEE channel
  indices (non-zero entries are applied). Replaces the hardcoded
  `map_Lee_ch[8]=1; map_Lee_ch[9]=1;` block that used to live inside
  `TLee::Set_Spectra_MatrixCov()`. Alternative size/value presets are
  available as comments in the header.
- `bool flag_Lee_scan_data = 0;` — enables the Feldman-Cousins data scan
  path added in `read_TLee_v20.cxx`.
- `bool flag_GOF = 0;` — goodness-of-fit toggle (consumed inside `TLee`).

**Default value changes:**

| Variable | Old value | New value |
|---|---|---|
| `syst_cov_flux_Xs_end`        | 19  | 17  |
| `syst_cov_mc_stat_end`        | 99  | 98  |
| `flag_display_graphics`       | 1   | 0   |
| `flag_GoF_output2file_default_0` | 1 | 0  |

**Input paths** — Shuyi's file referenced `/home/lsy/input_files/spectrum/run1-3/...`.
At adoption these were restored to the repo's existing relative paths:

```
spectra_file        = "./new_TLee_input_opendata5e19/merge.root"
flux_Xs_directory   = "./new_TLee_input_opendata5e19/flux_Xs/"
detector_directory  = "./new_TLee_input_opendata5e19/det/"
mc_directory        = "./new_TLee_input_opendata5e19/mc_stat/"
```

**`#pragma once`** — added at the top of the adopted file to preserve the
include-guard fix from commit `212e7f4` (B-05).

---

### `TLee.cxx` — new input-reading procedure (main improvement)

#### Input-reading procedure (before)

`Set_Spectra_MatrixCov()` populated `map_input_spectrum_ch_str` with
hardcoded assignments:
```cpp
map_input_spectrum_ch_str[1] = "nueCC_FC_norm";
map_input_spectrum_ch_str[2] = "nueCC_PC_norm";
// ... etc for all channels
map_Lee_ch[8] = 1;
map_Lee_ch[9] = 1;
```
Switching channel sets or adding LEE channels required editing `TLee.cxx`
and recompiling. Multiple alternative blocks for different schemes
(1u0p/1uNp, separate-nue, fake-data variants) were commented out inline.

#### Input-reading procedure (after)

Channels are auto-discovered from the spectra ROOT file:
```cpp
for(int ich=1; ich<=1000; ich++) {
  TH1F *h = (TH1F*)file_spectra->Get(TString::Format("histo_%d", ich));
  if( h == NULL ) break;
  map_input_spectrum_ch_str[ich] = h->GetTitle();  // name from histogram title
  delete h;
}
```
Same pattern for observation channels (`hdata_obsch_%d`). The number of
channels and their names are fully determined by the ROOT file — no recompile
needed when adding or renaming channels.

LEE-channel tagging is now pushed out to `read_TLee_v20.cxx` (see below),
using `config_Lee::array_LEE_ch[]`. `Set_Spectra_MatrixCov()` no longer
seeds `map_Lee_ch`.

#### Other changes in `TLee.cxx`

- **Reweight pathway fully removed** from `Set_Collapse`, `Set_POT_implement`,
  and `Plotting_systematics` (all `matrix_input_cov_reweight*`,
  `matrix_absolute_reweight*` accumulation and plotting code).
- **`Minimization_Lee_strength_FullCov`**: drops the
  `flag_Lee_minimization_after_constraint` branch (~140 lines of alternative
  chi2 with conditional YY|XX constraint). Now tracks `minimization_NDF`
  (set to `matrix_delta.GetNcols()`, decremented by 1 when Lee strength is
  floating).
- **Flux/Xs loop** in `Set_Spectra_MatrixCov`: collapses from 4 branches
  (flux ≤16 / Xs ==17 / reweight ==18 / reweight_cor ==19) to 2
  (flux ≤16 / else → Xs). Covariance files with `idx > 16` now all accumulate
  into `matrix_flux_Xs_frac`.
- Scan-progress printout in `Exe_Feldman_Cousins_Asimov` uses `num_scan-1`
  denominator (minor off-by-one fix).

---

### `draw.icc` — no change

The file is byte-identical to the version already in `src/`. Replaced for
completeness per the instruction.

---

### `read_TLee_v20.cxx` — CLI unchanged, major cleanup

**CLI:** unchanged (`-p <scaleF_POT>` and `-f <ifile>`, same as before).

**New LEE-channel bootstrap** (directly after `Set_config_file_directory`):
```cpp
int size_array_LEE_ch = sizeof(config_Lee::array_LEE_ch)
                       /sizeof(config_Lee::array_LEE_ch[0]);
for(int idx=0; idx<size_array_LEE_ch; idx++) {
  if( config_Lee::array_LEE_ch[idx]!=0 )
    Lee_test->map_Lee_ch[config_Lee::array_LEE_ch[idx]] = 1;
}
```
LEE channels are now declared in one place (`config_Lee::array_LEE_ch`) and
consumed by a generic loop — no more scattered per-channel assignments.

**Removed from `main`:**
- Assignments of `flag_syst_reweight`, `flag_syst_reweight_cor`,
  `flag_Lee_minimization_after_constraint` to `Lee_test->…`.
- Two large `if(0)` shape-only-covariance blocks (~80 lines).
- All `flag_both_numuCC` / `flag_CCpi0_*` / `flag_nueCC_*` per-channel GoF
  dispatch blocks.
- `flag_publicnote` canned GoF block.
- `flag_chi2_data_H0` block (chi2 of data under null hypothesis).
- `flag_dchi2_H0toH1` block (H0-vs-H1 likelihood-ratio test).
- TCanvas/TLine/TLatex code drawing `canv_gh_scan.png`.

**Added:**
- `if( config_Lee::flag_Lee_scan_data )` branch that calls
  `Set_measured_data()` + `Exe_Feldman_Cousins_Data()` (data FC scan).
- New `if(0)`-guarded `Exe_Goodness_of_fit_detailed` blocks (GoF on
  channel-range slices).
- End-of-run diagnostic printout: `map_Lee_ch` contents, MC-stat log range,
  `scaleF_POT` / `ifile`.

**Output gating:** all `tree_config` / `matrix_*.Write(...)` calls are now
inside `if( config_Lee::flag_GoF_output2file_default_0 ){ … }`. Previously
`file_collapsed_covariance_matrix.root` was always written. With the new
default `flag_GoF_output2file_default_0 = 0`, covariance output is opt-in.

`tree_config` branches for `flag_syst_reweight` and `flag_syst_reweight_cor`
are removed. Remaining branches: `flag_syst_flux_Xs`, `flag_syst_detector`,
`flag_syst_additional`, `flag_syst_mc_stat`, `user_Lee_strength_*`,
`user_scaleF_POT`, `vc_val_GOF`, `vc_val_GOF_NDF`.

---

## Scope of impact

Only two translation units reference the five changed files:
- `src/TLee.cxx` (includes `TLee.h` and `draw.icc`)
- `apps/read_TLee_v20.cxx` (includes `TLee.h` and `Configure_Lee.h`)

No other file in `src/` or `apps/` is affected.
