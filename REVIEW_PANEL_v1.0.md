# BandQPI-3D v1.0 — independent review before release

Three reviewers worked independently on the release candidate, each with a different mandate, without
access to the author's test scripts or conclusions. Each reported concrete findings with evidence
(commands, numbers, file locations). Every finding below is listed with what was done about it.

| reviewer | mandate | method |
|---|---|---|
| **R1 — computational physicist** | Are the numbers right? | Wrote an own HSX reader and Bloch-sum diagonaliser (Python/SciPy), validated it against SIESTA's `.EIG`, then recomputed MoS₂, HfCrBr₂ and Si: path eigenvalues, gaps, band edges, electron counts, moment, DOS integrals. |
| **R2 — new user / release engineer** | Does it install and run for someone else? | Unzipped the package on a clean folder, followed README literally, ran the test suite and two real systems, tried every option, broken inputs, spaces in paths, PATH/symlink installs. |
| **R3 — figure and page referee** | Are the figures and the page correct and honest? | Viewed every figure of six systems, drove every tab / selector / checkbox of the interactive pages in a browser (1536 × 754), decoded the QPI arrays and compared them pixel by pixel. |

---

## R1 — physics (verdict: no numerical error found)

| quantity | BandQPI-3D | R1, own code |
|---|---|---|
| reader check: own eigenvalues vs SIESTA `.EIG` | — | ≤ 1.1×10⁻⁴ eV (MoS₂), 1.2×10⁻⁵ (HfCrBr₂), 6.6×10⁻⁵ (Si) |
| path eigenvalues, 30 points, all bands | — | max \|ΔE\| 7.9×10⁻⁵ / 8.9×10⁻⁶ / 4.2×10⁻⁵ eV |
| MoS₂ gap (indirect Γ→K) | 1.594963 eV | 1.594963 eV |
| HfCrBr₂ spin-↓ direct gap | 3.218712 eV | 3.218712 eV |
| Si gap | 0.474994 eV | 0.474994 eV |
| HfCrBr₂ moment (own E_F) | 4.00001 μB | 3.99979 μB |
| HfCrBr₂ own E_F − SIESTA E_F | +0.006270 eV | +0.006278 eV |
| windowed DOS integrals | — | agree to ≤ 0.6 % (MoS₂), ≤ 0.07 % (Si, after grid extrapolation) |

Findings and actions:

| # | finding | action |
|---|---|---|
| 1 | Only the *direct* gap per spin was written; HfCrBr₂ spin-↓ fundamental gap is 3.218649 eV (VBM and CBM at inequivalent M points). | **Fixed**: `spinN_gap_eV` (fundamental gap per channel, 0 for a metallic channel) added to `summary.txt`. |
| 2 | VBM/CBM lines are printed for metallic channels, where they mean nothing. | **Documented** (README §8: valid only when `spinN_type` is gapped); the figures already use them only for gapped channels. |
| 3 | Occupied-band count rounded at SIESTA's E_F. | Kept; for metals it only labels which band is "VB"; no reported number depends on it. Noted. |
| 4 | Periodicity detected from the supercell / k-grid. | Kept on purpose: a direction with supercell 1 has no Hamiltonian coupling across the cell, so the bands do not disperse along it and treating it as non-periodic is exact. Tested with model J (vacuum as the first lattice vector, layer in the y–z plane): identical results. |
| 6 | Mid-gap E_F for insulators is a convention. | **Documented** in README §8. |
| — | MoS₂ input path has K written as (0.33, 0.33). | **Documented** (README §8); gaps come from the whole zone, not from the path (also tested by model I). |
| 10 | Spin–orbit branch not covered by R1's data. | Covered elsewhere: CrCl₃ (spin–orbit, nspin 8) agrees with SIESTA's `.bands` to 5×10⁻⁵ eV for every band (check B03). |

## R2 — installation and use

| # | severity | finding | action |
|---|---|---|---|
| 1 | blocker | `install.sh` crashed at the Python check (`importlib.util` not imported). | **Fixed**; fresh-folder `bash install.sh --test` now runs end to end. |
| 2 | blocker | No licence file; plotly.js MIT notice missing. | plotly.js notice **added** (`THIRD_PARTY_NOTICES.md`). **The licence of BandQPI-3D itself must be chosen by the author** before publication (see the note at the end). |
| 3 | major | README cited `VALIDATION_v1.0.md`, missing from the zip. | **Added** to the package. |
| 4 | major | Bad option values gave Fortran/Python tracebacks; failures exited 0. | **Fixed**: every option is checked first with a clear message (exit 2); `--ewin` checked before the computation; unknown options rejected; the program exits 1 if the figures or the audit fail. |
| 5 | major | `--qpi` ignored for 1D. | **Fixed**: 1D uses the user's energies; for molecules a note says it does not apply. |
| 6 | major | `python3 bandqpi3d.py` without a result folder deleted files in `./figures/`. | **Fixed**: the plotter refuses any folder without `summary.txt` (tested: user files untouched). |
| 7 | major | `--fast` exceeded its time budget. | **Documented** as approximate (the budget sets grid sizes). |
| 8 | major | README overstated ("nothing interpolated", "every number has an error", "second program", "exact DOS", T-matrix naming). | **Rewritten** to state exactly what is done. |
| 9 | major | Test harness ignored each model's own validation/audit verdict. | **Fixed**: 20 extra checks (validation FAIL = 0, audit FAIL = 0 per model); suite now 90 checks; `install.sh --test` exits 1 on any failure and also runs the negative controls. |
| 10 | major | Program not found when started through PATH / a symlink. | **Fixed** (executable folder resolved through `/proc/<pid>/exe`; tested with a symlink on PATH). |
| — | minor | `negative_controls.py` did not compile on Python 3.11. | **Fixed**. |
| — | minor | README said eight models / a few minutes; options undocumented; file mode 0600. | **Fixed**. |
| — | minor | Some terminal lines > 80 columns with long absolute paths; banner needs colour. | Open (cosmetic). |

## R3 — figures and interactive page

All six pages: **0 JavaScript errors**, no page-level scrollbars, v1.0 logo and author line present. Gap values
in every figure agree with `summary.txt`. No gap arrow for graphene or the one-band metal; HfCrBr₂'s
gapless spin channel is labelled "metallic" and the system is named "half-metal".

| # | severity | finding | action |
|---|---|---|---|
| 1 | blocker | Non-magnetic systems showed spin-conserving / spin-flip QPI panels (flat 50 %) in figure 6c. | **Fixed**: figures show spin QPI only for magnets and spin–orbit runs. The page keeps the options (the author asked for them) and labels them "non-magnetic: = JDOS" / "exactly 1/2"; the spin column of the peak table is hidden. |
| 3 | blocker | For collinear magnets panel "JDOS" equals the spin-conserving map, while the caption wrote J with the total A. | The data are right (a non-magnetic impurity cannot flip spin, so J = Σ_σ A_σ⋆A_σ); **caption corrected** to that formula, and explained in README §9. |
| 2 | blocker | One-band metal: Born QPI at E_F (= its van Hove energy) is noisy. | Physics of a logarithmic singularity at that exact energy, not a bug; the default energies for real metals avoid it. Open (noted). |
| 4 | major | Legends covering the largest DOS peaks (fat bands, overview). | **Fixed**: legends placed where they hide no data. |
| 5 | major | Spin up drawn on the right in overview panel (d), on the left in figure 1. | **Fixed**: spin ↑ left everywhere. |
| 6 | major | Overview (c) of magnets: spin-↓ grey lines and the arrow unexplained. | **Fixed** in the caption. |
| 10 | major | Fat-band figures used thin grey high-symmetry lines instead of the purple arrows. | **Fixed**. |
| 12 | major | Header wrapped to three lines for magnets, so side panels scrolled at 1536 × 754. | **Fixed**: header kept to two lines (full subtitle as tooltip). |
| 7–9, 11 | major | Overview contour panel (a) for one-band and Dirac models; colour of a minor shell in CrCl₃ 8(c); caption step text; benzene fat bands by element only. | Open (overview figure 8 only; the primary figures 1, 1b, 1e are correct). |
| 13–21 | minor | k in Bohr⁻¹ on the band axis vs Å⁻¹ elsewhere (SIESTA's convention; `--kunits=ang` available), small label overlaps, wording. | Mostly open (cosmetic). |
| 22 | — | MoS₂ comes out indirect (VBM at Γ). | A property of the user's SIESTA calculation (reproduced exactly by R1 and by SIESTA's own `.bands`), not of BandQPI-3D. |

---

## Release decision

* Numbers: **no error found** by an independent re-implementation, by SIESTA's own output, or by the
  90-check exact test suite; the audit catches every deliberately injected error.
* All blockers that concern the code are fixed and re-tested; remaining open items are cosmetic and concern
  the overview figure 8, not the results.
* **One decision is the author's**: the software licence (e.g. GPL-3.0 like SIESTA, or MIT / BSD-3). Add it as
  `LICENSE` and to `CITATION.cff` before uploading to GitHub / Zenodo.
