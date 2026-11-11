# BandQPI-3D v1.0 — validation record

Date: 11 November 2026. Every number below was produced by the released code (`siesta3d.f90`, `bandqpi3d.py`
of this package) and can be reproduced: part A with `bash install.sh --test`, part B by running the code on
the listed SIESTA/HONPAS outputs, part C with `tests/negative_controls.py`.

Accuracy rule used throughout: relative error ≤ 1 % (absolute limits in eV are given where a relative
error is meaningless, e.g. for a zero gap).

---

## A. Exact answers: ten synthetic tight-binding models (90 checks, 90 PASS)

`tests/make_models.py` writes genuine SIESTA HSX (version 2) files and fdf inputs. The reference values are
computed by `tests/check_models.py` from the model definition — its own Bloch sum and diagonalisation,
dense-grid state counting, and closed-form formulas — not from the HSX file or from BandQPI-3D.

| model | dim | key checks (code vs exact) |
|---|---|---|
| A square lattice, one band, t = −1 eV | 2D metal | path bands max \|ΔE\| 6×10⁻⁶ eV; band edges ±4\|t\|; E_F = 0; **DOS vs K(m)/(2π²t): median rel. error 1.2×10⁻⁴**, 95th percentile 8×10⁻⁴ |
| B staggered honeycomb, Δ = 1 eV | 2D insulator | gap 2.000000 eV (exact 2Δ), direct at K, VBM −1, CBM +1; integrated DOS error 9×10⁻⁵ |
| C graphene (Δ = 0), **automatic k-path** | 2D Dirac | gap 7.5×10⁻⁵ eV → "Dirac semimetal", **no gap arrow**; v_F = 8.92185×10⁵ m/s vs 3\|t\|d/2ħ = 8.92190×10⁵ (5×10⁻⁶) |
| D SSH chain, t₁ = −2, t₂ = −1.2 eV | 1D | gap 1.600000 eV = 2\|t₁−t₂\|; DOS 2.5×10⁻⁵ |
| E CsCl two-orbital crystal | 3D insulator | gap 1.600000 eV = 2Δ on the zone face; DOS 5.3×10⁻³ |
| F spin-split honeycomb | 2D half-metal | spin ↑ metallic, spin ↓ gap 3.000000 eV; electrons per spin and moment vs exact count < 1 % |
| G ferromagnetic insulator | 2D magnet | gaps ↑ 6.418861 / ↓ 1.200000 eV (exact); moment 1.000 μB; E_F inside the gap |
| H benzene Hückel ring | 0D | all six levels 2t cos(2πm/6) to 4×10⁻¹⁵ eV; HOMO–LUMO 5.000 eV |
| I = B with a path that misses K | 2D | gap still 2.000 eV (from the whole zone); figure states "edges off this path (on path 5.758 eV)" |
| J = B in the y–z plane, vacuum = first lattice vector | 2D | identical numbers to B: periodicity and plane are detected, not assumed |

In addition, each model's own validation report and independent audit must show zero FAIL (20 checks).
The full table (each check, value, exact value, error, limit) is written to `tests/run/model_check.txt`.
Fresh-folder run of `bash install.sh --test` on 11 Nov 2026: 90 PASS, 0 FAIL, negative controls all detected, 17.7 min on 2 cores.

## B. Real SIESTA / HONPAS calculations (seven systems)

Each run ends with `validation_report.txt` (engine and figure checks) and an independent audit.

| system | type | gap (eV) | vs SIESTA .bands (all bands) | DOS vs SIESTA .DOS | DOS k-grid error | validation | audit |
|---|---|---|---|---|---|---|---|
| MoS₂ monolayer | 2D, non-magnetic | 1.594963 indirect | 5×10⁻⁵ eV | no .DOS file | 0.04 % | ALL PASSED | 44 PASS 0 FAIL |
| graphene, HSE06 (HONPAS, old HSX) | 2D Dirac | 0.000 (no arrow) | 1.6×10⁻⁴ eV (3.4×10⁻⁴ relative, all bands) | 0.31 % | 0.03 % | ALL PASSED | 39 PASS 0 FAIL |
| TiO₂ nanowire | 1D | 0.821914 indirect | 5×10⁻⁵ eV | 0.43 % | 0.01 % | ALL PASSED | 41 PASS 0 FAIL |
| Si (8 atoms) | 3D | 0.474994 indirect | 6×10⁻⁵ eV | no .DOS file | 0.39 % | ALL PASSED | 42 PASS 0 FAIL |
| benzene | 0D | 5.131 HOMO–LUMO | levels vs .EIG 5×10⁻⁶ eV | — | — | ALL PASSED | 18 PASS 0 FAIL |
| HfCrBr₂ | 2D half-metal, 4.00 μB | ↑ metallic, ↓ 3.219 | 5×10⁻⁵ eV | 0.44 % | 0.001 % | ALL PASSED | 40 PASS 0 FAIL |
| CrCl₃ | 2D magnet + spin–orbit, 6.00 μB | 1.719 (↑-like), 3.182 (↓-like) | 5×10⁻⁵ eV | 0.34 % | 0.01 % | ALL PASSED | 48 PASS 0 FAIL |

## C. Negative controls (does the audit catch errors?)

`tests/negative_controls.py` copies a finished result, damages it, and requires the audit to fail.
On the honeycomb model and on MoS₂:

| damage | caught by |
|---|---|
| DOS scaled by 5 % | sum of PDOS ≠ total DOS, DOS vs SIESTA, electron count |
| 10 % of the occupied DOS removed | same |
| one band shifted by 50 meV on the path | path energies vs SIESTA .bands, all bands (B03); gap arrows (H02) |
| reported gap 3 % too large | gap on the periodic grid vs reported gap (C01) |
| Fermi level moved by 0.1 eV | E_F vs SIESTA .bands (B01), gap arrows (H02) |

The first MoS₂ run of this test found that a shift of the *lowest* (semicore) band went undetected because
the comparison with SIESTA was limited to E_F ± 10 eV; check B03 (all bands) was added and now catches it.

## D. Independent re-computation by a reviewer

A reviewer wrote a separate HSX reader and Bloch-sum diagonaliser (Python/SciPy), first confirmed it against
SIESTA's own `.EIG` (≤ 1.1×10⁻⁴ eV), then reproduced, for MoS₂, HfCrBr₂ and Si: the path eigenvalues
(≤ 8×10⁻⁵ eV), every VBM / CBM / gap to 10⁻⁶ eV, the electron counts and the magnetic moment
(4.00001 vs 3.99979 μB), and windowed DOS integrals (≤ 0.07 % after grid extrapolation). See
`REVIEW_PANEL_v1.0.md`.
