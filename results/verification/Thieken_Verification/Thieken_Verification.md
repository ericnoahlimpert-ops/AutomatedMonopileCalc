# Verification: Thieken et al. (2015)

Reference: Thieken, K., Achmus, M., Lemke, K. (2015). *A new static p-y
approach for piles with arbitrary dimensions in sand.* geotechnik 38(4),
267–288.

The paper gives its results only as figures. The reference values below were
read from the plots and are accurate to about ±5–10 % (marked ≈). The paper
states that its p-y approach lies within **0.96 … 1.04** of its 3D FE
results (Fig. 17).

**Result: the implementation reproduces the paper within the accuracy of the
figures. The secant stiffness of the monopile is 0.97 … 1.02 of the FE
values over the whole load range (2 mm … 500 mm), all load levels
converged.**

## Common settings

| Setting | Value |
|---|---|
| Spring model | `THIEKEN` |
| Load control | force |
| `DIAMETER_REFERENCE` | `"outer"` |
| `SHEAR_FACTOR` | 1.0 (as in IGtHPile, the program of the authors) |
| Element length | 0.5 m |
| Material | S355 (E = 210 GPa; paper: E = 210 GPa, ν = 0.27) |

## 1. Monopile reference system (sections 4.2 and 7.2 of the paper)

**Pile** (`pile_thieken_ref.csv`): D = 5 m, L = 30 m, t = 68.85 mm, load
eccentricity h = 25 m (`LEVER_ARM = 25`).

**Soil** (`soil_thieken_ref.csv`): very dense sand, Table 1 of the paper:

| γ' [kN/m³] | φ [°] | ν | E_s,ref [MPa] | λ_ES | G0,ref [MPa] | λ_G0 |
|---|---|---|---|---|---|---|
| 10.31 | 40 | 0.2 | 70 (κ = 700) | 0.5 | 106.3 (e = 0.60, eq. 10) | 0.5 |

### 1.1 Load-displacement behaviour (Fig. 17)

Runs `offshore_2.5MN` (0.5 … 2.5 MN) and `offshore_30MN` (3 … 30 MN).
Secant stiffness K = H/y at the mudline:

| H [MN] | y [mm] | K [MN/m] | K_FEM ≈ [MN/m] | K / K_FEM |
|---|---|---|---|---|
| 0.5 | 2.36 | 212.0 | 220 (extrapolated) | ≈ 0.96 |
| 1.0 | 5.34 | 187.2 | 192.6 | 0.972 |
| 1.5 | 8.69 | 172.6 | 174.9 | 0.987 |
| 2.0 | 12.33 | 162.2 | 162.1 | 1.001 |
| 2.5 | 16.21 | 154.2 | 152.7 | 1.010 |
| 3.0 | 20.32 | 147.6 | 146.2 | 1.010 |
| 6.0 | 48.97 | 122.5 | 120.6 | 1.016 |
| 9.0 | 83.74 | 107.5 | 106.9 | 1.006 |
| 12.0 | 124.29 | 96.5 | 96.8 | 0.997 |
| 15.0 | 170.58 | 87.9 | 88.6 | 0.993 |
| 18.0 | 222.64 | 80.8 | 81.5 | 0.992 |
| 21.0 | 280.71 | 74.8 | 75.3 | 0.993 |
| 24.0 | 345.02 | 69.6 | 69.9 | 0.996 |
| 27.0 | 416.10 | 64.9 | 64.9 | 1.000 |
| 30.0 | 494.10 | 60.7 | 60.3 | 1.007 |

K_FEM is interpolated (log-linear in y) from the values read from Fig. 17(1)
at φ = 40°: 220 / 195 / 155 / 120 / 92 / 60 MN/m at y = 0.0005 / 0.001 /
0.003 / 0.01 / 0.03 / 0.1·D. The comparison is produced by
`compare_thieken_ref.py`.

Iterations per level: 6 for the first level of each run, 4 for all others.
Additional check with Fig. 7: ≈ 160 mm at H = 15 MN (this implementation:
171 mm).

### 1.2 Profiles at H = 5 MN (Fig. 15)

Run `offshore_5MN` (2.5 and 5 MN).

| | This implementation | Fig. 15 (FE / Thieken) ≈ |
|---|---|---|
| Mudline displacement | 38.7 mm | 37–39 mm |
| Displacement at z = 10 m | 9.7 mm | 10 mm |
| Rotation point | 18.6 m | 18–20 m |
| Toe displacement | −2.8 mm | −3 mm |
| Max. bending moment | 140.5 MNm at z = 5.0 m | ≈ 140 MNm at ≈ 5 m |
| Moment at the mudline | 125.0 MNm (= 5 MN · 25 m) | 125 MNm |

Stretching factor (`pile_profile.csv`): 2.6 at the mudline, minimum 0.31 at
the rotation point (18.5 m), tip effect below L − 2D = 20 m up to 9.5 at the
toe — the shape of the SF panel in Fig. 15.

## 2. Onshore pile (section 7.3, Fig. 16)

Run `onshore_500kN`.

**Pile** (`pile_thieken_onshore.csv`): D = 1 m, L = 12.5 m, t = 16.35 mm,
no lever arm.

**Soil** (`soil_thieken_onshore.csv`): medium dense sand above the
groundwater table: γ = 19.76 kN/m³, φ = 35°, ν = 0.25, E_s,ref = 40 MPa
(κ = 400), λ_ES = 0.6, G0,ref = 89.4 MPa (e = 0.69), λ_G0 = 0.5.

| H [kN] | 100 | 200 | 300 | 400 | 500 |
|---|---|---|---|---|---|
| y [mm] | 1.56 | 4.21 | 7.57 | 11.54 | 16.06 |
| Iterations | 4 | 4 | 4 | 4 | 4 |

| At H = 500 kN | This implementation | Fig. 16 (FE / Thieken) ≈ |
|---|---|---|
| Mudline displacement | 16.1 mm | 15–17 mm |
| Rotation point | 6.7 m | 6.5–7 m |
| Max. bending moment | 960 kNm at 3.5 m | 900–950 kNm at ≈ 3.5 m |

## Input and result files

| Run | Content |
|---|---|
| `offshore_30MN` | 10 load levels 3 … 30 MN |
| `offshore_2.5MN` | 5 load levels 0.5 … 2.5 MN (small displacements) |
| `offshore_5MN` | profiles at 5 MN (Fig. 15) |
| `onshore_500kN` | 5 load levels 100 … 500 kN (Fig. 16) |

The input folder of every run is
`Inputs/Verification/Thieken_Verification/<run>` (`Input.py`,
`AdvancedSettings.py` and the CSV files); run it with
`run_abaqus.bat <input folder>`. The result folder `<run>/` contains
`progress.log` and `<pile>/<soil>/Thieken/` with the `.odb`, `.cae`,
`pile_head.csv`, `pile_profile.csv` and the plots `load_displacement.png`
and `pile_profile.png`. `python compare_thieken_ref.py` prints the table
of section 1.1 for `offshore_30MN` (pass the `pile_head.csv` of
`offshore_2.5MN` as argument for the small loads).

## Notes

- The paper's text gives y_A = 0.005·p_C/(E_s·SF), Fig. 14 an additional
  factor D; the Fig. 14 version is implemented.
- Section 7.4 of the paper gives t = 81.35 mm for the reference system,
  sections 7.2 and Fig. 17 give 68.85 mm; 68.85 mm is used.
- The shear factor of the pile is not stated in the paper; 1.0 is the value
  of IGtHPile.
