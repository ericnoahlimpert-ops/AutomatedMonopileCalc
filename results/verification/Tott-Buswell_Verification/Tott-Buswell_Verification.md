# Tott-Buswell model: Blessington piles

Plausibility check of the spring model `TOTT_BUSWELL` against Tott-Buswell,
Prendergast & Gavin (2024), *A CPT-based multi-spring model for lateral
monopile analysis under SLS conditions in sand*, Ocean Engineering 293,
116642.

**This is not an exact verification.** The paper gives neither tabulated
model results nor the CPT data ("The authors do not have permission to share
data"). The qc and G0 profiles were read approximately from Fig. 4, and the
reference loads from the small insets of Fig. 5. Agreement within about
±15 % is the best that can be expected.

## Cases

Piles from Table 2, dry dense sand, γ = 19.8 kN/m³, D_r ≈ 100 %,
q_b,res = 0.1·q_c,r (TOTT_BUSWELL_ALPHA = 0.1), all four spring components.
Force control in 20 levels, pile head at the mudline, moment H·h.

| Run | Pile | L [m] | D [m] | t [mm] | h [m] | Quarry | H_max |
|---|---|---|---|---|---|---|---|
| `Blessington_LP2` | LP2 | 1.50 | 0.51 | 10 | 1.0 | lower | 100 kN |
| `Blessington_LP3` | LP3 | 2.25 | 0.51 | 10 | 1.0 | lower | 200 kN |
| `Blessington_LP4` | LP4 | 3.00 | 0.51 | 10 | 1.0 | lower | 300 kN |
| `Blessington_UP1` | UP1 | 4.50 | 0.34 | 14 | 0.4 | upper | 250 kN |

The inputs are in `Inputs/Verification/Tott-Buswell_Verification/<run>`:
- `cpt_blessington_lower.csv` / `cpt_blessington_upper.csv`: average qc
  (Fig. 4a/b) and G0 (Fig. 4c) read from the figure.
- Steel E = 200 GPa (`PISA_STEEL`).
- 0.05 m elements, as in the paper.

## Results

Load at given mudline displacements, interpolated from `pile_head.csv`:

| Pile | H at 0.005 D | H at 0.01 D | H at 0.02 D | Fig. 5, multi-spring 0.1 q_c,r, at 0.01 D ≈ |
|---|---|---|---|---|
| LP2 (L/D = 3) | 38.1 kN | 62.4 kN | > 100 kN (0.0195 D at 100 kN) | 55–60 kN |
| LP3 (L/D = 4.5) | 60.1 kN | 99.5 kN | 163.6 kN | 100–110 kN |
| LP4 (L/D = 6) | 70.8 kN | 121.0 kN | 204.3 kN | 115–130 kN |
| UP1 (L/D = 13) | 33.4 kN | 57.4 kN | 98.2 kN | 60–75 kN |

- **Up to 0.01 D (SLS range of the paper):** the results lie within the
  reading accuracy of Fig. 5.
- **LP2 at larger displacements:** the response is stiffer than the paper's
  curve. The paper reaches 100 kN at roughly 0.03 D; this implementation
  reaches it at 0.02 D. The short pile depends strongly on qc near the
  surface and on the base springs, and both come from figure readings here.
- **UP1:** 250 kN at 0.067 D, close to the end of the multi-spring curve in
  Fig. 5d.

**Plausibility checks.** At the pile head, M = H·h exactly (e.g. LP4:
299.98 kNm at 300 kN). The m-θ curve, eq. (7), matches a numerical
integration of eqs. (2) and (6), and θ_f and m_f match Table 1.

## Notes

- For UP1 (L/D = 13.2), q_c,r is averaged over ±0.25·L/D = ±3.3 m around the
  tip, i.e. down to 7.8 m. The CPT ends at 5.5 m, so qc is taken as constant
  below it (warning in the console). The paper mentions that this range
  becomes disproportionately large for high L/D; the base springs are
  unimportant there.
- The paper also shows "p-y only" and "p-y + m-θ" curves. They can be
  reproduced with `TOTT_BUSWELL_COMPONENTS = ("p-y",)` or `("p-y", "m-theta")`;
  use a different `RESULTS_DIR` so the full model's results are not
  overwritten.

## Repeating

```
run_abaqus.bat Inputs\Verification\Tott-Buswell_Verification\Blessington_LP2
```

(about 10 s per run, plus the Abaqus/CAE start).
