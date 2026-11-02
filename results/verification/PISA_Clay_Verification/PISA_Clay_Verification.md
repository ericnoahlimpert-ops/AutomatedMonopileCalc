# PISA clay model: piles C1 and C4

Verification of the spring model `PISA_CLAY` (Cowden till calibration)
against Byrne et al. (2020), *PISA design model for monopiles for offshore
wind turbines: application to a stiff glacial clay till*, Géotechnique
70(11), 1030–1047.

## Reference

Table 5 of the paper gives exact results of the paper's own 1D model with
the second-stage parameters (Table 4) for the calibration piles C1 and C4,
using its finest meshes:

| Pile | D | L | t | h | H_ult at v_G = D/10 | H_sd at v_G = D/10 000 |
|---|---|---|---|---|---|---|
| C1 | 10 m | 20 m | 91 mm | 50 m | 12.9636 MN | 0.7395 MN |
| C4 | 10 m | 60 m | 91 mm | 50 m | 104.0325 MN | 0.9105 MN |

## Inputs

**The soil profile of the representative offshore Cowden till is only
given as figures**, so part of the input is read from them:

- **σ′_v:** exact, from the bulk unit weight 21.19 kN/m³ and hydrostatic
  pore pressure, γ′ = 11.38 kN/m³ (paper; Burd et al. 2020b, Table 13).
- **G0 = 1100·p′:** exact relation (Byrne et al. 2020, Table 2:
  G0* = 110 MPa at p′ref = 100 kPa, linear in p′, as in Zdravković et al.
  2020), with p′ = σ′_v·(1 + 2K0)/3.
- **K0:** read from Fig. 19 of Burd et al. (2020b), the clearest plot of
  the profile: 1.5 down to 5 m (the cap described by Zdravković et al.
  2020, Fig. 5), then linear to 0.98 at 13 m, 0.84 at 20 m, 0.75 at 30 m,
  0.70 at 40 m and 0.67 at 50 m.
- **su:** read from Fig. 4a of the paper: 50 kPa at the mudline, peak of
  155 kPa at 3.75 m, 105 kPa at 5.5 m, then linear to 250 kPa at 50 m and
  further at the same gradient. As the su profile of the paper is the
  output of the 3D constitutive model, it cannot be derived from the model
  parameters (Table 2).
- **soil.csv:** su and G0 as piecewise linear layers of 2.5 m (plus the
  kinks of su).

Model settings:

- `DIAMETER_REFERENCE = "mid"`, `SHEAR_FACTOR = 0.5`, steel E = 200 GPa
  (`PISA_STEEL`), as in the paper.
- 0.5 m elements.
- Displacement control to v_G = 1 m (5 levels) and to 1 mm (1 level).

## Results

| Pile | Value | This implementation | Table 5 | Difference |
|---|---|---|---|---|
| C1 | H_ult at v_G = 1 m | 12.440 MN | 12.964 MN | −4.0 % |
| C4 | H_ult at v_G = 1 m | 102.19 MN | 104.03 MN | −1.8 % |
| C1 | H_sd at v_G = 1 mm | 0.7203 MN | 0.7395 MN | −2.6 % |
| C4 | H_sd at v_G = 1 mm | 0.8903 MN | 0.9105 MN | −2.2 % |

The values are interpolated linearly to exactly v_G = 1 m / 1 mm from the
last load level, as for sand. All levels converged.

**Assessment.** All four values are 2–4 % low. The most likely cause is the
part of the soil profile that is read from figures, not the
implementation:

- **The same code reproduces the sand model within 0.2 %** (Burd et al.
  2020a, Table 7 and Fig. 13), where the soil input is given by formulas.
  That covers the conic function, the diameter, shear factor, springs,
  displacement control and evaluation. The clay model only replaces the
  normalisation (su, G0) and the parameter functions of Table 4, which were
  checked term by term.
- **H_sd does not depend on su at all:** in the linear range p = k_p·G0·v
  and m = k_m·G0·D²·ψ. With G0 = 1100·p′ exact, the remaining −2 % comes
  from K0, the only figure reading in G0.
- **H_ult is governed by su** (p_u = s_u·D·p̄_u). The −2 … −4 %
  corresponds to su read slightly too low in the upper part of the profile;
  C1 (L = 20 m) depends most on the near-surface peak.

**Effect of the K0 reading.** A first evaluation used K0 read from the
small Fig. 4b of the paper (dropping to 0.90 already at 10 m). G0 was
then about 15 % too low between 5 and 15 m, and the results were −7.1 %
(C1) and −5.1 % (C4) at 1 mm, and −4.7 % / −2.2 % at 1 m. The K0 profile
of Burd et al. (2020b), Fig. 19, roughly halved the small-displacement
deviation.

**Plausibility check.** M at the pile head equals H·h within 0.07 %; the
rest is the moment of the rotational spring at the head node (the head lies
at the mudline).

## Repeating

```
run_abaqus.bat Inputs\Verification\PISA_Clay_Verification\C1_ultimate
```

| Run | Duration |
|---|---|
| `C1_ultimate` | ≈ 3.5 min |
| `C4_ultimate` | ≈ 3 min |
| `C1_small_displacement` | ≈ 0.2 min |
| `C4_small_displacement` | ≈ 0.2 min |

References: Burd et al. (2020b), *Application of the PISA design model to
monopiles embedded in layered soils*, Géotechnique 70(11), 1067–1082;
Zdravković et al. (2020), *Ground characterisation for PISA pile testing
and analysis*, Géotechnique 70(11), 945–960.
