# PISA for layered soils

Verification of the spring model `PISA` (layered sand/clay profiles)
against Burd et al. (2020b), *Application of the PISA design model to
monopiles embedded in layered soils*, Géotechnique 70(11), 1067–1082.

## Reference values

The paper compares its 1D model with 3D finite-element analyses for 35
layered profiles. Table 4 gives the 3D loads H(3D) and the ratio
ρ = H(1D)/H(3D) at v_G = D/10 000 (small displacement, "sd") and D/10
(ultimate, "ult"). The paper's own 1D result therefore follows as
H(1D) = ρ·H(3D). ρ is given to two decimals, so the reference values carry
a rounding uncertainty of about ±0.5 %.

Two groups of cases were selected:

**Sand only:** inputs exactly defined by formulas.
- Dunkirk sand (GDSM), γ′ = 10.09 kN/m³, K0 = 0.4.
- G0 after Hardin & Black with B = 875 and e0 linear in D_R (0.741 at 45 %
  to 0.573 at 90 %, Burd et al. 2020a, Table 2).

**Sand over Cowden till:**
- The Cowden layer as in the PISA clay verification: G0 = 1100·p′ with K0
  from Fig. 19, su read from Byrne et al. (2020a), Fig. 4a.
- As in the paper, su is a function of p′: in the layered profile, su is
  taken from the homogeneous Cowden profile at the same mean effective
  stress.

| Case | Pile | Profile |
|---|---|---|
| A4 | D2t | Dunkirk D_R 45 % over 90 % (boundary at L/2) |
| E1A | D2t | D_R 45 % (L/4), 75 %, 90 % (from L/8 above the tip) |
| E1B | D2t | D_R 55 % (L/4), 70 %, 85 % (from L/8 above the tip) |
| F6 | C1t | E1A profile, L/D = 2 |
| A2 | D2t | Dunkirk D_R 45 % over Cowden till (boundary at L/2) |
| F1 | C1t | A2 profile, L/D = 2 |

Piles (Table 1): D2t: D = 8.75 m, t = 150 mm, L = 35 m, h = 87.5 m;
C1t: D = 10 m, t = 150 mm, L = 20 m, h = 50 m. Settings as in the other PISA
verifications: `DIAMETER_REFERENCE = "mid"`, `SHEAR_FACTOR = 0.5`,
`PISA_STEEL`, 0.5 m elements, displacement control (D/10 in 5 levels,
D/10 000 in 1 level), values interpolated to exactly D/10 and D/10 000.

## Results

| Case | H_sd this impl. | H_sd = ρ·H(3D) | Diff. | H_ult this impl. | H_ult = ρ·H(3D) | Diff. |
|---|---|---|---|---|---|---|
| A4 | 402.6 kN | 392.6 kN | +2.5 % | 38.72 MN | 38.48 MN | +0.6 % |
| E1A | 408.8 kN | 401.0 kN | +1.9 % | 44.90 MN | 44.70 MN | +0.4 % |
| E1B | 417.8 kN | 409.2 kN | +2.1 % | 43.97 MN | 43.81 MN | +0.4 % |
| F6 | 558.2 kN | 554.0 kN | +0.8 % | 24.83 MN | 24.78 MN | +0.2 % |
| A2 | 410.9 kN | 400.8 kN | +2.5 % | 26.16 MN | 26.40 MN | −0.9 % |
| F1 | 633.4 kN | 640.4 kN | −1.1 % | 13.33 MN | 13.52 MN | −1.5 % |

All runs finished, all levels converged.

## Assessment

- **Ultimate loads, sand only:** within 0.6 %, i.e. within the rounding of
  ρ plus the 0.2 % accuracy of the homogeneous verification. The local
  application of the sand curves per layer reproduces the paper's model.
- **Ultimate loads, sand over Cowden till:** −0.9 % and −1.5 %. This is
  consistent with the slightly too low Cowden su of the figure readings
  (PISA clay verification: −1.8 … −4 %), which here only acts in the lower
  half of the pile.
- **Small displacement, pile D2t (A4, E1A, E1B, A2):** about +2 %, more
  than the rounding of ρ. This is **unexplained**.
  - The same model reproduces Table 7 of the sand paper within 0.1 % at
    1 mm.
  - F6, the same profile on the C1t pile with the same wall thickness,
    deviates only +0.8 %.
  - Treating D = 8.75 m as the outer diameter would remove the deviation
    (A4 −0.3 %, E1A −0.9 % at small displacement). It is not adopted:
    with the same convention, the exact values of Table 7 for C1 and C4
    are missed by 1.2–1.6 % (see the PISA sand verification, "Pile
    diameter"). So no single diameter convention fits both papers.
- **Small displacement, F1 (sand over Cowden till, L/D = 2):** −1.1 %,
  within the uncertainty of the Cowden G0 (K0 from the figure).

Overall, the layered implementation reproduces the paper's 1D model within
1 % at ultimate displacement for profiles with exactly defined inputs, and
within about 2.5 % at small displacement.

## Repeating

Input folders `Inputs/Verification/PISA_Layered_Verification/<case>_ultimate`
and `<case>_small_displacement`, e.g.

```
run_abaqus.bat Inputs\Verification\PISA_Layered_Verification\E1A_ultimate
```

Durations: about 0.5 min per small-displacement run and 3–6.5 min per
ultimate run (iterative sand model).
