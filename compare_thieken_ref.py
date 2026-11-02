"""Compare a THIEKEN run of the reference system with Thieken et al. (2015).

Reference system (geotechnik 38(4), sections 4.2 and 7.2): D = 5 m,
L = 30 m, t = 68.85 mm, h = 25 m, very dense sand (Table 1). Input folders:
Inputs/Verification/Thieken_Verification/offshore_*.

Reference values: FE secant stiffness K = H / y at ground level, read from
Fig. 17(1) at phi = 40 deg (approx. +-5-10 %); the paper's p-y approach lies
within 0.96 ... 1.04 of these values.

Run with plain Python after the Abaqus run:
    python compare_thieken_ref.py [path to pile_head.csv]
(default: results/verification/Thieken_Verification/offshore_30MN/
pile_thieken_ref/soil_thieken_ref/Thieken/pile_head.csv).
Assumes LOAD_ANGLE = 0 and the pile head
at the mudline (MUDLINE_DEPTH = 0), so U1 is the ground-level displacement.
"""

import csv
import math
import os
import sys

D = 5.0
# (y / D, K_FEM [MN/m])
REFERENCE = ((0.0005, 220.0), (0.001, 195.0), (0.003, 155.0), (0.01, 120.0),
             (0.03, 92.0), (0.1, 60.0))


def reference_stiffness(y):
    """K_FEM [MN/m] at displacement y [m], log-linear interpolation in y;
    None outside the reference range."""
    points = [(math.log(r * D), k) for r, k in REFERENCE]
    x = math.log(y)
    for (x0, k0), (x1, k1) in zip(points[:-1], points[1:]):
        if x0 <= x <= x1:
            return k0 + (k1 - k0) * (x - x0) / (x1 - x0)
    return None


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "results", "verification",
        "Thieken_Verification", "offshore_30MN", "pile_thieken_ref",
        "soil_thieken_ref", "Thieken", "pile_head.csv")
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))

    print("%10s %10s %12s %12s %8s" % ("H [MN]", "y [mm]", "K [MN/m]",
                                       "K_FEM [MN/m]", "K/K_FEM"))
    for row in rows:
        h = float(row["H [N]"]) / 1e6
        y = float(row["U1 [m]"])
        if y <= 0.0:
            continue
        k = h / y
        k_ref = reference_stiffness(y)
        print("%10.3f %10.2f %12.1f %12s %8s"
              % (h, 1000.0 * y, k,
                 "-" if k_ref is None else "%.1f" % k_ref,
                 "-" if k_ref is None else "%.3f" % (k / k_ref)))
    print("\nPaper: Thieken approach / FEM = 0.96 ... 1.04 (Fig. 17); "
          "reference values read from the plot (+-5-10 %).")


main()
