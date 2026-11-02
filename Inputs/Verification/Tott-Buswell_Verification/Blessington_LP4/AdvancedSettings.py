"""Advanced settings of the automated monopile model builder.

Settings that are rarely changed. The settings of a single calculation are
in Input.py. An input folder without its own AdvancedSettings.py uses this
file (Inputs/AdvancedSettings.py).

Contents:
    1. Pile model
    2. Spring iteration
    3. Spring curve sampling
"""

# =========================================================================== #
# 1. Pile model
# =========================================================================== #
DIAMETER_REFERENCE = "outer"   # D in pile.csv: "outer" diameter or "mid"
                               # (mid-surface; use "mid" with the D of the
                               # PISA papers to reproduce their results)
SHEAR_FACTOR = 1.0             # transverse shear factor kappa of the pile
                               # (None = Abaqus default 0.53, PISA: 0.5,
                               # IGtHPile: 1.0)

# =========================================================================== #
# 2. Spring iteration (iterative spring models: PISA, THIEKEN)
# =========================================================================== #
ITERATION_TOLERANCE = 1e-3     # max. relative residual of v and psi along the
                               # pile (result vs. state of the springs)
RELAXATION = 0.5               # initial share of the new displacements used
                               # for the next spring update; adapted
                               # automatically (Aitken, 0.1 ... 1.0)

# =========================================================================== #
# 3. Spring curve sampling
# =========================================================================== #
# Every spring curve is handed to Abaqus as a table: x = 0 plus
# CURVE_POINTS - 1 logarithmically spaced points from GRID_START * x_max to
# x_max (x_max per model, see methods/spring_models), plus the kinks of the
# model curves. The table is mirrored for negative displacements.
CURVE_POINTS = 30              # points per curve, including x = 0 (>= 3)
GRID_START = 1e-4              # first point after 0 as a fraction of x_max
