# Verification case: Thieken_Verification/offshore_30MN
# Monopile reference system (D = 5 m, L = 30 m, h = 25 m, very dense sand),
# H = 3 ... 30 MN: load-displacement curve (Fig. 17)
# Reference: Thieken, Achmus & Lemke (2015), A new static p-y approach for piles with
# arbitrary dimensions in sand, geotechnik 38(4)
# Results and evaluation: results/verification/Thieken_Verification/Thieken_Verification.md

"""Input of the automated monopile model builder.

This is the file to edit for a calculation; rarely changed settings
(diameter reference, shear factor, iteration tolerance, relaxation, curve
sampling) are in AdvancedSettings.py. Both belong to an input folder
together with the pile and soil CSV files (default: Inputs; verification
cases: Inputs/Verification/...). Start the calculation with
run_abaqus.bat [input folder] (input folder default: INPUT_DIR in the .bat).

Contents:
    1. Files
    2. Pile model
    3. Soil springs
    4. Load
    5. Job and output
    6. Column numbers of soil.csv (per spring model)
"""

# =========================================================================== #
# 1. Files
# =========================================================================== #
# PILE_CSV and SOIL_CSV are relative to the folder of this Input.py (or
# absolute paths); RESULTS_DIR is relative to the project folder (main.py).
# PILE_CSV, SOIL_CSV and SPRING_MODEL (section 3) can each be a single value
# or a list, e.g. PILE_CSV = ["pile_C1.csv", "pile_C4.csv"]. Every
# combination pile x soil x spring model is calculated in its own folder
# RESULTS_DIR/<pile>/<soil>/<model> (file names without .csv), e.g.
# results/pile_thieken_onshore/soil_thieken_onshore/Thieken. With more than
# one combination, RESULTS_DIR/summary.csv lists the pile head results,
# max. M and Q of all of them.
PILE_CSV = "pile_thieken_ref.csv"  # pile geometry
SOIL_CSV = "soil_thieken_ref.csv"  # soil layers
RESULTS_DIR = "results/verification/Thieken_Verification/offshore_30MN" # output folder

# =========================================================================== #
# 2. Pile model
# =========================================================================== #
MODEL_NAME = "Monopile"
ELEMENT_COUNT = 60             # number of beam elements over the whole pile
                               # (head to tip); element size ~ pile length /
                               # ELEMENT_COUNT. Nodes are always placed at
                               # pile section and soil layer boundaries and
                               # at the mudline.
ELEMENT_TYPE = "B31"           # B31, B32, PIPE31, PIPE32
MUDLINE_DEPTH = 0.0            # depth of the mudline below the pile head [m]

# =========================================================================== #
# 3. Soil springs
# =========================================================================== #
SPRING_MODEL = "THIEKEN"
                               # "API", "KALLEHAVE", "SORENSEN", "THIEKEN"
                               # (p-y) or "PISA_SAND" (multi-spring), or a
                               # list of them
LOADING = "static"             # API, KALLEHAVE, SORENSEN: "static" or
                               # "cyclic"
MAX_ITERATIONS = 20            # iterative spring models (PISA, THIEKEN)

# =========================================================================== #
# 4. Load at pile head
# =========================================================================== #
LOAD_CONTROL = "force"         # "force": load levels of HORIZONTAL_LOAD
                               # "displacement": load levels of
                               # TARGET_DISPLACEMENT at the mudline; H is
                               # found iteratively for every level
HORIZONTAL_LOAD = 30e6         # [N] force control: maximum load;
                               # displacement control: estimate of H at
                               # TARGET_DISPLACEMENT (start value only)
TARGET_DISPLACEMENT = 0.5      # [m] displacement control: maximum mudline
                               # displacement (e.g. 0.1 * D)
LOAD_ANGLE = 0.0               # direction in the XY plane, from X towards Y
                               # [deg]
LEVER_ARM = 25.0               # height of the load above the pile head [m];
                               # applied as H plus moment H * LEVER_ARM
LOAD_INCREMENT = 0.1           # step width of the load levels as fraction
                               # of HORIZONTAL_LOAD / TARGET_DISPLACEMENT
MAX_INCREMENTS = 10            # maximum number of load increments

# Break criteria at the mudline (None = not used). The loading stops after
# the first load level that reaches a criterion; that level is kept, and the
# load at the criterion (interpolated between the last two levels) is
# written to progress.log and summary.csv.
MAX_DISPLACEMENT = None        # [m] mudline displacement, e.g. 0.1 * D
MAX_ROTATION = None            # [deg] mudline rotation, e.g. 0.5
VERTICAL_TANGENT_RATIO = None  # [-] |rotation at the toe| / |rotation at the
                               # mudline|, e.g. 0.01: the deflection line no
                               # longer has a vertical tangent at the toe

# =========================================================================== #
# 5. Job and output
# =========================================================================== #
# The Abaqus job of every combination is named <pile>_<soil>_<model>, e.g.
# pile_thieken_ref_soil_thieken_ref_Thieken (characters other than letters,
# digits and _ replaced by _, at most 80 characters).
# Every folder RESULTS_DIR/<pile>/<soil>/<model> contains after the run:
#   <job>.odb         results of the last converged load increment
#   <job>.cae         model with the springs and load of that increment
#   pile_head.csv     every increment: load, pile head displacements and
#                     rotations, mudline and toe displacement/rotation
#   pile_profile.csv  last converged increment, every pile node:
#                     displacements, rotations, section forces N, Q, M
#   load_displacement.png  H over mudline displacement and rotation (PLOTS)
#   pile_profile.png       last increment: deflection, Q and M over depth
# and, with several combinations, RESULTS_DIR/load_displacement_all.png.
NUM_CPUS = 1
PLOTS = True                   # create the PNG plots

# =========================================================================== #
# 6. Column numbers of soil.csv
# =========================================================================== #
# For every soil parameter, give the number of the column in soil.csv that
# contains it, counted from 1 (first column = 1, as in Excel: A = 1, B = 2).
# The first line of soil.csv is a header and is skipped; its names do not
# matter.
#
#     "parameter": 5      value in column 5
#     "parameter": None   not contained in soil.csv
#
# main.py uses COMMON_COLUMNS plus the block of the spring model selected
# with SPRING_MODEL (see MODEL_COLUMNS at the end). Columns of soil.csv that
# the selected model does not use are listed as a warning.
#
# Details on the parameters of each model are in the "Soil parameters
# (soil.csv)" section of the model in methods/spring_models.

# --------------------------------------------------------------------------- #
# All models (required)
# --------------------------------------------------------------------------- #
COMMON_COLUMNS = {
    "z_top": 1,         # top of layer, depth below mudline [m]
    "z_bot": 2,         # bottom of layer, depth below mudline [m]
    "soil_type": 3,     # "sand" or "clay"
    "gamma_eff": 4,     # effective unit weight [kN/m^3]
}

# --------------------------------------------------------------------------- #
# API p-y curves
# --------------------------------------------------------------------------- #
API_COLUMNS = {
    # sand
    "phi": 5,           # friction angle [deg]
    "k": 9,             # subgrade modulus [MN/m^3] (optional)
    # clay (Matlock)
    "su": 6,            # undrained shear strength [kPa]
    "eps50": 7,         # strain at 50 % of max. stress [-] (optional)
    "j": 8,             # Matlock constant J [-] (optional)
}

# --------------------------------------------------------------------------- #
# Kallehave et al. (2012): API p-y, E_py = k z0 (z/z0)^0.6 (D/D0)^0.5
# --------------------------------------------------------------------------- #
KALLEHAVE_COLUMNS = {
    # sand
    "phi": 5,           # friction angle [deg]
    "k": 9,             # subgrade modulus [MN/m^3] (optional)
    # clay (Matlock)
    "su": 6,            # undrained shear strength [kPa]
    "eps50": 7,         # strain at 50 % of max. stress [-] (optional)
    "j": 8,             # Matlock constant J [-] (optional)
}

# --------------------------------------------------------------------------- #
# Sorensen et al. (2010): API p-y, E_py = a (z/zref)^0.6 (D/Dref)^0.5 phi^3.6
# --------------------------------------------------------------------------- #
SORENSEN_COLUMNS = {
    # sand
    "phi": 5,           # friction angle [deg]
    # clay (Matlock)
    "su": 6,            # undrained shear strength [kPa]
    "eps50": 7,         # strain at 50 % of max. stress [-] (optional)
    "j": 8,             # Matlock constant J [-] (optional)
}

# --------------------------------------------------------------------------- #
# PISA sand model
# --------------------------------------------------------------------------- #
PISA_SAND_COLUMNS = {
    # sand
    "dr": 10,           # relative density [%]
    "phi": 5,           # friction angle [deg] (optional, k0 = 1 - sin phi
                        # for G0; without phi k0 = 0.4)
    # G0, one of: g0 / g0_ref + lambda_g0 / hb_b + e0
    "g0": 11,           # G0 directly [MPa], "top-bottom" allowed
    "g0_ref": None,     # power law: G0 at 100 kPa mean stress [MPa]
    "lambda_g0": None,  # power law: stress exponent of G0 [-]
    "hb_b": 12,         # Hardin & Black: constant B [-]
    "e0": 13,           # Hardin & Black: void ratio [-]
    # clay (API Matlock fallback)
    "su": 6,            # undrained shear strength [kPa]
    "eps50": 7,         # strain at 50 % of max. stress [-] (optional)
    "j": 8,             # Matlock constant J [-] (optional)
}

# --------------------------------------------------------------------------- #
# Thieken et al. (2015)
# --------------------------------------------------------------------------- #
THIEKEN_COLUMNS = {
    # sand
    "phi": 5,           # friction angle [deg]
    "nu": 14,           # Poisson's ratio [-]
    "es_ref": 15,       # Es at 100 kPa mean stress [MPa]
    "lambda_es": 16,    # stress exponent of Es [-]
    # G0, one of: g0 / g0_ref + lambda_g0 / hb_b + e0
    "g0": 11,           # G0 directly [MPa], "top-bottom" allowed
    "g0_ref": 17,       # power law: G0 at 100 kPa mean stress [MPa]
    "lambda_g0": 18,    # power law: stress exponent of G0 [-]
    "hb_b": 12,         # Hardin & Black: constant B [-]
    "e0": 13,           # Hardin & Black: void ratio [-]
    # clay (API Matlock fallback)
    "su": 6,            # undrained shear strength [kPa]
    "eps50": 7,         # strain at 50 % of max. stress [-] (optional)
    "j": 8,             # Matlock constant J [-] (optional)
}

# --------------------------------------------------------------------------- #
# Spring model (SPRING_MODEL in main.py) -> column block
# --------------------------------------------------------------------------- #
MODEL_COLUMNS = {
    "API": API_COLUMNS,
    "KALLEHAVE": KALLEHAVE_COLUMNS,
    "SORENSEN": SORENSEN_COLUMNS,
    "PISA_SAND": PISA_SAND_COLUMNS,
    "THIEKEN": THIEKEN_COLUMNS,
}
