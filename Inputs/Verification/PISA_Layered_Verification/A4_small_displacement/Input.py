# Verification case: PISA_Layered_Verification/A4_small_displacement
# Case A4, pile D2t: Dunkirk sand Dr 45 % over Dr 90 % (boundary at L/2),
# displacement control to v_G = 0.000875 m
# Reference: Burd et al. (2020b), Application of the PISA design model to
# monopiles embedded in layered soils, Geotechnique 70(11), Table 4:
# H(1D) = rho * H(3D).
# Results and evaluation:
# results/verification/PISA_Layered_Verification/PISA_Layered_Verification.md

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
PILE_CSV = "pile_pisa_D2t.csv"  # pile geometry
SOIL_CSV = "soil_A4.csv"  # soil layers
RESULTS_DIR = "results/verification/PISA_Layered_Verification/A4_small_displacement"        # output folder
CPT_CSV = None                 # CPT profile for TOTT_BUSWELL: columns
                               # z [m below mudline], qc [MPa], optional
                               # g0 [MPa]; one file for all soil files or
                               # {"soil.csv": "cpt.csv", ...}

# =========================================================================== #
# 2. Pile model
# =========================================================================== #
MODEL_NAME = "Monopile"
ELEMENT_COUNT = 70             # number of beam elements over the whole pile
                               # (head to tip); element size ~ pile length /
                               # ELEMENT_COUNT. Nodes are always placed at
                               # pile section and soil layer boundaries and
                               # at the mudline.
ELEMENT_TYPE = "B31"           # B31, B32, PIPE31, PIPE32
MUDLINE_DEPTH = 0.0            # depth of the mudline below the pile head [m]

# =========================================================================== #
# 3. Soil springs
# =========================================================================== #
SPRING_MODEL = "PISA"
                               # "API", "KALLEHAVE", "SORENSEN", "THIEKEN"
                               # (p-y), "PISA_SAND", "PISA_CLAY", "PISA"
                               # (layered sand/clay) or "TOTT_BUSWELL"
                               # (multi-spring), or a list of them
LOADING = "static"             # API, KALLEHAVE, SORENSEN: "static" or
                               # "cyclic"
MAX_ITERATIONS = 30            # iterative spring models (PISA_SAND,
                               # THIEKEN) and displacement control
TOTT_BUSWELL_ALPHA = 0.1       # TOTT_BUSWELL: residual base stress
                               # q_b,res = TOTT_BUSWELL_ALPHA * qc at the
                               # tip (paper: 0.1, fitted 0.04 ... 0.1)
TOTT_BUSWELL_COMPONENTS = ("p-y", "m-theta", "base")
                               # TOTT_BUSWELL: springs used, e.g.
                               # ("p-y",) or ("p-y", "m-theta") as in the
                               # comparisons of the paper

# =========================================================================== #
# 4. Load at pile head
# =========================================================================== #
LOAD_CONTROL = "displacement"         # "force": load levels of HORIZONTAL_LOAD
                               # "displacement": load levels of
                               # TARGET_DISPLACEMENT at the mudline; H is
                               # found iteratively for every level
HORIZONTAL_LOAD = 0.3926e6         # [N] force control: maximum load;
                               # displacement control: estimate of H at
                               # TARGET_DISPLACEMENT (start value only)
TARGET_DISPLACEMENT = 0.000875      # [m] displacement control: maximum mudline
                               # displacement (e.g. 0.1 * D)
LOAD_ANGLE = 0.0               # direction in the XY plane, from X towards Y
                               # [deg]
LEVER_ARM = 87.5               # height of the load above the pile head [m];
                               # applied as H plus moment H * LEVER_ARM
LOAD_INCREMENT = 1           # step width of the load levels as fraction
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
# PISA clay model (Cowden till)
# --------------------------------------------------------------------------- #
PISA_CLAY_COLUMNS = {
    # clay
    "su": 6,            # undrained shear strength [kPa], "top-bottom"
    # G0, one of: g0 / g0_ref + lambda_g0 / hb_b + e0
    "g0": 11,           # G0 directly [MPa], "top-bottom" allowed
    "g0_ref": 17,       # power law: G0 at 100 kPa mean stress [MPa]
                        # (paper: G0 = 1100 p', i.e. 110 with lambda 1)
    "lambda_g0": 18,    # power law: stress exponent of G0 [-]
    "hb_b": 12,         # Hardin & Black: constant B [-]
    "e0": 13,           # Hardin & Black: void ratio [-]
    "phi": 5,           # friction angle [deg] (optional, k0 for G0)
    "pisa_clay": None,  # calibration: "cowden" (default) or "bothkennar"
}

# --------------------------------------------------------------------------- #
# PISA for layered sand/clay profiles (Burd et al. 2020b): sand layers as
# PISA_SAND, clay layers as PISA_CLAY
# --------------------------------------------------------------------------- #
PISA_COLUMNS = {
    # sand
    "dr": 10,           # relative density [%]
    # clay
    "su": 6,            # undrained shear strength [kPa], "top-bottom"
    "pisa_clay": None,  # calibration: "cowden" (default) or "bothkennar"
    # G0 (all layers), one of: g0 / g0_ref + lambda_g0 / hb_b + e0
    "g0": 11,           # G0 directly [MPa], "top-bottom" allowed
    "g0_ref": 17,       # power law: G0 at 100 kPa mean stress [MPa]
    "lambda_g0": 18,    # power law: stress exponent of G0 [-]
    "hb_b": 12,         # Hardin & Black: constant B [-]
    "e0": 13,           # Hardin & Black: void ratio [-]
    "phi": 5,           # friction angle [deg] (optional, k0 for G0)
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
# Tott-Buswell et al. (2024), CPT-based multi-spring model; qc (and
# optionally G0) from CPT_CSV
# --------------------------------------------------------------------------- #
TOTT_BUSWELL_COLUMNS = {
    # sand: G0 only if the CPT file has no g0 column, one of:
    # g0 / g0_ref + lambda_g0 / hb_b + e0
    "g0": 11,           # G0 directly [MPa], "top-bottom" allowed
    "g0_ref": 17,       # power law: G0 at 100 kPa mean stress [MPa]
    "lambda_g0": 18,    # power law: stress exponent of G0 [-]
    "hb_b": 12,         # Hardin & Black: constant B [-]
    "e0": 13,           # Hardin & Black: void ratio [-]
    "phi": 5,           # friction angle [deg] (optional, k0 for G0)
    "dr": 10,           # relative density [%] (optional, base shear;
                        # default 75 %)
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
    "PISA_CLAY": PISA_CLAY_COLUMNS,
    "PISA": PISA_COLUMNS,
    "THIEKEN": THIEKEN_COLUMNS,
    "TOTT_BUSWELL": TOTT_BUSWELL_COLUMNS,
}
