"""Automated 1D monopile model builder for Abaqus/CAE.

The settings of a calculation are in an input folder (default: Inputs):
Input.py, AdvancedSettings.py (rarely changed settings; if the folder has
none, the one of Inputs is used) and the pile/soil CSV files. The folder is
given by the environment variable MONOPILE_INPUT_DIR, which run_abaqus.bat
and run_abaqus_gui.bat set (INPUT_DIR setting or first argument):
    run_abaqus.bat [input folder]       without GUI
    run_abaqus_gui.bat [input folder]   in Abaqus/CAE
CSV paths in Input.py are relative to the input folder, RESULTS_DIR is
relative to the project folder (the folder of this script).

Every combination of pile, soil and spring model (PILE_CSV, SOIL_CSV and
SPRING_MODEL in Input.py may be lists) is calculated in its own folder
RESULTS_DIR/<pile>/<soil>/<model>, e.g.
results/pile_thieken_onshore/soil_thieken_onshore/Thieken. With more than
one combination, RESULTS_DIR/summary.csv lists the results of all of them.
"""

import itertools
import math
import os
import re
import shutil
import sys
import time
import traceback

# Abaqus does not always add the script folder to the import path
try:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    BASE_DIR = os.getcwd()
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from abaqus import mdb

import importlib.util

DEFAULT_INPUT_DIR = os.path.join(BASE_DIR, "Inputs")
INPUT_DIR = os.path.abspath(os.path.join(
    BASE_DIR, os.environ.get("MONOPILE_INPUT_DIR") or DEFAULT_INPUT_DIR))


def load_settings(name):
    """Load <name>.py from the input folder (AdvancedSettings: from Inputs
    if the input folder has none). Loaded anew on every run, as an open CAE
    session keeps imported modules."""
    path = os.path.join(INPUT_DIR, name + ".py")
    if not os.path.isfile(path) and name == "AdvancedSettings":
        path = os.path.join(DEFAULT_INPUT_DIR, name + ".py")
    if not os.path.isfile(path):
        raise IOError("%s.py not found in the input folder %s"
                      % (name, INPUT_DIR))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[name] = module
    return module


Input = load_settings("Input")
AdvancedSettings = load_settings("AdvancedSettings")
from Input import (COMMON_COLUMNS, MODEL_COLUMNS, MODEL_NAME, ELEMENT_COUNT,
                   ELEMENT_TYPE, MUDLINE_DEPTH, LOADING, MAX_ITERATIONS,
                   LOAD_CONTROL, HORIZONTAL_LOAD, TARGET_DISPLACEMENT,
                   LOAD_ANGLE, LEVER_ARM, LOAD_INCREMENT, MAX_INCREMENTS,
                   MAX_DISPLACEMENT, MAX_ROTATION, VERTICAL_TANGENT_RATIO,
                   NUM_CPUS, PLOTS)
from AdvancedSettings import (DIAMETER_REFERENCE, SHEAR_FACTOR,
                              ITERATION_TOLERANCE, RELAXATION, CURVE_POINTS,
                              GRID_START)
from methods.csv_reader import (read_cpt_csv, read_pile_csv, read_soil_csv,
                                unused_soil_columns)
from methods.geometry import PileGeometryBuilder
from methods.spring_models import SPRING_MODELS, create_spring_model
from methods.analysis import add_boundary_conditions, add_loading_step
from methods.solver import (DISPLACEMENT, FORCE, BreakCriteria,
                            IncrementalSolver, IterationSettings,
                            load_factors)
from methods.progress import ProgressLog, duration
from methods.postprocessing import (write_pile_head_csv, read_pile_profile,
                                    write_pile_profile_csv, profile_extremes,
                                    write_summary_csv)

RESULTS_DIR = os.path.join(BASE_DIR, Input.RESULTS_DIR)
SUMMARY_CSV = os.path.join(RESULTS_DIR, "summary.csv")
PROGRESS_LOG = os.path.join(RESULTS_DIR, "progress.log")
SUMMARY_PLOT = os.path.join(RESULTS_DIR, "load_displacement_all.png")
# TOTT_BUSWELL settings (optional in Input.py)
CPT_CSV = getattr(Input, "CPT_CSV", None)
TOTT_BUSWELL_ALPHA = getattr(Input, "TOTT_BUSWELL_ALPHA", 0.1)
TOTT_BUSWELL_COMPONENTS = getattr(Input, "TOTT_BUSWELL_COMPONENTS",
                                  ("p-y", "m-theta", "base"))
CPT_MODELS = ("TOTT_BUSWELL",)   # spring models that need a CPT file
JOB_NAME_LENGTH = 80        # maximum length of an Abaqus/CAE job name
MAX_PATH_LENGTH = 259       # Windows MAX_PATH (260) without terminator


class OutputFiles(object):
    """Output files of one calculation in a results folder."""

    def __init__(self, folder, job_name):
        self.folder = folder
        self.job_name = job_name
        self.odb = os.path.join(folder, job_name + ".odb")
        self.cae = os.path.join(folder, job_name + ".cae")
        self.head_csv = os.path.join(folder, "pile_head.csv")
        self.profile_csv = os.path.join(folder, "pile_profile.csv")
        self.load_plot = os.path.join(folder, "load_displacement.png")
        self.profile_plot = os.path.join(folder, "pile_profile.png")
        self.work_dir = os.path.join(folder, "_work")  # temporary job files

    def prepare(self):
        """Check the path length, remove the output of a previous run and
        the temporary files."""
        # Abaqus writes <job>.<ext> (extensions up to ~10 characters) to the
        # work directory; Windows paths are limited to 260 characters
        longest = os.path.join(self.work_dir, self.job_name + ".SMABulk_x")
        if len(longest) > MAX_PATH_LENGTH:
            raise ValueError(
                "Path of the Abaqus job files too long (%d characters, "
                "Windows allows %d): %s. Use shorter file names for pile.csv "
                "and soil.csv, a shorter RESULTS_DIR or move the project to "
                "a shorter folder."
                % (len(longest), MAX_PATH_LENGTH, longest))
        if not os.path.isdir(self.folder):
            os.makedirs(self.folder)
        shutil.rmtree(self.work_dir, ignore_errors=True)
        for path in (self.odb, self.cae, self.head_csv, self.profile_csv,
                     self.load_plot, self.profile_plot,
                     os.path.splitext(self.cae)[0] + ".jnl"):
            try:
                if os.path.exists(path):
                    os.remove(path)
            except OSError:
                # Still open in this CAE session (repeated run); overwritten
                pass

    def clean_up(self):
        """Delete the temporary job files and the CAE journal (.jnl)."""
        shutil.rmtree(self.work_dir, ignore_errors=True)
        journal = os.path.splitext(self.cae)[0] + ".jnl"
        if os.path.exists(journal):
            os.remove(journal)


class Case(object):
    """One combination of pile.csv, soil.csv and spring model."""

    def __init__(self, name, pile_csv, soil_csv, spring_model, folder,
                 cpt_csv=None):
        self.name = name
        self.pile_csv = pile_csv
        self.soil_csv = soil_csv
        self.spring_model = spring_model
        self.cpt_csv = cpt_csv              # Tott-Buswell model only
        self.files = OutputFiles(folder, job_name(name))
        self.results = []                   # LoadLevelResult after the run


def job_name(case_name):
    """Abaqus job name <pile>_<soil>_<model> from the case name
    <pile>/<soil>/<model>: Abaqus/CAE accepts letters, digits and _ only
    (no spaces, dots) and at most 80 characters."""
    return re.sub(r"[^A-Za-z0-9_]", "_", case_name)[:JOB_NAME_LENGTH]


def as_list(value):
    return list(value) if isinstance(value, (list, tuple)) else [value]


def model_label(model):
    """Folder name of a spring model (e.g. THIEKEN -> Thieken)."""
    cls = SPRING_MODELS.get(model)
    return cls.label if cls is not None else model


def cpt_file(soil_csv):
    """CPT file of a soil file: CPT_CSV is one file for all soil files or
    a dict {soil file name: CPT file}. Paths relative to the input folder."""
    if isinstance(CPT_CSV, dict):
        name = os.path.basename(soil_csv)
        if name not in CPT_CSV:
            raise ValueError("CPT_CSV in Input.py has no entry for %s." % name)
        path = CPT_CSV[name]
    else:
        path = CPT_CSV
    if path is None:
        raise ValueError("CPT_CSV in Input.py is needed for the CPT spring "
                         "model.")
    return os.path.normpath(os.path.join(INPUT_DIR, path))


def build_cases():
    """All combinations of PILE_CSV x SOIL_CSV x SPRING_MODEL, each with its
    folder RESULTS_DIR/<pile>/<soil>/<model>."""
    # CSV paths relative to the input folder (absolute paths stay as they are)
    piles = [os.path.normpath(os.path.join(INPUT_DIR, p))
             for p in as_list(Input.PILE_CSV)]
    soils = [os.path.normpath(os.path.join(INPUT_DIR, s))
             for s in as_list(Input.SOIL_CSV)]
    models = [m.upper() for m in as_list(Input.SPRING_MODEL)]

    def stem(path):
        return os.path.splitext(os.path.basename(path))[0]

    cases, names = [], set()
    for pile, soil, model in itertools.product(piles, soils, models):
        base = "/".join((stem(pile), stem(soil), model_label(model)))
        name, i = base, 2
        while name in names:                # same file name in two folders
            name, i = "%s_%d" % (base, i), i + 1
        names.add(name)
        cases.append(Case(name, pile, soil, model,
                          os.path.join(RESULTS_DIR, *name.split("/")),
                          cpt_file(soil) if model in CPT_MODELS else None))
    return cases


def remove_empty_folders(folder):
    """Remove folder and its parents up to RESULTS_DIR while they are
    empty (failed case)."""
    while os.path.normcase(folder) != os.path.normcase(RESULTS_DIR):
        try:
            os.rmdir(folder)                # only if it is empty
        except OSError:
            return
        folder = os.path.dirname(folder)


def main():
    progress = ProgressLog(PROGRESS_LOG)
    try:
        run_all(progress)
    except Exception as err:
        progress.write("ABORTED after %s: %s" % (progress.elapsed(), err))
        raise


def make_plot(plot, report=None):
    """Run plot(plots module); a plotting error is only reported, it does
    not fail the calculation."""
    try:
        from methods import plots
        plot(plots)
    except Exception as err:
        print("WARNING: plot not created: %s" % err)
        if report:
            report("WARNING: plot not created: %s" % err)


def check_settings():
    """Check the load settings of Input.py before the first calculation."""
    if LOAD_CONTROL not in (FORCE, DISPLACEMENT):
        raise ValueError("LOAD_CONTROL must be '%s' or '%s', got %r."
                         % (FORCE, DISPLACEMENT, LOAD_CONTROL))
    if not HORIZONTAL_LOAD > 0.0:
        raise ValueError("HORIZONTAL_LOAD must be positive.")
    if LOAD_CONTROL == DISPLACEMENT and not TARGET_DISPLACEMENT > 0.0:
        raise ValueError("TARGET_DISPLACEMENT must be positive.")
    # Limits of the break criteria are checked by BreakCriteria
    BreakCriteria(MAX_DISPLACEMENT, MAX_ROTATION, VERTICAL_TANGENT_RATIO)


def run_all(progress):
    check_settings()
    cases = build_cases()
    missing = sorted(set(path for case in cases
                         for path in (case.pile_csv, case.soil_csv,
                                      case.cpt_csv)
                         if path is not None and not os.path.isfile(path)))
    if missing:
        raise IOError("Input file(s) not found (CSV paths are relative to "
                      "the input folder %s): %s"
                      % (INPUT_DIR, ", ".join(missing)))
    progress.write("Start: %d calculation(s) (pile x soil x spring model), "
                   "input folder %s" % (len(cases),
                                        os.path.relpath(INPUT_DIR, BASE_DIR)))
    if len(cases) > 1:
        print("%d combinations of pile.csv x soil.csv x spring model."
              % len(cases))
        for path in (SUMMARY_CSV, SUMMARY_PLOT):
            if os.path.exists(path):
                os.remove(path)

    summary = []
    for i, case in enumerate(cases, start=1):
        prefix = "[%d/%d] %s:" % (i, len(cases), case.name)

        def report(message, prefix=prefix):
            progress.write("%s %s" % (prefix, message))

        if len(cases) > 1:
            print("\n" + "=" * 79)
            print("Case %d/%d: %s" % (i, len(cases), case.name))
            print("=" * 79)
        case_start = time.time()
        report("started")
        try:
            summary.append(run_case(case, report))
            report("FINISHED (%s) after %s, results in %s"
                   % (summary[-1]["status"],
                      duration(time.time() - case_start),
                      os.path.relpath(case.files.folder, BASE_DIR)))
        except Exception as err:
            traceback.print_exc()
            print("Case %s FAILED: %s" % (case.name, err))
            report("FAILED after %s: %s"
                   % (duration(time.time() - case_start), err))
            summary.append(case_summary(case, "failed: %s" % err))
            case.files.clean_up()
            remove_empty_folders(case.files.folder)
        if len(cases) > 1:
            write_summary_csv(summary, SUMMARY_CSV)   # updated after each case

    if len(cases) > 1:
        print("\nSummary of all combinations: %s" % SUMMARY_CSV)
        for row in summary:
            print("  %-50s %s" % (row["case"], row["status"]))
        if PLOTS:
            make_plot(lambda plots: plots.plot_load_displacement_all(
                [(case.name, case.results) for case in cases], SUMMARY_PLOT),
                progress.write)
    if len(cases) == 1:
        print("Results in %s" % cases[0].files.folder)
    else:
        print("Results in %s" % RESULTS_DIR)
    # len() instead of sum(): Abaqus/CAE replaces sum() in the namespace
    # of this script by a version that does not accept generators
    ok = len([row for row in summary if row["status"].startswith("ok")])
    progress.write("All calculations finished after %s: %d ok, %d with "
                   "problems%s" % (progress.elapsed(), ok, len(summary) - ok,
                                   ", see summary.csv" if len(cases) > 1
                                   else ""))


def run_case(case, report=None):
    """Build, solve and postprocess one combination. Returns its summary
    row."""
    files = case.files
    files.prepare()

    sections = read_pile_csv(case.pile_csv)
    print("Read %d pile section(s) from %s" % (len(sections), case.pile_csv))
    for s in sections:
        print("  %r" % s)

    # Column numbers of soil.csv for the spring model (Input.py)
    soil_columns = model_soil_columns(case.spring_model)
    soil_layers = read_soil_csv(case.soil_csv, soil_columns)
    print("Read %d soil layer(s) from %s (columns of %s in Input.py)"
          % (len(soil_layers), case.soil_csv, case.spring_model))
    unused = unused_soil_columns(case.soil_csv, soil_columns)
    if unused:
        print("  WARNING: soil.csv column(s) not used by %s, ignored: %s"
              % (case.spring_model, ", ".join(unused)))
    for layer in soil_layers:
        print("  %r" % layer)

    # Nodes at the mudline and at all soil layer boundaries (depths below
    # the pile head), in addition to the pile section boundaries
    layer_depths = [MUDLINE_DEPTH] + [MUDLINE_DEPTH + z for layer
                                      in soil_layers
                                      for z in (layer.z_top, layer.z_bot)]
    builder = PileGeometryBuilder(sections, model_name=MODEL_NAME,
                                  element_count=ELEMENT_COUNT,
                                  element_type=ELEMENT_TYPE,
                                  diameter_reference=DIAMETER_REFERENCE,
                                  shear_factor=SHEAR_FACTOR,
                                  fixed_depths=layer_depths)
    builder.build()
    lengths = [b - a for a, b in zip(builder.node_depths[:-1],
                                     builder.node_depths[1:])]
    print("Pile meshed with %d elements, length %.3f-%.3f m (D = %s "
          "diameter, shear factor %s)."
          % (len(builder.part.elements), min(lengths), max(lengths),
             DIAMETER_REFERENCE,
             "Abaqus default" if SHEAR_FACTOR is None else SHEAR_FACTOR))
    if len(lengths) > ELEMENT_COUNT:
        print("  WARNING: more pile/soil boundaries than ELEMENT_COUNT, %d "
              "elements used." % len(lengths))

    model, instance = builder.model, builder.instance
    add_boundary_conditions(model, instance)
    add_loading_step(model)
    factors = load_factors(LOAD_INCREMENT, MAX_INCREMENTS)
    if LOAD_CONTROL == DISPLACEMENT:
        levels = ("mudline displacement up to %.4g m"
                  % (factors[-1] * TARGET_DISPLACEMENT))
    else:
        levels = "H up to %.4g MN" % (factors[-1] * HORIZONTAL_LOAD / 1e6)
    print("Boundary conditions and step 'Loading' added. %s control, load "
          "at %g deg, lever arm %g m, %d load increment(s): %s."
          % (LOAD_CONTROL, LOAD_ANGLE, LEVER_ARM, len(factors), levels))
    if factors[-1] < 1.0:
        print("WARNING: MAX_INCREMENTS * LOAD_INCREMENT < 1, the full load is "
              "not reached.")
    if report:
        report("model built (%d elements), %d load levels, %s control, %s"
               % (len(lengths), len(factors), LOAD_CONTROL, levels))

    cpt_profile = None
    if case.cpt_csv is not None:
        cpt_profile = read_cpt_csv(case.cpt_csv)
        print("Read %d CPT point(s) from %s (%.2f-%.2f m below mudline%s)"
              % (len(cpt_profile.depths), case.cpt_csv,
                 cpt_profile.depths[0], cpt_profile.z_max,
                 ", with G0" if cpt_profile.has_g0 else ""))

    # Height of the load above the mudline (load eccentricity h)
    spring_model = create_spring_model(
        case.spring_model, sections, soil_layers, MUDLINE_DEPTH,
        loading=LOADING, load_height=LEVER_ARM + MUDLINE_DEPTH,
        diameter_reference=DIAMETER_REFERENCE, curve_points=CURVE_POINTS,
        grid_start=GRID_START, cpt_profile=cpt_profile,
        tott_buswell_alpha=TOTT_BUSWELL_ALPHA,
        tott_buswell_components=TOTT_BUSWELL_COMPONENTS)
    print("Spring model: %s" % spring_model.name)
    for note in spring_model.notes:
        print("  %s" % note)

    solver = IncrementalSolver(
        model, instance, builder.nodes(), spring_model, builder.depth_to_z,
        files.job_name, files.work_dir, files.odb, NUM_CPUS,
        IterationSettings(MAX_ITERATIONS, ITERATION_TOLERANCE, RELAXATION),
        report=report, mudline_depth=MUDLINE_DEPTH,
        criteria=BreakCriteria(MAX_DISPLACEMENT, MAX_ROTATION,
                               VERTICAL_TANGENT_RATIO))
    if LOAD_CONTROL == DISPLACEMENT:
        results = solver.run_displacement(TARGET_DISPLACEMENT, LOAD_ANGLE,
                                          LEVER_ARM, factors, HORIZONTAL_LOAD)
    else:
        results = solver.run(HORIZONTAL_LOAD, LOAD_ANGLE, LEVER_ARM, factors)

    # Postprocessing
    write_pile_head_csv(results, files.head_csv)
    print("Pile head displacements of all increments: %s" % files.head_csv)
    if not results:
        print("WARNING: no load level converged, no .odb/.cae written.")
        files.clean_up()
        return case_summary(case, "failed: no load level converged")

    last = results[-1]
    mdb.saveAs(pathName=files.cae)
    profile = read_pile_profile(files.odb, instance.name, LOAD_ANGLE,
                                MUDLINE_DEPTH)
    write_pile_profile_csv(profile, files.profile_csv, solver.spring_set,
                           spring_model.output_columns)
    print("Displacements and section forces of the last increment: %s"
          % files.profile_csv)
    print("Last load level %d: H = %g N (%s)"
          % (last.increment, last.load,
             "converged" if last.converged else "NOT converged"))
    print("Pile head: U = (%.4f, %.4f, %.4f) m, "
          "UR = (%.5f, %.5f, %.5f) rad" % (last.head_u + last.head_ur))
    m_max, q_max = profile_extremes(profile)
    print("Max. |M| = %.4g Nm, max. |Q| = %.4g N" % (m_max, q_max))
    # Plausibility of the Abaqus section forces at the pile head: M must
    # equal H * lever arm, Q is H minus the spring force at the head node
    head = profile[0]
    print("Check pile head (Abaqus vs. load): M = %.4g Nm (H * lever arm "
          "= %.4g Nm), Q = %.4g N (H = %.4g N)"
          % (head.m, last.moment, head.q, last.load))
    case.results = results
    if PLOTS:
        make_plot(lambda plots: plots.plot_load_displacement(
            results, files.load_plot, case.name, solver.break_info), report)
        make_plot(lambda plots: plots.plot_pile_profile(
            profile, files.profile_plot, case.name, last.load), report)
    files.clean_up()

    not_converged = [r.increment for r in results if not r.converged]
    if not_converged:
        status = "iteration not converged (level %s)" % ", ".join(
            str(k) for k in not_converged)
    elif solver.break_info is not None:
        status = "ok (break criterion)"
    elif last.load_factor < factors[-1] - 1e-9:
        status = "stopped at %.0f %% of the load" % (100.0 * last.load_factor)
    else:
        status = "ok"
    return case_summary(case, status, last, m_max, q_max, solver.break_info)


def case_summary(case, status, last=None, m_max=None, q_max=None,
                 break_info=None):
    """Summary row of one case (values None if it failed)."""
    def relative(path):
        return os.path.relpath(path, BASE_DIR)

    row = {"case": case.name, "pile_csv": relative(case.pile_csv),
           "soil_csv": relative(case.soil_csv),
           "cpt_csv": relative(case.cpt_csv) if case.cpt_csv else None,
           "spring_model": case.spring_model, "status": status,
           "load_factor": None, "H": None, "v_head": None, "psi_head": None,
           "v_mudline": None, "psi_mudline": None,
           "M_max": m_max, "Q_max": q_max,
           "criterion": None, "H_criterion": None, "v_criterion": None,
           "folder": relative(case.files.folder)}
    if last is not None:
        alpha = math.radians(LOAD_ANGLE)
        (ux, uy, _), (rx, ry, _) = last.head_u, last.head_ur
        row.update(load_factor=last.load_factor, H=last.load,
                   v_head=ux * math.cos(alpha) + uy * math.sin(alpha),
                   psi_head=-rx * math.sin(alpha) + ry * math.cos(alpha),
                   v_mudline=last.mudline[0], psi_mudline=last.mudline[1])
    if break_info is not None:
        row.update(criterion="%s %s" % (break_info.name,
                                        break_info.limit_text),
                   H_criterion=break_info.load,
                   v_criterion=break_info.displacement)
    return row


def model_soil_columns(spring_model):
    """Column numbers of soil.csv for a spring model: COMMON_COLUMNS plus
    the model block of MODEL_COLUMNS in Input.py."""
    key = spring_model.upper()
    if key not in MODEL_COLUMNS:
        raise KeyError("Input.py has no column block for spring model '%s' "
                       "(MODEL_COLUMNS: %s)."
                       % (spring_model, ", ".join(sorted(MODEL_COLUMNS))))
    columns = dict(COMMON_COLUMNS)
    columns.update(MODEL_COLUMNS[key])
    return columns


main()
