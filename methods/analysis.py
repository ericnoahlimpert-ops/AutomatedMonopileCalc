"""Boundary conditions, steps, loads and job for the monopile model.

Step sequence:
    Initial  (Abaqus default) boundary conditions at the pile tip
    Loading  static step with a horizontal force at the pile head and the
             moment from its lever arm above the head

Units: N, m.

Must be run inside Abaqus/CAE.
"""

import math
import os

from abaqus import mdb
from abaqusConstants import (OFF, ON, ANALYSIS, PERCENTAGE, SINGLE, DEFAULT,
                             LAST_INCREMENT)
import step  # noqa: F401  (registers step methods on the model)
import load  # noqa: F401  (registers load/BC methods on the model)
import job  # noqa: F401
from odbAccess import openOdb

LOADING_STEP = "Loading"
# RF and SF include reaction moments (RM) and section moments (SM)
FIELD_OUTPUT = ("U", "UR", "RF", "SF", "CTF", "CU")


def add_boundary_conditions(model, instance):
    """Pile tip: no vertical displacement (U3) and no twist about the pile
    axis (UR3). Lateral support comes from the soil springs."""
    model.DisplacementBC(name="BC-Tip", createStepName="Initial",
                         region=instance.sets["PILE_TIP"], u3=0.0, ur3=0.0)


def add_loading_step(model, name=LOADING_STEP, initial_inc=0.1, min_inc=1e-6,
                     max_num_inc=1000, nlgeom=False):
    """Static step for the nonlinear soil springs.

    initial_inc = 1.0 (full load in one increment) was tested: only ~1.5 %
    faster for C4, because the job time is dominated by the Abaqus start-up,
    so the more robust 0.1 is kept."""
    model.StaticStep(name=name, previous="Initial",
                     description="Horizontal load at pile head",
                     timePeriod=1.0, initialInc=initial_inc, minInc=min_inc,
                     maxInc=1.0, maxNumInc=max_num_inc,
                     nlgeom=ON if nlgeom else OFF)
    # Only the end of the step is evaluated: write output at the last
    # increment only (smaller odb, less writing)
    model.fieldOutputRequests["F-Output-1"].setValues(
        variables=FIELD_OUTPUT, frequency=LAST_INCREMENT)
    for request in model.historyOutputRequests.values():
        request.setValues(frequency=LAST_INCREMENT)


def add_horizontal_load(model, instance, magnitude, angle=0.0,
                        lever_arm=0.0, step_name=LOADING_STEP):
    """Horizontal load H [N] acting lever_arm [m] above the pile head, applied
    as force H plus moment M = H * lever_arm at the pile head, so the tower
    does not have to be modelled. angle [deg] is measured from the global X
    axis towards Y. Returns the moment magnitude [Nm].

    Calling it again updates the magnitude of the existing loads (used for
    incremental loading)."""
    angle_r = math.radians(angle)
    fx = magnitude * math.cos(angle_r)
    fy = magnitude * math.sin(angle_r)
    region = instance.sets["PILE_HEAD"]
    if "Load-H" in model.loads:
        model.loads["Load-H"].setValues(cf1=fx, cf2=fy)
    else:
        model.ConcentratedForce(name="Load-H", createStepName=step_name,
                                region=region, cf1=fx, cf2=fy)

    if abs(lever_arm) > 0.0:
        # M = r x F with r = (0, 0, lever_arm) pointing up from the head
        cm1, cm2 = -lever_arm * fy, lever_arm * fx
        if "Load-M" in model.loads:
            model.loads["Load-M"].setValues(cm1=cm1, cm2=cm2)
        else:
            model.Moment(name="Load-M", createStepName=step_name,
                         region=region, cm1=cm1, cm2=cm2)
    return abs(magnitude * lever_arm)


def run_job(model, job_name, work_dir, num_cpus=1):
    """Create the job, run it in work_dir and wait for it to finish.
    Returns the path of the output database."""
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    if job_name in mdb.jobs:
        del mdb.jobs[job_name]
    analysis_job = mdb.Job(name=job_name, model=model.name, type=ANALYSIS,
                           numCpus=num_cpus, numDomains=num_cpus,
                           memory=90, memoryUnits=PERCENTAGE,
                           explicitPrecision=SINGLE,
                           nodalOutputPrecision=SINGLE,
                           multiprocessingMode=DEFAULT)

    old_dir = os.getcwd()
    os.chdir(work_dir)  # Abaqus writes the job files to the cwd
    try:
        analysis_job.submit(consistencyChecking=OFF)
        analysis_job.waitForCompletion()
    finally:
        os.chdir(old_dir)

    sta_file = os.path.join(work_dir, job_name + ".sta")
    if not _completed(sta_file):
        raise RuntimeError("Job '%s' did not complete successfully, see the "
                           ".msg/.dat files in %s." % (job_name, work_dir))
    return os.path.join(work_dir, job_name + ".odb")


def _completed(sta_file):
    if not os.path.isfile(sta_file):
        return False
    with open(sta_file) as f:
        return "COMPLETED SUCCESSFULLY" in f.read()


def read_head_results(odb_path, instance_name, step_name=LOADING_STEP):
    """Displacement and rotation of the pile head at the end of the step."""
    odb = openOdb(odb_path, readOnly=True)
    try:
        frame = odb.steps[step_name].frames[-1]
        head = odb.rootAssembly.instances[instance_name.upper()] \
            .nodeSets["PILE_HEAD"]
        u = frame.fieldOutputs["U"].getSubset(region=head).values[0].data
        ur = frame.fieldOutputs["UR"].getSubset(region=head).values[0].data
        return tuple(float(v) for v in u), tuple(float(v) for v in ur)
    finally:
        odb.close()


def read_head_results_all_frames(odb_path, instance_name,
                                 step_name=LOADING_STEP):
    """Pile head displacement and rotation of every output frame after the
    start of the step: [(step time, (U1, U2, U3), (UR1, UR2, UR3))]."""
    odb = openOdb(odb_path, readOnly=True)
    try:
        head = odb.rootAssembly.instances[instance_name.upper()] \
            .nodeSets["PILE_HEAD"]
        results = []
        for frame in odb.steps[step_name].frames:
            if frame.frameValue <= 0.0:
                continue
            u = frame.fieldOutputs["U"].getSubset(region=head).values[0].data
            ur = frame.fieldOutputs["UR"].getSubset(region=head).values[0].data
            results.append((float(frame.frameValue),
                            tuple(float(v) for v in u),
                            tuple(float(v) for v in ur)))
        return results
    finally:
        odb.close()


def set_output_time_points(model, times, name="LoadLevels",
                           step_name=LOADING_STEP):
    """Write the field output at the given step times (e.g. one frame per
    load level) instead of only at the end of the step. Replaces time
    points of the same name from an earlier call."""
    if name in model.timePoints:
        del model.timePoints[name]
    model.TimePoint(name=name, points=tuple((t,) for t in times))
    model.fieldOutputRequests["F-Output-1"].setValues(timePoint=name)


def _frame_node_results(frame, instance):
    """Signed v [m] and psi [rad] of every node of one frame, see
    read_node_results."""
    u = frame.fieldOutputs["U"].getSubset(region=instance).values
    ur = frame.fieldOutputs["UR"].getSubset(region=instance).values
    displacements = dict((val.nodeLabel, (val.data[0], val.data[1]))
                         for val in u)
    rotations = dict((val.nodeLabel, (val.data[0], val.data[1]))
                     for val in ur)
    ux, uy = max(displacements.values(), key=lambda d: math.hypot(*d))
    length = math.hypot(ux, uy) or 1.0
    ex, ey = ux / length, uy / length
    return dict((label, (float(dx * ex + dy * ey),
                         float(-ey * rotations[label][0]
                               + ex * rotations[label][1])))
                for label, (dx, dy) in displacements.items())


def read_node_results(odb_path, instance_name, step_name=LOADING_STEP):
    """Signed lateral displacement v [m] and rotation psi [rad] of every pile
    node at the end of the step: {node label: (v, psi)}.

    Both are projected onto the direction e of the largest horizontal
    displacement: v = U.e, psi = rotation about the horizontal axis normal
    to e (positive when the pile head moves in +e)."""
    odb = openOdb(odb_path, readOnly=True)
    try:
        frame = odb.steps[step_name].frames[-1]
        instance = odb.rootAssembly.instances[instance_name.upper()]
        return _frame_node_results(frame, instance)
    finally:
        odb.close()


def read_node_results_all_frames(odb_path, instance_name,
                                 step_name=LOADING_STEP):
    """read_node_results for every output frame after the start of the
    step: [(step time, {node label: (v, psi)})]."""
    odb = openOdb(odb_path, readOnly=True)
    try:
        instance = odb.rootAssembly.instances[instance_name.upper()]
        return [(float(frame.frameValue),
                 _frame_node_results(frame, instance))
                for frame in odb.steps[step_name].frames
                if frame.frameValue > 0.0]
    finally:
        odb.close()
