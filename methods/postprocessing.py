"""Postprocessing: result tables of the incremental analysis.

pile_head.csv      pile head displacements and rotations of every converged
                   load increment
pile_profile.csv   final load increment, one row per pile node:
                   displacements, rotations and section forces N, Q, M

Section forces come from the Abaqus field outputs SF (section forces) and
SM (section moments) at the integration points of the beam elements; the
mean over the integration points of an element is its value at
mid-element. All values written are Abaqus values (nothing is replaced by
the applied load), so the table can be used to check the Abaqus results:
    inner nodes   N, Q, M linearly interpolated between the mid-element
                  values of the two adjacent elements
    head and tip  only one adjacent element: N and Q of that element (both
                  constant within an element, as the springs act at the
                  nodes only), M = M_mid + Q (z_node - z_mid) along the
                  element (M linear with dM/dz = Q, z downwards). This
                  replaces the linear extrapolation over two elements,
                  which is distorted by the spring force at the node
                  in between.
    Q_above/below shear of the element above/below the node (empty above
                  the head and below the tip); the jump Q_above - Q_below
                  is the spring force at the node.
For a plausibility check, Q at the head (shear just below the head node,
i.e. H minus the spring force at the head node) and M at the head
(= H * lever arm) can be compared with the applied load, see main.py.

Sign conventions (load angle alpha from global X towards Y):
    v    = U1 cos(alpha) + U2 sin(alpha)      displacement in load direction
    psi  = -UR1 sin(alpha) + UR2 cos(alpha)   rotation in the load plane,
                                              positive when the head moves
                                              in load direction
    N    = SF1                                axial force (tension positive)
    Q    = SF2 sin(alpha) - SF3 cos(alpha)    shear force in the load plane,
                                              Q = +H at the pile head
    M    = SM1 sin(alpha) + SM2 cos(alpha)    bending moment in the load
                                              plane, M = +H * lever arm at
                                              the pile head
The SF/SM mapping follows from the beam orientation n1 = global X used in
geometry.py (checked with loads at 0, 45 and 90 deg).

Must be run inside Abaqus (odbAccess).
"""

import csv
import math

from odbAccess import openOdb

from .analysis import LOADING_STEP


# --------------------------------------------------------------------------- #
# Pile head, every load increment
# --------------------------------------------------------------------------- #
def write_pile_head_csv(results, path):
    """One row per load level (results: list of solver.LoadLevelResult):
    load, pile head displacements and rotations (global components) and the
    displacement v / rotation psi in the main deflection direction at the
    mudline and the pile toe. load_factor is the fraction of HORIZONTAL_LOAD
    (force control) or TARGET_DISPLACEMENT (displacement control)."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["increment", "load_factor", "H [N]", "M_head [Nm]",
                         "U1 [m]", "U2 [m]", "U3 [m]",
                         "UR1 [rad]", "UR2 [rad]", "UR3 [rad]",
                         "v_mudline [m]", "psi_mudline [rad]",
                         "v_toe [m]", "psi_toe [rad]",
                         "iterations", "converged"])
        for r in results:
            writer.writerow([r.increment, "%.4f" % r.load_factor,
                             "%.6e" % r.load, "%.6e" % r.moment]
                            + ["%.6e" % v for v in r.head_u + r.head_ur]
                            + ["" if v is None else "%.6e" % v
                               for v in r.mudline + r.toe]
                            + [r.iterations, r.converged])


# --------------------------------------------------------------------------- #
# Pile profile, final load increment
# --------------------------------------------------------------------------- #
class ProfileRow(object):
    def __init__(self, node, z_pile, z_mudline, u, ur, v, psi, n, q, m,
                 q_above=None, q_below=None):
        self.node = node
        self.z_pile = z_pile          # depth below pile head [m]
        self.z_mudline = z_mudline    # depth below mudline [m]
        self.u = u                    # (U1, U2, U3) [m]
        self.ur = ur                  # (UR1, UR2, UR3) [rad]
        self.v = v                    # displacement in load direction [m]
        self.psi = psi                # rotation in the load plane [rad]
        self.n = n                    # axial force [N]
        self.q = q                    # shear force in the load plane [N]
        self.m = m                    # bending moment in the load plane [Nm]
        self.q_above = q_above        # shear just above the node [N]
        self.q_below = q_below        # shear just below the node [N]


class _ElementForces(object):
    """Constant N, Q and linear M of one beam element in the load plane."""

    def __init__(self, z_top, z_bot, n, q, m_mid):
        self.z_top = z_top
        self.z_bot = z_bot
        self.n = n
        self.q = q
        self.m_mid = m_mid

    def moment_at(self, z):
        return self.m_mid + self.q * (z - 0.5 * (self.z_top + self.z_bot))


def _mean(values):
    return sum(values) / len(values)


def _section_forces(node_depth, above, below):
    """(N, Q, M) at a node from the Abaqus values of the adjacent
    elements, see the module docstring."""
    if above and below and above is not below:
        z0 = 0.5 * (above.z_top + above.z_bot)
        z1 = 0.5 * (below.z_top + below.z_bot)
        w = (node_depth - z0) / (z1 - z0)
        return tuple(a + w * (b - a) for a, b in
                     ((above.n, below.n), (above.q, below.q),
                      (above.m_mid, below.m_mid)))
    # pile head or tip, or a node inside an element (B32 mid-node)
    element = above or below
    return element.n, element.q, element.moment_at(node_depth)


def read_pile_profile(odb_path, instance_name, load_angle, mudline_depth,
                      step_name=LOADING_STEP):
    """Displacements, rotations and section forces of all pile nodes at the
    end of the step, sorted from pile head to tip."""
    alpha = math.radians(load_angle)
    ca, sa = math.cos(alpha), math.sin(alpha)

    odb = openOdb(odb_path, readOnly=True)
    try:
        frame = odb.steps[step_name].frames[-1]
        inst = odb.rootAssembly.instances[instance_name.upper()]
        depth = dict((n.label, 0.0 - n.coordinates[2]) for n in inst.nodes)
        u = dict((val.nodeLabel, tuple(float(x) for x in val.data)) for val in
                 frame.fieldOutputs["U"].getSubset(region=inst).values)
        ur = dict((val.nodeLabel, tuple(float(x) for x in val.data)) for val
                  in frame.fieldOutputs["UR"].getSubset(region=inst).values)
        # all integration points of an element (one for B31, two for B32)
        sf, sm = {}, {}
        for output, values in ((sf, frame.fieldOutputs["SF"]),
                               (sm, frame.fieldOutputs["SM"])):
            for val in values.getSubset(region=inst).values:
                output.setdefault(val.elementLabel, []).append(
                    [float(x) for x in val.data])
        connectivity = dict((e.label, e.connectivity) for e in inst.elements)
    finally:
        odb.close()

    # Element forces in the load plane; the mean over the integration
    # points is the value at mid-element (N, Q constant, M linear)
    elements = []
    node_elements = {}
    for label, values in sf.items():
        ends = sorted(depth[n] for n in (connectivity[label][0],
                                         connectivity[label][-1]))
        element = _ElementForces(
            ends[0], ends[1], _mean([v[0] for v in values]),
            _mean([v[1] * sa - v[2] * ca for v in values]),
            _mean([v[0] * sa + v[1] * ca for v in sm[label]]))
        elements.append(element)
        for node in connectivity[label]:
            node_elements.setdefault(node, []).append(element)

    rows = []
    for node in sorted(u, key=lambda n: depth[n]):
        d = depth[node]
        adjacent = node_elements[node]
        above = next((e for e in adjacent if e.z_top < d - 1e-6), None)
        below = next((e for e in adjacent if e.z_bot > d + 1e-6), None)
        n, q, m = _section_forces(d, above, below)

        ux, uy, _ = u[node]
        rx, ry, _ = ur[node]
        rows.append(ProfileRow(
            node, d, d - mudline_depth, u[node], ur[node],
            ux * ca + uy * sa, -rx * sa + ry * ca, n, q, m,
            above.q if above else None, below.q if below else None))
    return rows


def write_summary_csv(rows, path):
    """One row per calculated combination (batch runs, see main.py)."""
    columns = (("case", "case", None), ("pile_csv", "pile_csv", None),
               ("soil_csv", "soil_csv", None),
               ("cpt_csv", "cpt_csv", None),
               ("spring_model", "spring_model", None),
               ("status", "status", None),
               ("load_factor", "load_factor", "%.4f"),
               ("H [N]", "H", "%.6e"), ("v_head [m]", "v_head", "%.6e"),
               ("psi_head [rad]", "psi_head", "%.6e"),
               ("v_mudline [m]", "v_mudline", "%.6e"),
               ("psi_mudline [rad]", "psi_mudline", "%.6e"),
               ("M_max [Nm]", "M_max", "%.6e"), ("Q_max [N]", "Q_max", "%.6e"),
               ("break criterion", "criterion", None),
               ("H_criterion [N]", "H_criterion", "%.6e"),
               ("v_criterion [m]", "v_criterion", "%.6e"),
               ("folder", "folder", None))
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([title for title, _, _ in columns])
        for row in rows:
            writer.writerow([
                "" if row[key] is None
                else (fmt % row[key] if fmt else row[key])
                for _, key, fmt in columns])


def profile_extremes(rows):
    """(max |M| [Nm], max |Q| [N]) along the pile."""
    return (max(abs(r.m) for r in rows), max(abs(r.q) for r in rows))


def write_pile_profile_csv(rows, path, spring_set=None, columns=()):
    """One row per pile node: depth, displacements, rotations, N, Q, M.

    columns: extra model specific columns (SpringModel.output_columns),
    taken from the node outputs of spring_set (the springs of the written
    load level); empty for nodes without springs (above the mudline)."""
    outputs = {}
    if spring_set is not None:
        outputs = dict((node.node_label, node.output)
                       for node in spring_set.nodes)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["node", "z_pile [m]", "z_mudline [m]",
                         "U1 [m]", "U2 [m]", "U3 [m]",
                         "UR1 [rad]", "UR2 [rad]", "UR3 [rad]",
                         "v [m]", "psi [rad]",
                         "N [N]", "Q [N]", "M [Nm]",
                         "Q_above [N]", "Q_below [N]"] + list(columns))
        for r in rows:
            output = outputs.get(r.node, {})
            writer.writerow([r.node, "%.4f" % r.z_pile, "%.4f" % r.z_mudline]
                            + ["%.6e" % x for x in r.u + r.ur]
                            + ["%.6e" % x for x in (r.v, r.psi, r.n, r.q,
                                                    r.m)]
                            + ["" if x is None else "%.6e" % x
                               for x in (r.q_above, r.q_below)]
                            + ["" if output.get(c) is None
                               else "%.4f" % output[c] for c in columns])
