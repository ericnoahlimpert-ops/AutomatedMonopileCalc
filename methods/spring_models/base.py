"""Common framework for soil spring models.

A spring model describes the soil reaction at a point along the pile
(distributed lateral load p(v) and, optionally, distributed moment m(psi))
and at the pile base (base shear H_B(v) and base moment M_B(psi)). This
module turns those into discrete springs at the pile nodes:

- every node gets the half element above and below it as tributary length
  (clipped at the mudline), each half evaluated at its own mid-depth,
- the contributions are summed and sampled into force-displacement tables.

To add a new model, subclass SpringModel and implement distributed_springs()
and, if needed, base_springs(). Set iterative = True if the curves depend on
the pile response (the analysis is then repeated until convergence, see
main.py) and use node_result = (v, psi) from the previous run. v and psi
are signed: projected onto the direction of the largest pile deflection.
Models that need the whole deflection line override prepare().
Document the soil.csv parameters a model needs in a "Soil parameters
(soil.csv)" section of its module docstring.

Units: soil parameters in kN, m, kPa. All curves handed out by this module
are in N / Nm versus m / rad.
"""

import csv
import math

TOLERANCE = 1e-9
# Default curve sampling (CURVE_POINTS, GRID_START in AdvancedSettings.py)
N_CURVE_POINTS = 30
GRID_START = 1e-4  # first non-zero grid point as fraction of x_max


# --------------------------------------------------------------------------- #
# Curves
# --------------------------------------------------------------------------- #
class Contribution(object):
    """One term of a node spring: force(x) [N or Nm] with x [m or rad].

    x_max       displacement up to which the curve should be sampled
    breakpoints displacements that must be part of the sampling grid
                (e.g. kinks of bilinear curves)
    """

    def __init__(self, function, x_max, breakpoints=()):
        self.function = function
        self.x_max = x_max
        self.breakpoints = tuple(breakpoints)


def curve_contribution(curve, length, x_max, breakpoints=()):
    """Lateral spring term [N vs m] of a p-y curve object (curve.p(y) in
    kN/m, y in m) acting over the pile length [m]."""
    return Contribution(lambda y: 1000.0 * length * curve.p(y), x_max,
                        breakpoints)


class ReactionCurve(object):
    """Positive branch of a spring curve, x >= 0 and force >= 0."""

    def __init__(self, x, force):
        self.x = list(x)
        self.force = list(force)

    @classmethod
    def from_contributions(cls, contributions, n=N_CURVE_POINTS,
                           start=GRID_START):
        contributions = [c for c in contributions if c is not None]
        if not contributions:
            return None
        x_max = max(c.x_max for c in contributions)
        grid = set(geometric_grid(x_max, n, start))
        for c in contributions:
            grid.update(b for b in c.breakpoints if 0.0 < b < x_max)
        grid = sorted(grid)
        force = [sum(c.function(x) for c in contributions) for x in grid]
        curve = cls(grid, force)
        return None if curve.is_zero else curve

    @property
    def is_zero(self):
        return max(self.force) <= TOLERANCE

    def symmetric_table(self):
        """((force, x), ...) sorted by x and mirrored about the origin, as
        needed for nonlinear connector elasticity."""
        positive = [(f, x) for x, f in zip(self.x, self.force)
                    if x > TOLERANCE]
        negative = [(-f, -x) for f, x in reversed(positive)]
        return tuple(negative + [(0.0, 0.0)] + positive)


def geometric_grid(x_max, n=N_CURVE_POINTS, start=GRID_START):
    """n points: 0 plus geometrically spaced points from start * x_max to
    x_max."""
    x_min = start * x_max
    ratio = (x_max / x_min) ** (1.0 / (n - 2))
    return [0.0] + [x_min * ratio ** i for i in range(n - 1)]


# --------------------------------------------------------------------------- #
# Results of a spring model
# --------------------------------------------------------------------------- #
class NodeSprings(object):
    """Springs acting at one pile node. lateral acts on U1 and U2,
    rotational on UR1 and UR2 (both directions use the same curve)."""

    def __init__(self, node_label, z_pile, z_soil, tributary_length,
                 lateral=None, rotational=None):
        self.node_label = node_label
        self.z_pile = z_pile              # depth below pile head [m]
        self.z_soil = z_soil              # depth below mudline [m]
        self.tributary_length = tributary_length
        self.lateral = lateral            # ReactionCurve [N vs m] or None
        self.rotational = rotational      # ReactionCurve [Nm vs rad] or None
        self.output = {}                  # SpringModel.node_output()

    @property
    def is_empty(self):
        return self.lateral is None and self.rotational is None


class SpringSet(object):
    """All springs of one spring model evaluation."""

    def __init__(self, nodes, base=None):
        self.nodes = nodes                # list of NodeSprings
        self.base = base                  # NodeSprings at pile tip or None


class SoilPoint(object):
    """Soil and pile state at one integration point along the pile."""

    def __init__(self, z_pile, z_soil, layer, sigma_v, diameter, length,
                 thickness=None):
        self.thickness = thickness  # pile wall thickness [m]
        self.z_pile = z_pile      # depth below pile head [m]
        self.z_soil = z_soil      # depth below mudline [m]
        self.layer = layer        # SoilLayer
        self.sigma_v = sigma_v    # vertical effective stress [kPa]
        self.diameter = diameter  # pile diameter [m]
        self.length = length      # pile length represented by this point [m]


# --------------------------------------------------------------------------- #
# Base class
# --------------------------------------------------------------------------- #
class SpringModel(object):
    """Base class for soil spring models."""

    name = "base"
    label = "base"               # name of the results folder
    iterative = False            # curves depend on the pile response
    has_rotational_springs = False
    curve_points = N_CURVE_POINTS  # sampling of the curves, see
    grid_start = GRID_START        # geometric_grid()
    output_columns = ()          # extra pile_profile.csv columns, see
                                 # node_output()

    def __init__(self, pile_sections, soil_layers, mudline_depth):
        self.pile_sections = pile_sections
        self.soil_layers = soil_layers
        self.mudline_depth = float(mudline_depth)
        self.warnings = []
        self.notes = []       # information printed once (e.g. G0 method)

        pile_top = pile_sections[0].z_top
        self.pile_tip = pile_sections[-1].z_bot
        if not pile_top <= self.mudline_depth < self.pile_tip:
            raise ValueError("Mudline depth %g lies outside the pile (%g to "
                             "%g)." % (mudline_depth, pile_top, self.pile_tip))
        self.embedded_length = self.pile_tip - self.mudline_depth
        if soil_layers[-1].z_bot < self.embedded_length - 1e-6:
            raise ValueError("Soil profile ends at %g m below mudline, but the "
                             "pile reaches %g m."
                             % (soil_layers[-1].z_bot, self.embedded_length))

    # -- hooks for subclasses ----------------------------------------------- #
    def distributed_springs(self, point, node_result):
        """Return (lateral, rotational) Contribution (or None) for a soil
        point. node_result is (v, psi) of the node from the previous
        analysis, or None in the first run."""
        raise NotImplementedError

    def base_springs(self, point, node_result):
        """Return (shear, moment) Contribution (or None) at the pile base."""
        return None, None

    def prepare(self, nodes, node_results):
        """Called before every evaluation with the sorted node list
        [(label, depth below pile head)] and the results of the previous
        run (or None)."""
        pass

    def node_output(self, z_soil, diameter):
        """Model specific values at a spring node (depth z_soil below
        mudline) for pile_profile.csv: {column in output_columns: value}.
        Evaluated with the springs, so they belong to the same state."""
        return {}

    # -- framework ---------------------------------------------------------- #
    def warn(self, message):
        if message not in self.warnings:
            self.warnings.append(message)

    def clay_fallback(self, point, loading="static"):
        """API Matlock p-y springs (lateral only) for clay layers in models
        that are formulated for sand only, with a warning."""
        self.warn("%s: clay layers use API Matlock p-y curves, lateral only "
                  "(%s is formulated for sand)." % (self.name, self.name))
        return matlock_contribution(point, loading), None

    def compute(self, nodes, node_results=None):
        """Evaluate the springs.

        nodes         list of (node label, depth below pile head)
        node_results  {node label: (v, psi)} from the previous analysis
                      (signed lateral displacement [m] and rotation [rad]
                      in the main deflection direction) or None
        """
        nodes = sorted(nodes, key=lambda n: n[1])
        self.prepare(nodes, node_results)
        node_results = node_results or {}
        springs = []
        for i, (label, z) in enumerate(nodes):
            z_above = 0.5 * (nodes[i - 1][1] + z) if i > 0 else z
            z_below = 0.5 * (z + nodes[i + 1][1]) if i < len(nodes) - 1 else z

            lateral, rotational, length = [], [], 0.0
            for z0, z1 in ((z_above, z), (z, z_below)):
                z0 = max(z0, self.mudline_depth)  # only below the mudline
                if z1 - z0 <= TOLERANCE:
                    continue
                point = self.soil_point(0.5 * (z0 + z1), z1 - z0)
                lat, rot = self.distributed_springs(point,
                                                    node_results.get(label))
                lateral.append(lat)
                rotational.append(rot)
                length += z1 - z0

            if length > 0.0:
                node = NodeSprings(label, z, z - self.mudline_depth, length,
                                   self._curve(lateral),
                                   self._curve(rotational))
                if self.output_columns:
                    node.output = self.node_output(
                        node.z_soil,
                        section_at(self.pile_sections, z).diameter_at(z))
                if not node.is_empty:
                    springs.append(node)

        tip_label, tip_z = nodes[-1]
        base = None
        # Base reactions use the soil just above the tip (layer above a
        # boundary that coincides with the tip)
        point = self.soil_point(tip_z, 0.0, z_layer=tip_z - 1e-6)
        shear, moment = self.base_springs(point, node_results.get(tip_label))
        if shear is not None or moment is not None:
            base = NodeSprings(tip_label, tip_z, tip_z - self.mudline_depth,
                               0.0, self._curve([shear]),
                               self._curve([moment]))
            if base.is_empty:
                base = None
        return SpringSet(springs, base)

    def _curve(self, contributions):
        return ReactionCurve.from_contributions(
            contributions, self.curve_points, self.grid_start)

    def soil_point(self, z_pile, length, z_layer=None):
        z_soil = z_pile - self.mudline_depth
        z_layer_soil = (z_layer if z_layer is not None else z_pile) \
            - self.mudline_depth
        section = section_at(self.pile_sections, z_pile)
        return SoilPoint(z_pile, z_soil,
                         layer_at(self.soil_layers, z_layer_soil),
                         effective_stress(self.soil_layers, z_soil),
                         section.diameter_at(z_pile), length,
                         section.thickness_at(z_pile))

    def mean_embedded_diameter(self, n=200, minus_thickness=False):
        """Average pile diameter over the embedded length. With
        minus_thickness, D - t is averaged (mid-surface diameter of a pile
        whose D is the outer diameter)."""
        dz = self.embedded_length / n
        depths = [self.mudline_depth + (i + 0.5) * dz for i in range(n)]
        total = 0.0
        for z in depths:
            section = section_at(self.pile_sections, z)
            total += section.diameter_at(z)
            if minus_thickness:
                total -= section.thickness_at(z)
        return total / n


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def layer_at(layers, z):
    """Soil layer containing depth z below mudline (lower layer on a
    boundary), or None if z lies outside the profile."""
    for layer in layers:
        if layer.z_top - TOLERANCE <= z < layer.z_bot - TOLERANCE:
            return layer
    if layers and abs(z - layers[-1].z_bot) <= 1e-6:
        return layers[-1]
    return None


def effective_stress(layers, z):
    """Vertical effective stress s'v [kPa] at depth z below mudline."""
    sigma = 0.0
    for layer in layers:
        if z <= layer.z_top:
            break
        sigma += layer.gamma_eff * (min(z, layer.z_bot) - layer.z_top)
    return sigma


def section_at(sections, z):
    for s in sections:
        if s.z_top - 1e-6 <= z <= s.z_bot + 1e-6:
            return s
    raise ValueError("Depth %g is outside the pile." % z)


def missing_parameters(soil_layers, names, soil_type="sand"):
    """One "soil layer i: a, b missing" message per layer of soil_type in
    which any of the soil.csv parameters names is not given."""
    errors = []
    for i, layer in enumerate(soil_layers, start=1):
        if layer.soil_type != soil_type:
            continue
        missing = [name for name in names if getattr(layer, name) is None]
        if missing:
            errors.append("soil layer %d: %s missing"
                          % (i, ", ".join(missing)))
    return errors


def write_springs_csv(spring_set, path):
    """Write all spring curves to a CSV (one row per curve point)."""
    rows = []
    for node in spring_set.nodes:
        rows.append((node, "lateral", node.lateral))
        rows.append((node, "rotational", node.rotational))
    if spring_set.base is not None:
        rows.append((spring_set.base, "base_shear", spring_set.base.lateral))
        rows.append((spring_set.base, "base_moment",
                     spring_set.base.rotational))

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["node", "z_pile [m]", "z_mudline [m]",
                         "tributary_length [m]", "spring",
                         "x [m or rad]", "F [N or Nm]"])
        for node, kind, curve in rows:
            if curve is None:
                continue
            for x, force in zip(curve.x, curve.force):
                writer.writerow([node.node_label, "%.4f" % node.z_pile,
                                 "%.4f" % node.z_soil,
                                 "%.4f" % node.tributary_length, kind,
                                 "%.6e" % x, "%.6e" % force])


# --------------------------------------------------------------------------- #
# Clay: API Matlock (1970) soft clay p-y curves (all models)
# --------------------------------------------------------------------------- #
#     pu = min(3*su + s'v + J*su*X/D, 9*su) * D
#     yc = 2.5 * eps50 * D
#     static: p/pu = 0.5 * (y/yc)^(1/3) up to 8*yc, then 1.0
#     cyclic: as static up to 3*yc (p/pu = 0.72); beyond that 0.72 for
#             X >= XR, otherwise linear decrease to 0.72*X/XR at 15*yc
#
# s'v is the vertical effective stress from the layered gamma_eff, which
# replaces gamma'*X in the original single-layer formula. p in kN/m, y in m.

class ClayPY(object):
    def __init__(self, layer, z, sigma_v, diameter, loading):
        su = layer.su_at(z)
        self.pu = min(3.0 * su + sigma_v + layer.j * su * z / diameter,
                      9.0 * su) * diameter
        self.yc = 2.5 * layer.eps50 * diameter
        self.loading = loading
        # Transition depth for constant su: XR = 6D / (gamma' D / su + J)
        self.xr = 6.0 * diameter / (layer.gamma_eff * diameter / su + layer.j)
        self.z = z

    @property
    def y_ref(self):
        return (15.0 if self.loading == "cyclic" else 8.0) * self.yc

    def p(self, y):
        ratio = y / self.yc
        if self.loading != "cyclic":
            return self.pu * min(0.5 * ratio ** (1.0 / 3.0), 1.0)
        if ratio <= 3.0:
            return self.pu * 0.5 * ratio ** (1.0 / 3.0)
        if self.z >= self.xr:
            return self.pu * 0.72
        p_end = 0.72 * self.z / self.xr
        if ratio >= 15.0:
            return self.pu * p_end
        return self.pu * (0.72 + (p_end - 0.72) * (ratio - 3.0) / 12.0)


def matlock_contribution(point, loading="static"):
    """Lateral spring term [N vs m] of a clay soil point (Matlock)."""
    curve = ClayPY(point.layer, point.z_soil, point.sigma_v, point.diameter,
                   loading)
    return curve_contribution(curve, point.length,
                              max(0.1 * point.diameter, curve.y_ref))


# --------------------------------------------------------------------------- #
# Small-strain shear modulus G0 (PISA, THIEKEN, TOTT_BUSWELL)
# --------------------------------------------------------------------------- #
# The method is chosen per soil layer from the parameters given in soil.csv,
# in this order of priority:
#
# 1. g0                   G0 directly [MPa], "top-bottom" allowed
# 2. g0_ref + lambda_g0   POWER_LAW, stress-dependent power law of the Thieken
#                         model (Thieken et al. 2015, eq. 10 written with a
#                         reference value):
#                             G0 = G0_ref * (s_m / s_ref)^lambda_G0
#                             s_ref = 100 kPa
# 3. hb_b + e0            HARDIN_BLACK, Hardin & Black (1968) as used in the
#                         PISA sand model (Burd et al. 2020, eq. 7):
#                             G0 = B p_ref / (0.3 + 0.7 e0^2) * sqrt(s_m/p_ref)
#                             B = hb_b (soil specific), e0 = void ratio,
#                             p_ref = 101.3 kPa
#                         Dunkirk sand (PISA): B = 875; e0 = 0.741 / 0.685 /
#                         0.629 / 0.573 for Dr = 45 / 60 / 75 / 90 % (Table 2)
#
# Mean effective stress (both methods and Es of the Thieken model):
#     s_m = s'v (1 + 2 k0)/3,  k0 = 1 - sin phi  (phi given in soil.csv)
#                              k0 = 0.4          (fallback without phi;
#                                                 Dunkirk value of PISA)
#
# All stresses and moduli in kPa.

DIRECT = "DIRECT"
POWER_LAW = "POWER_LAW"
HARDIN_BLACK = "HARDIN_BLACK"

DESCRIPTION = {
    DIRECT: "g0 from soil.csv",
    POWER_LAW: "power law (g0_ref, lambda_g0)",
    HARDIN_BLACK: "Hardin & Black (hb_b, e0)",
}

K0_FALLBACK = 0.4    # at-rest coefficient without phi (Dunkirk, PISA)

# Hardin & Black
HB_P_REF = 101.3     # [kPa]

# Power law (Thieken)
SIGMA_REF = 100.0    # [kPa]


def at_rest_k0(phi):
    """Earth pressure coefficient at rest: 1 - sin phi, or K0_FALLBACK if
    phi [deg] is None."""
    if phi is None:
        return K0_FALLBACK
    return 1.0 - math.sin(math.radians(phi))


def mean_stress(sigma_v, phi):
    """Mean effective stress s_m = s'v (1 + 2 k0)/3 [kPa], k0 from
    at_rest_k0(phi)."""
    return sigma_v * (1.0 + 2.0 * at_rest_k0(phi)) / 3.0


def hardin_black_g0(sigma_v, phi, b, e0):
    """G0 [kPa] after Hardin & Black, eq. (7) of Burd et al. (2020).
    b: soil specific constant B [-], e0: void ratio [-]."""
    p_mean = mean_stress(sigma_v, phi)
    return (b * HB_P_REF / (0.3 + 0.7 * e0 ** 2)
            * math.sqrt(max(p_mean, 0.0) / HB_P_REF))


def power_law_g0(sigma_v, phi, g0_ref, lambda_g0):
    """G0 [kPa] = G0_ref (s_m / 100 kPa)^lambda_G0; g0_ref in kPa."""
    return g0_ref * (mean_stress(sigma_v, phi) / SIGMA_REF) ** lambda_g0


def g0_method(layer):
    """G0 method of a sand layer from its soil.csv parameters.

    Returns (method, problem): method is DIRECT, POWER_LAW, HARDIN_BLACK or
    None; problem describes missing or incomplete parameters (None if the
    method is usable)."""
    if layer.g0_top is not None:
        return DIRECT, None
    for method, name, parameters in (
            (POWER_LAW, "power law", (("g0_ref", layer.g0_ref),
                                      ("lambda_g0", layer.lambda_g0))),
            (HARDIN_BLACK, "Hardin & Black", (("hb_b", layer.hb_b),
                                              ("e0", layer.e0)))):
        if any(value is not None for _, value in parameters):
            missing = [p for p, value in parameters if value is None]
            if missing:
                return method, ("%s for G0 incomplete, missing %s"
                                % (name, ", ".join(missing)))
            return method, None
    return None, "no G0 input: give g0, g0_ref + lambda_g0, or hb_b + e0"


def small_strain_modulus(layer, z_soil, sigma_v):
    """G0 [kPa] of a sand layer at depth z_soil below mudline, with the
    method chosen by g0_method()."""
    method, problem = g0_method(layer)
    if problem:
        raise ValueError(problem)
    if method == DIRECT:
        return layer.g0_at(z_soil)
    if method == POWER_LAW:
        return power_law_g0(sigma_v, layer.phi, 1000.0 * layer.g0_ref,
                            layer.lambda_g0)
    return hardin_black_g0(sigma_v, layer.phi, layer.hb_b, layer.e0)


def check_g0_inputs(model_name, soil_layers, soil_type="sand"):
    """Check the G0 input of all layers of soil_type ("sand" or "clay").
    Returns (errors, notes): one error per layer with missing/incomplete
    input, one note per checked layer naming the G0 method used."""
    errors, notes = [], []
    for i, layer in enumerate(soil_layers, start=1):
        if layer.soil_type != soil_type:
            continue
        method, problem = g0_method(layer)
        if problem:
            errors.append("soil layer %d: %s" % (i, problem))
        elif method == DIRECT:
            notes.append("%s: soil layer %d: G0 from %s"
                         % (model_name, i, DESCRIPTION[method]))
        else:
            k0 = ("1 - sin phi = %.3f" % at_rest_k0(layer.phi)
                  if layer.phi is not None
                  else "%g (fallback, no phi)" % K0_FALLBACK)
            notes.append("%s: soil layer %d: G0 from %s, k0 = %s"
                         % (model_name, i, DESCRIPTION[method], k0))
    return errors, notes
