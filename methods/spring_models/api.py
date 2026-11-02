"""Lateral soil springs from API p-y curves (API RP 2GEO / RP 2A).

Sand:  p = A * pu * tanh(k * H * y / (A * pu))
       pu = min((C1*H + C2*D) * s'v, C3 * D * s'v)
       A  = max(3.0 - 0.8 * H / D, 0.9) static, 0.9 cyclic
Clay:  Matlock (1970) soft clay, ClayPY in base.py (shared by all models)

s'v is the vertical effective stress from the layered gamma_eff, which
replaces gamma'*H in the original single-layer formulas.

Units inside the curve classes: kN, m, kPa.

Soil parameters (soil.csv):
    all layers  gamma_eff  effective unit weight [kN/m^3], required
                           -> vertical effective stress s'v
    sand        phi        friction angle [deg], required
                           -> C1, C2, C3 and the ultimate resistance pu
                k          initial modulus of subgrade reaction [MN/m^3],
                           optional (default: coarse fit of the API chart
                           for submerged sand from phi, API_SAND_K_TABLE)
    clay        su         undrained shear strength [kPa], required,
                           "top-bottom" for a linear variation
                eps50      strain at 50 % of max. stress [-], optional
                           (default 0.01) -> yc = 2.5 eps50 D
                J          Matlock constant [-], optional (default 0.5)
                           -> pu and transition depth XR
"""

import math

from .base import (SpringModel, curve_contribution, matlock_contribution,
                   missing_parameters)

TOLERANCE = 1e-9

# Initial modulus of subgrade reaction k [MN/m^3] for submerged sand vs. phi.
# Reese et al. (1974) values used in API: loose 5.4, medium 16.3, dense 33.9.
# Interpolated linearly, clamped outside the table. This is a coarse fit of
# the API chart: give k explicitly in soil.csv when you have a better value.
API_SAND_K_TABLE = ((29.0, 5.4), (33.0, 16.3), (36.0, 33.9))


def api_sand_k(phi):
    """Initial modulus of subgrade reaction [kN/m^3] for submerged sand."""
    table = API_SAND_K_TABLE
    if phi <= table[0][0]:
        k = table[0][1]
    elif phi >= table[-1][0]:
        k = table[-1][1]
    else:
        for (p0, k0), (p1, k1) in zip(table[:-1], table[1:]):
            if p0 <= phi <= p1:
                k = k0 + (phi - p0) / (p1 - p0) * (k1 - k0)
                break
    return k * 1000.0


def api_sand_coefficients(phi):
    """API coefficients C1, C2, C3 as a function of friction angle [deg]."""
    phi_r = math.radians(phi)
    alpha = phi_r / 2.0
    beta = math.radians(45.0) + phi_r / 2.0
    k0 = 0.4
    ka = math.tan(math.radians(45.0) - phi_r / 2.0) ** 2

    tan_b = math.tan(beta)
    tan_bp = math.tan(beta - phi_r)
    c1 = (tan_b ** 2 * math.tan(alpha) / tan_bp
          + k0 * (math.tan(phi_r) * math.sin(beta)
                  / (math.cos(alpha) * tan_bp)
                  + tan_b * (math.tan(phi_r) * math.sin(beta)
                             - math.tan(alpha))))
    c2 = tan_b / tan_bp - ka
    c3 = ka * (tan_b ** 8 - 1.0) + k0 * math.tan(phi_r) * tan_b ** 4
    return c1, c2, c3


# --------------------------------------------------------------------------- #
# p-y curves (p in kN/m, y in m)
# --------------------------------------------------------------------------- #
class SandPY(object):
    def __init__(self, layer, z, sigma_v, diameter, loading):
        c1, c2, c3 = api_sand_coefficients(layer.phi)
        self.pu = min((c1 * z + c2 * diameter) * sigma_v,
                      c3 * diameter * sigma_v)
        self.a = 0.9 if loading == "cyclic" else max(3.0 - 0.8 * z / diameter,
                                                     0.9)
        k = layer.k * 1000.0 if layer.k is not None else api_sand_k(layer.phi)
        self.kz = k * z

    @property
    def y_ref(self):
        """Displacement at which the curve is practically fully mobilised."""
        if self.kz <= TOLERANCE:
            return 0.0
        return 3.0 * self.a * self.pu / self.kz

    def p(self, y):
        if self.pu <= TOLERANCE or self.kz <= TOLERANCE:
            return 0.0
        apu = self.a * self.pu
        return apu * math.tanh(self.kz * y / apu)


def api_lateral_contribution(point, loading="static", sand_curve=SandPY):
    """Lateral spring term [N vs m] of a soil point from the API curves.
    sand_curve: curve class for sand (SandPY or a variant with another
    initial stiffness, e.g. kallehave.py); clay uses Matlock (base.py)."""
    if point.layer.soil_type != "sand":
        return matlock_contribution(point, loading)
    curve = sand_curve(point.layer, point.z_soil, point.sigma_v,
                       point.diameter, loading)
    return curve_contribution(curve, point.length,
                              max(0.1 * point.diameter, curve.y_ref))


# --------------------------------------------------------------------------- #
# Spring model
# --------------------------------------------------------------------------- #
class ApiPyModel(SpringModel):
    """Lateral p-y springs only, no iteration. Variants with another sand
    curve (kallehave.py, sorensen.py) subclass it and set sand_curve."""

    name = "API"
    label = "API"
    sand_curve = SandPY

    def __init__(self, pile_sections, soil_layers, mudline_depth,
                 loading="static", **_):
        super(ApiPyModel, self).__init__(pile_sections, soil_layers,
                                         mudline_depth)
        if loading not in ("static", "cyclic"):
            raise ValueError("loading must be 'static' or 'cyclic'.")
        self.loading = loading
        errors = missing_parameters(soil_layers, ("phi",))
        if errors:
            raise ValueError("%s: soil parameters in soil.csv (%s)."
                             % (self.name, "; ".join(errors)))

    def distributed_springs(self, point, node_result):
        return api_lateral_contribution(point, self.loading,
                                        self.sand_curve), None
