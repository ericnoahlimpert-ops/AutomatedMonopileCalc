"""Thieken et al. (2015) p-y approach for piles of arbitrary dimensions in sand.

Thieken, K., Achmus, M., Lemke, K. (2015). A new static p-y approach for
piles with arbitrary dimensions in sand. geotechnik 38, No. 4, 267-288.
https://doi.org/10.1002/gete.201400036

Ultimate resistance (DIN 4085, eqs 17-19, 29):
    K_pgh   = ((1+sin phi)/(1-sin phi) * (1 - 0.53 d_p)^(0.26+5.96 phi))
              / sqrt(1 + tan^2 d_p),   d_p = -2/3 phi   (angles in rad)
    pC_basic = 11/16 * gamma' * z^1.5 * K_pgh * (1 + 2 tan phi) * sqrt(D)
    pC       = pC_basic * [3 - 2 (z / 2.5D)^0.25] >= pC_basic

Basic p-y curve with stretching factor SF (eqs 20-28, Fig. 14):
    yB = pC / (Es SF),  yC = 4 yB,  yA = 0.005 D yB,  pB = RB pC
    RB = 0.63 - 0.33 (z/D / 7.5)^0.5 >= 0.35
    y <= yB:       p = pB (y/yB)^(1/n) <= E0py * y
    yB < y < yC:   linear from pB to pC
    y >= yC:       p = pC
    n  = nA                          for y >= yA
       = n0 + (nA - n0)(y/yA)^0.25   for y <  yA
    nA = 2.4 - 0.08 (phi[deg] - 30)
    n0 = log(1e-5 D / yB) / log(E0py 1e-5 D / pB)
    E0py = Esd * SF = G0 * 2(1-nu^2)/(1-nu-2nu^2) * SF

Soil stiffness (power laws of eqs 10, 14 with reference values), s_ref =
100 kPa:
    G0 = G0_ref * (s_m/s_ref)^lambda_G0     (power law; g0 directly or
                                             Hardin & Black (hb_b, e0) also
                                             possible, see base.py)
    Es = Es_ref * (s_m/s_ref)^lambda_ES
    s_m = s'v (1 + 2 k0)/3,  k0 = 1 - sin phi
In the paper's notation Es_ref = kappa * s_ref and
G0_ref = 6900 (2.17-e)^2/(1+e) * s_ref^lambda_G0 [kPa].

Stretching factor SF = SF_bend + SF_tip (eqs 30-35), from the deflection
line of the previous run (z below mudline, L embedded length):
    z <= z_y0:            SF_bend = 2.3 y/y_max + 0.3
    z_y0 < z <= z_ymin:   SF_bend = 1.2 y/y_min + 0.3
    z > z_ymin:           SF_bend = 1.2 y/y_min + 0.3 >= 1.0
    z >= L - 2D:          SF_tip  = 5 ((z-L)/(2D) + 1)^5
    z >= L - 0.1D:        SF_tip  = 5 ((z-L)/(2D) + 1)^5 + 3
The first run uses the basic curves (SF = 1), then the analysis is repeated
until the pile head displacement converges.

Layered soil: gamma' * z is replaced by the vertical effective stress s'v,
i.e. pC_basic = 11/16 * s'v * sqrt(z) * ...; phi and the stiffness
parameters are taken from the local layer. The approach was derived for
homogeneous sand. Clay layers fall back to API Matlock curves.

Note on yA: the text gives yA = 0.005 pC/(Es SF), Fig. 14 shows an
additional factor D. The Fig. 14 version yA = 0.005 D pC/(Es SF) is used
(D in m).

Soil parameters (soil.csv, no defaults):
    all layers  gamma_eff  effective unit weight [kN/m^3], required
                           -> vertical effective stress s'v (pC, s_m)
    sand        phi        friction angle [deg], required
                           -> K_pgh and pC, exponent nA, k0 = 1 - sin phi
                nu         Poisson's ratio [-], required
                           -> Esd = G0 * 2(1-nu^2)/(1-nu-2nu^2), i.e. E0py
                es_ref     oedometric stiffness Es at s_m = 100 kPa [MPa],
                           required -> yB, yC, yA (displacement scale)
                lambda_es  stress exponent of Es [-], required
                G0 -> initial stiffness E0py, chosen from the given
                columns (base.py), priority:
                g0         G0 directly [MPa], "top-bottom" allowed
                g0_ref     G0 at s_m = 100 kPa [MPa]  power law (Thieken
                lambda_g0  stress exponent of G0 [-]  formulation)
                hb_b       Hardin & Black constant B  Hardin & Black as in
                e0         void ratio [-]             PISA
    clay        su, eps50, J as for API Matlock (Thieken is for sand only,
                API p-y curves are used)

For orientation, Table 1 of the paper lists the values used in the
authors' FE study (North Sea sands), converted to Es_ref and G0_ref:

    description   nu     Es_ref  lambda_ES  G0_ref  lambda_G0  phi
                         [MPa]              [MPa]
    very dense    0.200  70.0    0.50       106.3   0.50       40.0
    dense         0.225  50.0    0.55        96.6   0.50       37.5
    medium dense  0.250  40.0    0.60        89.4   0.50       35.0
    loose         0.275  32.5    0.65        81.1   0.50       32.5
    very loose    0.300  25.0    0.70        71.9   0.50       30.0
"""

import math

from .base import (SIGMA_REF, TOLERANCE, SpringModel, check_g0_inputs,
                   curve_contribution, mean_stress, missing_parameters,
                   small_strain_modulus)

# soil.csv columns every sand layer needs (plus those of the G0 method)
REQUIRED_PARAMETERS = ("phi", "nu", "es_ref", "lambda_es")


def k_pgh(phi_deg):
    """Passive earth pressure coefficient after DIN 4085, eq. (18)."""
    phi = math.radians(phi_deg)
    delta = -2.0 / 3.0 * phi
    kp = ((1.0 + math.sin(phi)) / (1.0 - math.sin(phi))
          * (1.0 - 0.53 * delta) ** (0.26 + 5.96 * phi))
    return kp / math.sqrt(1.0 + math.tan(delta) ** 2)


class ThiekenPY(object):
    """p-y curve (p in kN/m, y in m) at one depth with stretching factor."""

    def __init__(self, layer_params, z, sigma_v, diameter, sf=1.0):
        # es_ref and g0 (G0 at this depth) in kPa
        phi, nu, es_ref, lambda_es, g0 = layer_params
        z = max(z, 1e-3)
        sigma_v = max(sigma_v, 1e-3)
        phi_r = math.radians(phi)

        # Ultimate resistance, eqs (19), (29)
        pc_basic = (11.0 / 16.0 * sigma_v * math.sqrt(z) * k_pgh(phi)
                    * (1.0 + 2.0 * math.tan(phi_r)) * math.sqrt(diameter))
        self.pc = pc_basic * max(3.0 - 2.0 * (z / (2.5 * diameter)) ** 0.25,
                                 1.0)

        # Soil stiffness, eqs (10), (11), (14) with reference values
        sigma_m = mean_stress(sigma_v, phi)  # k0 = 1 - sin phi
        es = es_ref * (sigma_m / SIGMA_REF) ** lambda_es
        self.e0py = g0 * 2.0 * (1.0 - nu ** 2) / (1.0 - nu - 2.0 * nu ** 2) \
            * sf

        # Supporting points, eqs (20)-(28)
        self.yb = self.pc / (es * sf)
        self.yc = 4.0 * self.yb
        self.ya = 0.005 * diameter * self.yb
        rb = max(0.63 - 0.33 * (z / diameter / 7.5) ** 0.5, 0.35)
        self.pb = rb * self.pc
        self.n_a = 2.4 - 0.08 * (phi - 30.0)
        self.y_0 = 1e-5 * diameter
        self.n_0 = (math.log(self.y_0 / self.yb)
                    / math.log(self.e0py * self.y_0 / self.pb))

    def p(self, y):
        if y <= 0.0:
            return 0.0
        if y >= self.yc:
            return self.pc
        if y > self.yb:
            return self.pb + (self.pc - self.pb) * (y - self.yb) \
                / (self.yc - self.yb)
        if y >= self.ya:
            n = self.n_a
        else:
            n = self.n_0 + (self.n_a - self.n_0) * (y / self.ya) ** 0.25
        return min(self.pb * (y / self.yb) ** (1.0 / n), self.e0py * y)


class ThiekenSandModel(SpringModel):
    """Iterative Thieken et al. (2015) p-y springs (lateral only)."""

    name = "THIEKEN"
    label = "Thieken"
    iterative = True
    has_rotational_springs = False
    output_columns = ("SF_bend [-]", "SF_tip [-]", "SF [-]")

    def __init__(self, pile_sections, soil_layers, mudline_depth, **_):
        super(ThiekenSandModel, self).__init__(pile_sections, soil_layers,
                                               mudline_depth)
        self.diameter_ref = self.mean_embedded_diameter()
        self.deflection = None
        if len(soil_layers) > 1:
            self.warn("THIEKEN: derived for homogeneous sand; the layered "
                      "profile is evaluated layer by layer.")
        errors, self.notes = check_g0_inputs("THIEKEN", soil_layers)
        errors += missing_parameters(soil_layers, REQUIRED_PARAMETERS)
        if errors:
            raise ValueError(
                "THIEKEN: soil parameters in soil.csv (%s). Every sand layer "
                "needs %s and G0 (g0, g0_ref + lambda_g0, or hb_b + e0)."
                % ("; ".join(errors), ", ".join(REQUIRED_PARAMETERS)))

    # -- soil parameters ---------------------------------------------------- #
    def _layer_params(self, layer, z_soil, sigma_v):
        """(phi, nu, Es_ref [kPa], lambda_ES, G0 [kPa]) at a soil point."""
        g0 = small_strain_modulus(layer, z_soil, max(sigma_v, 1e-3))
        return (layer.phi, layer.nu, 1000.0 * layer.es_ref, layer.lambda_es,
                g0)

    # -- deflection line ---------------------------------------------------- #
    def prepare(self, nodes, node_results):
        """Extract y_max, z_y0, y_min, z_ymin from the deflection line of the
        previous run (depths below mudline)."""
        if not node_results:
            self.deflection = None
            return
        line = [(z - self.mudline_depth, node_results[label][0])
                for label, z in nodes
                if z >= self.mudline_depth - TOLERANCE
                and label in node_results]
        y_max = max(y for _, y in line)
        if y_max <= 0.0:
            self.deflection = None
            return
        z_y0 = None
        for (z0, y0), (z1, y1) in zip(line[:-1], line[1:]):
            if y0 > 0.0 >= y1:
                z_y0 = z0 + (z1 - z0) * y0 / (y0 - y1)
                break
        below = [(z, y) for z, y in line if z_y0 is not None and z > z_y0]
        if below:
            z_ymin, y_min = min(below, key=lambda zy: zy[1])
        else:
            z_ymin, y_min = None, None
        self.deflection = (line, y_max, z_y0, y_min, z_ymin)

    def _deflection_at(self, z):
        line = self.deflection[0]
        if z <= line[0][0]:
            return line[0][1]
        for (z0, y0), (z1, y1) in zip(line[:-1], line[1:]):
            if z0 <= z <= z1:
                return y0 + (y1 - y0) * (z - z0) / (z1 - z0)
        return line[-1][1]

    def stretching_factor(self, z, diameter):
        """SF = SF_bend + SF_tip at depth z below mudline (1.0 in the first
        run)."""
        return max(sum(self.stretching_components(z, diameter)), 1e-3)

    def stretching_components(self, z, diameter):
        """(SF_bend, SF_tip) at depth z below mudline; (1.0, 0.0) in the
        first run (basic curves)."""
        if self.deflection is None:
            return 1.0, 0.0
        _, y_max, z_y0, y_min, z_ymin = self.deflection
        y = self._deflection_at(z)
        if z_y0 is None or z <= z_y0:
            sf_bend = 2.3 * y / y_max + 0.3
        elif y_min is None or y_min >= 0.0:
            sf_bend = 0.3
        elif z <= z_ymin:
            sf_bend = 1.2 * y / y_min + 0.3
        else:
            sf_bend = max(1.2 * y / y_min + 0.3, 1.0)

        length = self.embedded_length
        sf_tip = 0.0
        if z >= length - 2.0 * diameter:
            sf_tip = 5.0 * ((z - length) / (2.0 * diameter) + 1.0) ** 5
            if z >= length - 0.1 * diameter:
                sf_tip += 3.0
        return sf_bend, sf_tip

    def node_output(self, z_soil, diameter):
        """Stretching factors at the node depth. The springs use SF at the
        mid-depths of the half elements above and below the node."""
        sf_bend, sf_tip = self.stretching_components(z_soil, diameter)
        return {"SF_bend [-]": sf_bend, "SF_tip [-]": sf_tip,
                "SF [-]": max(sf_bend + sf_tip, 1e-3)}

    # -- soil reactions ----------------------------------------------------- #
    def distributed_springs(self, point, node_result):
        if point.layer.soil_type != "sand":
            return self.clay_fallback(point)

        sf = self.stretching_factor(point.z_soil, point.diameter)
        curve = ThiekenPY(self._layer_params(point.layer, point.z_soil,
                                             point.sigma_v),
                          point.z_soil, point.sigma_v, point.diameter, sf)
        return curve_contribution(curve, point.length, 1.2 * curve.yc,
                                  (curve.y_0, curve.ya, curve.yb,
                                   curve.yc)), None
