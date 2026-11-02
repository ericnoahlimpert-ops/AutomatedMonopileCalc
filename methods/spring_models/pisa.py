"""PISA design model for monopiles in sand (multi-spring model).

Burd, H. J. et al. (2020). PISA design model for monopiles for offshore wind
turbines: application to a marine sand. Geotechnique 70, No. 11, 1048-1066.
https://doi.org/10.1680/jgeot.18.P.277

Four soil reaction components (Fig. 2):
    p(v)      distributed lateral load       -> lateral node springs
    m(v, psi) distributed moment             -> rotational node springs
    H_B(v)    base horizontal force          -> lateral base spring
    M_B(psi)  base moment                    -> rotational base spring

All curves use the conic function, eqs (1)-(5), in the normalised form of
Table 1, with the parameters of the General Dunkirk Sand Model (GDSM),
Table 6 (depth variation forms from Table 5).

The distributed moment is normalised by the current lateral load,
m = |p(v)| * D * m_bar(psi_bar) (eq. 6). It is therefore evaluated with the
displacement v of the previous analysis run, so the model is iterative. The
first run uses the ultimate lateral load p_u instead (upper bound of m, so
short piles that rely on m do not fail in the first run); every following
run updates m with the latest displacements until the response converges.

Pile diameter: the PISA papers use one diameter D (their pile tables) both
for the thin-walled section properties (I = pi D^3 t/8, A = pi D t) and in
all formulas (normalisations, z/D, L/D). In thin-walled theory this is the
mean (mid-surface) diameter of the wall, so this implementation uses
D_mid = D - t if DIAMETER_REFERENCE is "outer" (D in pile.csv = outer
diameter) and D_mid = D if it is "mid" (D in pile.csv = the D of the
papers). With "mid" and the D of the papers, Table 7 of Burd et al. (2020)
is reproduced within 0.2 %; treating that D as the outer diameter gives
1.2-1.6 % lower loads (see the PISA verification).

Calibration space: 2 <= L/D <= 6, 5 <= h/D <= 15, 45 % <= Dr <= 90 %,
5 m <= D <= 10 m, homogeneous sand, monotonic drained loading. Outside of
this space the normalised depths z/D and L/D and the relative density are
clamped to the calibration limits (with a warning), so the parameter
functions stay positive. Results there are an extrapolation.

Soil parameters (soil.csv):
    all layers  gamma_eff  effective unit weight [kN/m^3], required
                           -> vertical effective stress s'v used in all
                              normalisations (p, v, H_B, M_B, psi)
    sand        dr         relative density [%], required
                           -> GDSM parameters (Table 6), clamped to 45-90 %
                G0 -> normalised displacement/rotation v_bar, psi_bar,
                chosen from the given columns (base.py), priority:
                g0         G0 directly [MPa], "top-bottom" allowed
                g0_ref     G0 at s_m = 100 kPa [MPa]  power law (Thieken)
                lambda_g0  stress exponent [-]
                hb_b       Hardin & Black constant B [-]  Hardin & Black, eq. 7
                e0         void ratio [-]
                           Dunkirk sand: B = 875, e0 = 0.741 / 0.685 /
                           0.629 / 0.573 for Dr = 45 / 60 / 75 / 90 %
                phi        friction angle [deg], optional -> k0 = 1 - sin phi
                           for the mean stress s_m of the G0 calculation
                           (without phi: k0 = 0.4)
    clay        su, eps50, J as for API Matlock (PISA clay not implemented,
                API p-y curves without distributed moment are used)
"""

import math

from .base import (Contribution, SpringModel, check_g0_inputs,
                   missing_parameters, small_strain_modulus)

# Calibration space (Burd et al. 2020, Table 6)
DR_RANGE = (0.45, 0.90)
L_OVER_D_RANGE = (2.0, 6.0)
H_OVER_D_RANGE = (5.0, 15.0)
D_RANGE = (5.0, 10.0)


def conic(x, k, n, x_u, y_u):
    """Normalised soil reaction from the four-parameter conic function,
    eqs (2)-(5). x >= 0, returns y >= 0."""
    if x <= 0.0:
        return 0.0
    if x >= x_u:
        return y_u
    a = 1.0 - 2.0 * n
    b = 2.0 * n * x / x_u - (1.0 - n) * (1.0 + x * k / y_u)
    c = x * k / y_u * (1.0 - n) - n * x * x / (x_u * x_u)
    denominator = -b + math.sqrt(max(b * b - 4.0 * a * c, 0.0))
    if abs(denominator) < 1e-15:
        return y_u * x / x_u
    return y_u * 2.0 * c / denominator


class GdsmParameters(object):
    """Soil reaction curve parameters of the General Dunkirk Sand Model.
    Each method returns (k, n, x_u, y_u) for the conic function."""

    def __init__(self, dr):
        self.dr = dr  # relative density as a decimal

    def lateral(self, z_over_d, z_over_l):
        dr = self.dr
        v_pu = 146.1 - 92.11 * dr
        k_p = (8.731 - 0.6982 * dr) - 0.9178 * z_over_d
        n_p = 0.917 + 0.06193 * dr
        p_u = (0.3667 + 25.89 * dr) + (0.3375 - 8.900 * dr) * z_over_l
        return k_p, n_p, v_pu, p_u

    def moment(self, z_over_l):
        dr = self.dr
        k_m = 17.00
        m_u = 0.2605 + (-0.1989 + 0.2019 * dr) * z_over_l
        return k_m, 0.0, m_u / k_m, m_u

    def base_shear(self, l_over_d):
        dr = self.dr
        v_hu = (0.5150 + 2.883 * dr) + (0.1695 - 0.7018 * dr) * l_over_d
        k_h = (6.505 - 2.985 * dr) + (-0.007969 - 0.4299 * dr) * l_over_d
        n_h = (0.09978 + 0.7974 * dr) + (0.004994 - 0.07005 * dr) * l_over_d
        h_u = (0.09952 + 0.7996 * dr) + (0.03988 - 0.1606 * dr) * l_over_d
        return k_h, n_h, v_hu, h_u

    def base_moment(self, l_over_d):
        dr = self.dr
        psi_mu = 44.89
        k_m = 0.3515
        n_m = 0.300 + 0.4986 * dr
        m_u = (0.09981 + 0.3710 * dr) + (0.01998 - 0.09041 * dr) * l_over_d
        return k_m, n_m, psi_mu, m_u


def _clamp(value, bounds):
    return min(max(value, bounds[0]), bounds[1])


class PisaModel(SpringModel):
    """Common part of the PISA models (sand and clay): mid-surface
    diameter, L/D and the checks of the calibration space."""

    has_rotational_springs = True
    model_text = "PISA"         # prefix of the messages
    warn_layered = True         # warn for profiles with several layers

    def __init__(self, pile_sections, soil_layers, mudline_depth,
                 load_height=None, diameter_reference="outer"):
        super(PisaModel, self).__init__(pile_sections, soil_layers,
                                        mudline_depth)
        if diameter_reference not in ("outer", "mid"):
            raise ValueError("diameter_reference must be 'outer' or 'mid'.")
        # PISA uses one D for the thin-walled section and the formulas: the
        # mid-surface diameter (see module docstring)
        self.d_minus_t = diameter_reference == "outer"

        self.diameter_ref = self.mean_embedded_diameter(
            minus_thickness=self.d_minus_t)
        l_over_d = self.embedded_length / self.diameter_ref
        self.l_over_d = _clamp(l_over_d, L_OVER_D_RANGE)
        self._check_calibration_space(l_over_d, load_height)

    # -- calibration space -------------------------------------------------- #
    def _check_calibration_space(self, l_over_d, load_height):
        name = self.model_text
        if not L_OVER_D_RANGE[0] <= l_over_d <= L_OVER_D_RANGE[1]:
            self.warn("%s: L/D = %.2f outside calibration range %g-%g; "
                      "z/D and L/D are clamped to the range."
                      % ((name, l_over_d) + L_OVER_D_RANGE))
        if not D_RANGE[0] <= self.diameter_ref <= D_RANGE[1]:
            self.warn("%s: D = %.2f m outside calibration range %g-%g m."
                      % ((name, self.diameter_ref) + D_RANGE))
        if load_height is not None:
            h_over_d = load_height / self.diameter_ref
            if not H_OVER_D_RANGE[0] <= h_over_d <= H_OVER_D_RANGE[1]:
                self.warn("%s: h/D = %.2f outside calibration range "
                          "%g-%g." % ((name, h_over_d) + H_OVER_D_RANGE))
        if len(self.soil_layers) > 1 and self.warn_layered:
            self.warn("%s: the model is calibrated for a homogeneous soil; "
                      "the layered profile is evaluated layer by layer."
                      % name)

    # -- soil properties ---------------------------------------------------- #
    def mid_diameter(self, point):
        """Mid-surface diameter D_mid used in all PISA formulas: D - t if D
        in pile.csv is the outer diameter, otherwise D itself."""
        if self.d_minus_t:
            return point.diameter - point.thickness
        return point.diameter


class PisaSandModel(PisaModel):
    """PISA GDSM sand model with distributed moment and base springs.

    Clay layers are not covered by the sand formulation; they fall back to
    API Matlock p-y curves without distributed moment (with a warning).
    """

    name = "PISA_SAND"
    label = "PISA_Sand"
    iterative = True

    def __init__(self, pile_sections, soil_layers, mudline_depth,
                 load_height=None, diameter_reference="outer", **_):
        errors, notes = check_g0_inputs("PISA", soil_layers)
        errors += missing_parameters(soil_layers, ("dr",))
        if errors:
            raise ValueError("PISA: soil parameters in soil.csv (%s)."
                             % "; ".join(errors))
        super(PisaSandModel, self).__init__(pile_sections, soil_layers,
                                            mudline_depth, load_height,
                                            diameter_reference)
        self.notes.extend(notes)
        for layer in soil_layers:
            if layer.soil_type == "sand" and \
                    not 100 * DR_RANGE[0] <= layer.dr <= 100 * DR_RANGE[1]:
                self.warn("PISA: Dr = %g%% outside calibration range "
                          "45-90%%; clamped." % layer.dr)

    # -- soil properties ---------------------------------------------------- #
    def _sand_state(self, point):
        """GDSM parameters, s'v [kPa] and G0 [kPa] at a soil point."""
        layer = point.layer
        dr = _clamp(layer.dr / 100.0, DR_RANGE)
        sigma_v = max(point.sigma_v, 1e-3)
        g0 = small_strain_modulus(layer, point.z_soil, sigma_v)
        return GdsmParameters(dr), sigma_v, g0

    # -- soil reactions ----------------------------------------------------- #
    def distributed_springs(self, point, node_result):
        if point.layer.soil_type != "sand":
            return self.clay_fallback(point)
        return self.sand_distributed_springs(point, node_result)

    def base_springs(self, point, node_result):
        if point.layer.soil_type != "sand":
            self.warn("PISA: pile tip in clay, no base springs applied.")
            return None, None
        return self.sand_base_springs(point)

    def sand_distributed_springs(self, point, node_result):
        """(lateral, rotational) of a sand soil point."""
        params, sigma_v, g0 = self._sand_state(point)
        d, length = self.mid_diameter(point), point.length
        z_over_d = _clamp(point.z_soil / d, (0.0, L_OVER_D_RANGE[1]))
        z_over_l = _clamp(point.z_soil / self.embedded_length, (0.0, 1.0))

        # Distributed lateral load: p = s'v D p_bar(v_bar), v_bar = v G0/(D s'v)
        k_p, n_p, v_pu, p_u = params.lateral(z_over_d, z_over_l)
        v_scale = d * sigma_v / g0          # v = v_bar * v_scale

        def p(v):                           # [kN/m]
            return sigma_v * d * conic(v / v_scale, k_p, n_p, v_pu, p_u)

        lateral = Contribution(lambda v: 1000.0 * length * p(v),
                               1.2 * v_pu * v_scale, (v_pu * v_scale,))

        # Distributed moment: m = |p(v)| D m_bar(psi_bar), psi_bar = psi G0/s'v
        # First run: p = p_u (upper bound), afterwards p at the last v
        rotational = None
        if node_result is None:
            p_current = sigma_v * d * p_u
        else:
            p_current = p(abs(node_result[0]))
        if p_current > 0.0:
            k_m, n_m, psi_mu, m_u = params.moment(z_over_l)
            psi_scale = sigma_v / g0
            rotational = Contribution(
                lambda psi: 1000.0 * length * p_current * d
                * conic(psi / psi_scale, k_m, n_m, psi_mu, m_u),
                1.2 * psi_mu * psi_scale, (psi_mu * psi_scale,))
        return lateral, rotational

    def sand_base_springs(self, point):
        """(shear, moment) at the pile base in a sand layer."""
        params, sigma_v, g0 = self._sand_state(point)
        d = self.mid_diameter(point)

        # H_B = s'v D^2 H_bar(v_bar)
        k_h, n_h, v_hu, h_u = params.base_shear(self.l_over_d)
        v_scale = d * sigma_v / g0
        shear = Contribution(
            lambda v: 1000.0 * sigma_v * d ** 2
            * conic(v / v_scale, k_h, n_h, v_hu, h_u),
            1.2 * v_hu * v_scale, (v_hu * v_scale,))

        # M_B = s'v D^3 M_bar(psi_bar)
        k_mb, n_mb, psi_mu, m_u = params.base_moment(self.l_over_d)
        psi_scale = sigma_v / g0
        moment = Contribution(
            lambda psi: 1000.0 * sigma_v * d ** 3
            * conic(psi / psi_scale, k_mb, n_mb, psi_mu, m_u),
            1.2 * psi_mu * psi_scale, (psi_mu * psi_scale,))
        return shear, moment
