"""PISA design model for monopiles in clay (multi-spring model).

Byrne, B. W. et al. (2020). PISA design model for monopiles for offshore
wind turbines: application to a stiff glacial clay till. Geotechnique 70,
No. 11, 1030-1047. https://doi.org/10.1680/jgeot.18.P.255

Four soil reaction components (Fig. 2), as in the sand model (pisa.py):
    p(v)      distributed lateral load       -> lateral node springs
    m(psi)    distributed moment             -> rotational node springs
    H_B(v)    base horizontal force          -> lateral base spring
    M_B(psi)  base moment                    -> rotational base spring

All curves use the conic function, eqs (17)-(21), normalised with the
undrained shear strength su and the small-strain shear modulus G0 (Table 3):

    p_bar = p / (su D)        v_bar   = v G0 / (su D)
    m_bar = m / (su D^2)      psi_bar = psi G0 / su
    H_bar = H_B / (su D^2)    M_bar   = M_B / (su D^3)

with the second-stage parameters of the Cowden till model (Table 4), which
the paper recommends for design. Alternatively the Bothkennar clay (soft
clay) calibration of Burd et al. (2020b), Table 12, is available
(AWG 2018; soil.csv column pisa_clay = "bothkennar"):

Burd, H. J. et al. (2020b). Application of the PISA design model to
monopiles embedded in layered soils. Geotechnique 70, No. 11, 1067-1082.
https://doi.org/10.1680/jgeot.20.PISA.009 Unlike in sand, the distributed moment
does not depend on p, so the springs depend on the depth only and the model
is not iterative. su and G0 are taken at the depth of each half element
(base springs: at the pile tip).

D in all formulas is the mid-surface diameter, as in the sand model (see
pisa.py): D - t if DIAMETER_REFERENCE is "outer", D if it is "mid".

Calibration space: 2 <= L/D <= 6, 5 <= h/D <= 15, 5 m <= D <= 10 m,
undrained loading of the clay with the su and G0 profiles of the reference
site (Cowden till: Fig. 4; Bothkennar clay: Burd et al. 2020b, Fig. 3). Outside of
it, z/D and L/D are clamped (with a warning). The parameters are specific
to these profiles; other clays or strongly different su/G0 profiles need a
separate calibration (paper, Discussion).

Soil parameters (soil.csv):
    clay        su         undrained shear strength (triaxial compression)
                           [kPa], "top-bottom" allowed, required
                pisa_clay  calibration: "cowden" (default) or "bothkennar"
                G0 -> normalised displacement/rotation, chosen from the
                given columns (base.py), priority:
                g0         G0 directly [MPa], "top-bottom" allowed
                g0_ref     G0 at s_m = 100 kPa [MPa]  power law; the paper
                lambda_g0  stress exponent [-]       uses G0 = 1100 p', i.e.
                                                     g0_ref = 110, lambda = 1
                hb_b, e0   Hardin & Black
                phi        optional, k0 = 1 - sin phi for the mean stress of
                           the G0 calculation (without phi: k0 = 0.4)
    all layers  gamma_eff  effective unit weight [kN/m^3], required by
                           soil.csv (used only for the G0 calculation)
    sand        not covered by the clay model (error)
"""

import math

from .base import Contribution, check_g0_inputs, small_strain_modulus
from .pisa import L_OVER_D_RANGE, PisaModel, _clamp, conic


class CowdenParameters(object):
    """Soil reaction curve parameters of the Cowden till model, second-stage
    calibration (Table 4). Each method returns (k, n, x_u, y_u) for the
    conic function."""

    @staticmethod
    def lateral(z_over_d):
        k_p = 10.60 - 1.650 * z_over_d
        n_p = 0.9390 - 0.03345 * z_over_d
        p_u = 10.70 - 7.101 * math.exp(-0.3085 * z_over_d)
        return k_p, n_p, 241.4, p_u

    @staticmethod
    def moment(z_over_d):
        # Bilinear (n = 0); ultimate rotation given by m_u / k_m
        k_m = 1.420 - 0.09643 * z_over_d
        m_u = 0.2899 - 0.04775 * z_over_d
        return k_m, 0.0, m_u / k_m, m_u

    @staticmethod
    def base_shear(l_over_d):
        k_h = 2.717 - 0.3575 * l_over_d
        n_h = 0.8793 - 0.03150 * l_over_d
        h_u = 0.4038 + 0.04812 * l_over_d
        return k_h, n_h, 235.7, h_u

    @staticmethod
    def base_moment(l_over_d):
        k_m = 0.2146 - 0.002132 * l_over_d
        n_m = 1.079 - 0.1087 * l_over_d
        m_u = 0.8192 - 0.08588 * l_over_d
        return k_m, n_m, 173.1, m_u


class BothkennarParameters(object):
    """Soil reaction curve parameters of the Bothkennar clay model (soft
    clay, AWG 2018), as given in Burd et al. (2020b), Table 12."""

    @staticmethod
    def lateral(z_over_d):
        k_p = 12.05 - 1.547 * z_over_d
        n_p = 0.7204 - 0.002679 * z_over_d
        p_u = 7.743 - 3.945 * math.exp(-0.8456 * z_over_d)
        return k_p, n_p, 173.8, p_u

    @staticmethod
    def moment(z_over_d):
        k_m = 1.698 - 0.1576 * z_over_d
        m_u = 0.4862 - 0.05674 * z_over_d
        return k_m, 0.0, m_u / k_m, m_u

    @staticmethod
    def base_shear(l_over_d):
        k_h = 3.008 - 0.2701 * l_over_d
        n_h = 0.3113 + 0.04263 * l_over_d
        h_u = 0.5279 + 0.06864 * l_over_d
        return k_h, n_h, 291.5, h_u

    @staticmethod
    def base_moment(l_over_d):
        k_m = 0.3409 - 0.01995 * l_over_d
        n_m = 0.6990 - 0.1155 * l_over_d
        m_u = 0.8756 - 0.09195 * l_over_d
        return k_m, n_m, 187.0, m_u


# Clay calibrations, selected per layer with the soil.csv column pisa_clay
CLAY_CALIBRATIONS = {
    "cowden": CowdenParameters,
    "bothkennar": BothkennarParameters,
}
DEFAULT_CLAY = "cowden"
CLAY_DESCRIPTION = {
    "cowden": "Cowden till (Byrne et al. 2020, Table 4, second stage)",
    "bothkennar": "Bothkennar clay (Burd et al. 2020b, Table 12)",
}


def clay_calibration_name(layer):
    """Name of the PISA clay calibration of a clay layer."""
    return layer.pisa_clay or DEFAULT_CLAY


def check_clay_calibrations(model_name, soil_layers):
    """(errors, notes) for the pisa_clay entries of all clay layers."""
    errors, notes = [], []
    for i, layer in enumerate(soil_layers, start=1):
        if layer.soil_type != "clay":
            continue
        name = clay_calibration_name(layer)
        if name not in CLAY_CALIBRATIONS:
            errors.append("soil layer %d: pisa_clay '%s' unknown (use %s)"
                          % (i, name, " or ".join(sorted(CLAY_CALIBRATIONS))))
        else:
            notes.append("%s: soil layer %d: %s"
                         % (model_name, i, CLAY_DESCRIPTION[name]))
    return errors, notes


class PisaClayModel(PisaModel):
    """PISA clay model (Cowden till) with distributed moment and base
    springs. Sand layers are not covered (error)."""

    name = "PISA_CLAY"
    label = "PISA_Clay"
    iterative = False
    model_text = "PISA_CLAY"

    def __init__(self, pile_sections, soil_layers, mudline_depth,
                 load_height=None, diameter_reference="outer", **_):
        sands = [i for i, layer in enumerate(soil_layers, start=1)
                 if layer.soil_type == "sand"]
        if sands:
            raise ValueError("PISA_CLAY: soil layer(s) %s are sand; the "
                             "clay model covers clay layers only (use "
                             "PISA_SAND for sand)."
                             % ", ".join(str(i) for i in sands))
        errors, notes = check_g0_inputs("PISA_CLAY", soil_layers, "clay")
        calibration_errors, calibration_notes = check_clay_calibrations(
            "PISA_CLAY", soil_layers)
        errors += calibration_errors
        notes += calibration_notes
        if errors:
            raise ValueError("PISA_CLAY: soil parameters in soil.csv (%s)."
                             % "; ".join(errors))
        super(PisaClayModel, self).__init__(pile_sections, soil_layers,
                                            mudline_depth, load_height,
                                            diameter_reference)
        self.notes.extend(notes)

    # -- soil properties ---------------------------------------------------- #
    def _clay_state(self, point):
        """su [kPa] and G0 [kPa] at a soil point."""
        su = max(point.layer.su_at(point.z_soil), 1e-3)
        g0 = max(small_strain_modulus(point.layer, point.z_soil,
                                      point.sigma_v), 1e-3)
        return su, g0

    # -- soil reactions ----------------------------------------------------- #
    def distributed_springs(self, point, node_result):
        return self.clay_distributed_springs(point)

    def base_springs(self, point, node_result):
        return self.clay_base_springs(point)

    def clay_distributed_springs(self, point):
        """(lateral, rotational) of a clay soil point."""
        params = CLAY_CALIBRATIONS[clay_calibration_name(point.layer)]
        su, g0 = self._clay_state(point)
        d, length = self.mid_diameter(point), point.length
        z_over_d = _clamp(point.z_soil / d, (0.0, L_OVER_D_RANGE[1]))

        # p = su D p_bar(v_bar), v_bar = v G0 / (su D)
        k_p, n_p, v_pu, p_u = params.lateral(z_over_d)
        v_scale = d * su / g0               # v = v_bar * v_scale
        lateral = Contribution(
            lambda v: 1000.0 * length * su * d
            * conic(v / v_scale, k_p, n_p, v_pu, p_u),
            1.2 * v_pu * v_scale, (v_pu * v_scale,))

        # m = su D^2 m_bar(psi_bar), psi_bar = psi G0 / su
        k_m, n_m, psi_mu, m_u = params.moment(z_over_d)
        psi_scale = su / g0
        rotational = Contribution(
            lambda psi: 1000.0 * length * su * d ** 2
            * conic(psi / psi_scale, k_m, n_m, psi_mu, m_u),
            1.2 * psi_mu * psi_scale, (psi_mu * psi_scale,))
        return lateral, rotational

    def clay_base_springs(self, point):
        """(shear, moment) at the pile base in a clay layer."""
        params = CLAY_CALIBRATIONS[clay_calibration_name(point.layer)]
        su, g0 = self._clay_state(point)
        d = self.mid_diameter(point)

        # H_B = su D^2 H_bar(v_bar)
        k_h, n_h, v_hu, h_u = params.base_shear(self.l_over_d)
        v_scale = d * su / g0
        shear = Contribution(
            lambda v: 1000.0 * su * d ** 2
            * conic(v / v_scale, k_h, n_h, v_hu, h_u),
            1.2 * v_hu * v_scale, (v_hu * v_scale,))

        # M_B = su D^3 M_bar(psi_bar)
        k_mb, n_mb, psi_mu, m_u = params.base_moment(self.l_over_d)
        psi_scale = su / g0
        moment = Contribution(
            lambda psi: 1000.0 * su * d ** 3
            * conic(psi / psi_scale, k_mb, n_mb, psi_mu, m_u),
            1.2 * psi_mu * psi_scale, (psi_mu * psi_scale,))
        return shear, moment
