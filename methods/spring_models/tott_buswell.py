"""Tott-Buswell et al. (2024): CPT-based multi-spring model for
monopiles in sand (SLS).

Tott-Buswell, J., Prendergast, L. J., Gavin, K. (2024). A CPT-based
multi-spring model for lateral monopile analysis under SLS conditions in
sand. Ocean Engineering 293, 116642.
https://doi.org/10.1016/j.oceaneng.2023.116642

Four soil reaction components (Fig. 1), all informed by the cone resistance
qc of a CPT profile (cpt.csv, see csv_reader.read_cpt_csv):

    p(y)        distributed lateral load, Li et al. (2014), eq. (1):
                    p = 3.6 D (g' D) (qc / (g' D))^0.72 (y / D)^0.66
    m(theta)    distributed moment from the vertical shaft friction (axial
                tau-w curve of Lehane et al. 2020 integrated around the
                circumference), eqs. (2)-(7):
                    m = a theta - b theta^2  (theta < theta_f), m_f beyond
                    a = pi D^2 G0 / 16,   b = G0^2 D^2 / (96 tau_f)
                    theta_f = 3 pi tau_f / G0,  m_f = 3 pi^2 D^2 tau_f / 32
                with the shaft friction of the unified CPT method, eq. (3):
                    tau_f = (ft/fc) (s'rc + ds'rd) tan delta_f
                    s'rc  = (qc/44) Are^0.3 max(1, H/D)^-0.4
                    ds'rd = (qc/10) (qc/s'v)^-0.33 (d_CPT/D)
                    Are = 1 - PLR (Di/D)^2,  PLR = tanh(0.3 (Di/d_CPT)^0.5)
                H = distance from the pile tip, ft/fc = 0.8, delta_f = 29 deg,
                d_CPT = 35.7 mm
    M_B(theta)  base moment, bilinear, eqs. (8), (9):
                    M_B,f = q_b,res pi D^2/8 * 2D/(3 pi) = q_b,res D^3 / 12
                    theta_B,f = 44.98 s'v,b / G0
    H_B(y)      base shear, bilinear, eqs. (10)-(12):
                    H_B,f = pi D^2/4 q_b,res tan 35 deg
                    y_B,f = 2 y_bar D s'v,b / G0
                    y_bar = (0.52 + 2.88 Dr) + (0.17 - 0.70 Dr) L/D,
                    L/D limited to 2 ... 6

The residual base stress is q_b,res = alpha * qc,r, qc,r = mean qc from
0.25 L/D [m] above to 0.25 L/D [m] below the pile tip;
alpha = TOTT_BUSWELL_ALPHA (Input.py, paper: 0.1 recommended, 0.04 ... 0.1
fitted to field tests).

The springs depend on the depth only (no iteration). qc is averaged over the
pile length of every spring half element (the paper averages qc over its
0.05 m elements); with strongly varying qc use short elements.

Validity (paper): small to medium displacements (SLS, about 0.01 D at the
mudline). The power-law p-y curve has no ultimate value, so the model does
not predict the ultimate capacity. p-y curves are sampled up to
y = 0.5 D (constant beyond).

Deviations from the paper and choices of this implementation:
- D is the outer diameter (D + t if DIAMETER_REFERENCE = "mid"); the inner
  diameter is Di = D - 2t (the paper writes D - t).
- s'v is the vertical effective stress of soil.csv (layered), g' in eq. (1)
  the gamma_eff of the layer.
- Dr for the base shear (eq. 12): dr of the soil layer at the tip, 75 %
  (paper) if not given.
- Clay layers use the API Matlock p-y curves (lateral only, with a
  warning); no base springs if the tip is in clay.

Soil parameters (soil.csv):
    all layers  gamma_eff  effective unit weight [kN/m^3], required
                           -> s'v, g' in eq. (1)
    sand        G0 only if the CPT file has no g0 column; chosen from the
                given columns (base.py), priority:
                g0         G0 directly [MPa], "top-bottom" allowed
                g0_ref     G0 at s_m = 100 kPa [MPa]  power law (Thieken)
                lambda_g0  stress exponent [-]
                hb_b       Hardin & Black constant B [-]
                e0         void ratio [-]
                phi        friction angle [deg], optional -> k0 = 1 - sin phi
                           for the mean stress of the G0 calculation
                dr         relative density [%], optional (base shear,
                           default 75 %)
    clay        su, eps50, J as for API Matlock
CPT file (CPT_CSV in Input.py): z [m below mudline], qc [MPa],
optional g0 [MPa] (replaces the G0 of soil.csv).
"""

import math

from .base import (Contribution, SpringModel, check_g0_inputs,
                   effective_stress, layer_at, section_at,
                   small_strain_modulus)

# p-y curve, Li et al. (2014)
PY_COEFFICIENT = 3.6
PY_QC_EXPONENT = 0.72
PY_Y_EXPONENT = 0.66
PY_RANGE = 0.5             # curves sampled up to y = PY_RANGE * D

# Shaft friction, unified CPT method (Lehane et al. 2020)
FT_FC = 0.8                # loading configuration ratio ft/fc
DELTA_F = 29.0             # interface friction angle [deg]
D_CPT = 0.0357             # diameter of the CPT probe [m]

# Base springs
TIP_RANGE = 0.25           # qc,r averaged over +-TIP_RANGE * L/D [m]
THETA_B_BAR = 44.98        # normalised failure rotation of the base moment
BASE_PHI = 35.0            # friction angle in the annulus [deg]
DR_DEFAULT = 75.0          # relative density for the base shear [%]
L_OVER_D_RANGE = (2.0, 6.0)
DEFAULT_ALPHA = 0.1        # q_b,res / qc,r

COMPONENTS = ("p-y", "m-theta", "base")


def bilinear(x, x_f, f_f):
    """Bilinear curve: linear up to (x_f, f_f), constant beyond."""
    return f_f * min(x / x_f, 1.0)


def shaft_friction(qc, sigma_v, diameter, inner_diameter, h):
    """Ultimate shaft friction tau_f [kPa] of eq. (3). qc, sigma_v in kPa,
    h = distance from the pile tip [m]."""
    plr = math.tanh(0.3 * math.sqrt(inner_diameter / D_CPT))
    are = 1.0 - plr * (inner_diameter / diameter) ** 2
    sigma_rc = qc / 44.0 * are ** 0.3 * max(1.0, h / diameter) ** -0.4
    sigma_rd = (qc / 10.0 * (qc / max(sigma_v, 1e-3)) ** -0.33
                * D_CPT / diameter)
    return FT_FC * (sigma_rc + sigma_rd) * math.tan(math.radians(DELTA_F))


class TottBuswellModel(SpringModel):
    """Tott-Buswell CPT-based multi-spring model (p-y, m-theta, base
    moment, base shear)."""

    name = "TOTT_BUSWELL"
    label = "Tott-Buswell"
    iterative = False
    has_rotational_springs = True
    output_columns = ("qc [MPa]", "G0 [MPa]", "tau_f [kPa]")

    def __init__(self, pile_sections, soil_layers, mudline_depth,
                 cpt_profile=None, tott_buswell_alpha=DEFAULT_ALPHA,
                 tott_buswell_components=COMPONENTS,
                 diameter_reference="outer", **_):
        super(TottBuswellModel, self).__init__(pile_sections, soil_layers,
                                               mudline_depth)
        if cpt_profile is None:
            raise ValueError("TOTT_BUSWELL needs a CPT profile: set "
                             "CPT_CSV in Input.py.")
        if diameter_reference not in ("outer", "mid"):
            raise ValueError("diameter_reference must be 'outer' or 'mid'.")
        if not tott_buswell_alpha > 0.0:
            raise ValueError("TOTT_BUSWELL_ALPHA must be positive, got %r."
                             % (tott_buswell_alpha,))
        components = tuple(c.lower() for c in tott_buswell_components)
        unknown = [c for c in components if c not in COMPONENTS]
        if unknown or not components:
            raise ValueError("TOTT_BUSWELL_COMPONENTS: use one or more of %s, "
                             "got %r." % (", ".join(COMPONENTS),
                                          tott_buswell_components))
        self.cpt = cpt_profile
        self.alpha = float(tott_buswell_alpha)
        self.components = components
        self.plus_thickness = diameter_reference == "mid"
        self.has_rotational_springs = ("m-theta" in components
                                       or "base" in components)

        if not self.cpt.has_g0:
            errors, notes = check_g0_inputs("TOTT_BUSWELL", soil_layers)
            if errors:
                raise ValueError("TOTT_BUSWELL: G0 needed (g0 column in "
                                 "the CPT file or G0 parameters in "
                                 "soil.csv): %s." % "; ".join(errors))
            self.notes.extend(notes)
        self._prepare_base()
        self.notes.append(
            "TOTT_BUSWELL: CPT %s (%.2f m), G0 from %s, components %s"
            % (self.cpt.path, self.cpt.z_max,
               "the CPT file" if self.cpt.has_g0 else "soil.csv",
               ", ".join(self.components)))
        self.notes.append(
            "TOTT_BUSWELL: q_b,res = %.3g * qc,r = %.4g MPa (qc,r = mean "
            "qc over %.3f ... %.3f m)" % (self.alpha, self.qb_res / 1000.0,
                                          self.tip_range[0],
                                          self.tip_range[1]))
        self.notes.append("TOTT_BUSWELL: valid for small to medium "
                          "displacements (SLS), no ultimate capacity.")
        if self.cpt.z_max < self.embedded_length - 1e-6:
            self.warn("TOTT_BUSWELL: the CPT ends at %.2f m, above the "
                      "pile tip (%.2f m); qc below is taken constant."
                      % (self.cpt.z_max, self.embedded_length))

    # -- soil properties ---------------------------------------------------- #
    def outer_diameter(self, z_pile):
        """Outer diameter D and wall thickness t at a depth below the pile
        head."""
        section = section_at(self.pile_sections, z_pile)
        d, t = section.diameter_at(z_pile), section.thickness_at(z_pile)
        return (d + t if self.plus_thickness else d), t

    def shear_modulus(self, layer, z0, z1, sigma_v):
        """G0 [kPa]: mean of the CPT g0 between z0 and z1, or from the soil
        layer at the mid-depth."""
        if self.cpt.has_g0:
            return self.cpt.average(z0, z1, "g0")
        return small_strain_modulus(layer, 0.5 * (z0 + z1), sigma_v)

    def _prepare_base(self):
        """Residual base stress q_b,res and L/D (base springs)."""
        tip = self.pile_tip
        d, _ = self.outer_diameter(tip)
        self.l_over_d = self.embedded_length / d
        half_range = TIP_RANGE * self.l_over_d
        self.tip_range = (max(self.embedded_length - half_range, 0.0),
                          self.embedded_length + half_range)
        self.qb_res = self.alpha * self.cpt.average(*self.tip_range)
        if self.cpt.z_max < self.tip_range[1] - 1e-6 \
                and "base" in self.components:
            self.warn("TOTT_BUSWELL: the CPT ends at %.2f m, qc,r needs "
                      "qc down to %.2f m; qc below the CPT is taken "
                      "constant." % (self.cpt.z_max, self.tip_range[1]))
        if not L_OVER_D_RANGE[0] <= self.l_over_d <= L_OVER_D_RANGE[1] \
                and "base" in self.components:
            self.notes.append("TOTT_BUSWELL: L/D = %.2f limited to %g-%g "
                              "for the base shear displacement (eq. 12)"
                              % ((self.l_over_d,) + L_OVER_D_RANGE))

    # -- soil reactions ----------------------------------------------------- #
    def distributed_springs(self, point, node_result):
        if point.layer.soil_type != "sand":
            return self.clay_fallback(point)

        d, t = self.outer_diameter(point.z_pile)
        z0 = point.z_soil - 0.5 * point.length
        z1 = point.z_soil + 0.5 * point.length
        qc = self.cpt.average(z0, z1)
        length = point.length

        lateral = None
        if "p-y" in self.components:
            gamma_d = point.layer.gamma_eff * d
            p_ref = (PY_COEFFICIENT * d * gamma_d
                     * (qc / gamma_d) ** PY_QC_EXPONENT)       # [kN/m]
            lateral = Contribution(
                lambda y: 1000.0 * length * p_ref * (y / d) ** PY_Y_EXPONENT,
                PY_RANGE * d)

        rotational = None
        if "m-theta" in self.components:
            sigma_v = max(point.sigma_v, 1e-3)
            g0 = self.shear_modulus(point.layer, z0, z1, sigma_v)
            tau_f = shaft_friction(qc, sigma_v, d, d - 2.0 * t,
                                   self.embedded_length - point.z_soil)
            if g0 > 0.0 and tau_f > 0.0:
                a = math.pi * d ** 2 * g0 / 16.0
                b = g0 ** 2 * d ** 2 / (96.0 * tau_f)
                theta_f = 3.0 * math.pi * tau_f / g0
                m_f = a * a / (4.0 * b)

                def moment(theta):              # [Nm per pile length]
                    m = a * theta - b * theta ** 2 if theta < theta_f else m_f
                    return 1000.0 * length * m

                rotational = Contribution(moment, 1.2 * theta_f, (theta_f,))
        return lateral, rotational

    def base_springs(self, point, node_result):
        if "base" not in self.components:
            return None, None
        if point.layer.soil_type != "sand":
            self.warn("TOTT_BUSWELL: pile tip in clay, no base springs "
                      "applied.")
            return None, None

        d, _ = self.outer_diameter(point.z_pile)
        sigma_v = max(point.sigma_v, 1e-3)
        g0 = self.shear_modulus(point.layer, point.z_soil, point.z_soil,
                                sigma_v)
        dr = (point.layer.dr if point.layer.dr is not None
              else DR_DEFAULT) / 100.0
        l_over_d = min(max(self.l_over_d, L_OVER_D_RANGE[0]),
                       L_OVER_D_RANGE[1])

        # Base shear: H_B,f = pi D^2/4 q_b,res tan 35, y_B,f eq. (11)
        h_f = (math.pi * d ** 2 / 4.0 * self.qb_res
               * math.tan(math.radians(BASE_PHI)))
        y_bar = (0.52 + 2.88 * dr) + (0.17 - 0.70 * dr) * l_over_d
        y_f = 2.0 * y_bar * d * sigma_v / g0
        shear = Contribution(lambda y: 1000.0 * bilinear(y, y_f, h_f),
                             1.2 * y_f, (y_f,))

        # Base moment: M_B,f = q_b,res D^3 / 12, theta_B,f eq. (9)
        m_f = self.qb_res * d ** 3 / 12.0
        theta_f = THETA_B_BAR * sigma_v / g0
        moment = Contribution(lambda th: 1000.0 * bilinear(th, theta_f, m_f),
                              1.2 * theta_f, (theta_f,))
        return shear, moment

    def node_output(self, z_soil, diameter):
        """qc, G0 and tau_f at the node depth (the springs use the mean qc
        of the half elements above and below the node)."""
        z_pile = z_soil + self.mudline_depth
        layer = layer_at(self.soil_layers, z_soil)
        if layer is None or layer.soil_type != "sand":
            return {"qc [MPa]": self.cpt.value_at(z_soil) / 1000.0}
        d, t = self.outer_diameter(z_pile)
        qc = self.cpt.value_at(z_soil)
        sigma_v = max(effective_stress(self.soil_layers, z_soil), 1e-3)
        g0 = self.shear_modulus(layer, z_soil, z_soil, sigma_v)
        tau_f = shaft_friction(qc, sigma_v, d, d - 2.0 * t,
                               self.embedded_length - z_soil)
        return {"qc [MPa]": qc / 1000.0, "G0 [MPa]": g0 / 1000.0,
                "tau_f [kPa]": tau_f}
