"""PISA design model for layered soils (sand and clay).

Burd, H. J. et al. (2020b). Application of the PISA design model to
monopiles embedded in layered soils. Geotechnique 70, No. 11, 1067-1082.
https://doi.org/10.1680/jgeot.20.PISA.009

The soil reaction curves calibrated for homogeneous soils are applied
locally at each level of the pile: the soil type of the layer selects the
framework (Table 11),
    sand  GDSM, Burd et al. (2020a), see pisa.py:
          normalised with s'v and G0, distributed moment m = |p| D m_bar
    clay  Cowden till (Byrne et al. 2020a) or Bothkennar clay (Burd et al.
          2020b, Table 12), soil.csv column pisa_clay, see pisa_clay.py:
          normalised with su and G0
and the local s'v, su and G0 of the layer are used. The base springs use the
soil layer at the pile tip. Layers below the pile tip have no influence (as
in the paper). The model is iterative if the profile contains sand (the
sand distributed moment depends on p).

Validity (paper): the 1D model matches 3D finite-element analyses for most
layered profiles (small displacements: eta_sd >= 0.9 in all cases). It is
less accurate at large displacements for combinations of very soft clay and
very dense sand, especially for L/D = 2 (warning).

Soil parameters (soil.csv): per layer as for PISA_SAND (sand layers) and
PISA_CLAY (clay layers):
    sand  dr, G0 input (g0 / g0_ref + lambda_g0 / hb_b + e0), phi optional
    clay  su, G0 input, pisa_clay ("cowden" default or "bothkennar")
    all   gamma_eff
"""

from .base import check_g0_inputs, missing_parameters
from .pisa import DR_RANGE, PisaModel, PisaSandModel
from .pisa_clay import (PisaClayModel, check_clay_calibrations,
                        clay_calibration_name)

# Combinations the paper found difficult (eta_ult < 0.9): soft Bothkennar
# clay with very dense sand
VERY_DENSE_DR = 85.0          # [%]


class PisaLayeredModel(PisaSandModel, PisaClayModel):
    """PISA model for layered profiles of sand and clay layers."""

    name = "PISA"
    label = "PISA"
    model_text = "PISA"
    warn_layered = False

    def __init__(self, pile_sections, soil_layers, mudline_depth,
                 load_height=None, diameter_reference="outer", **_):
        errors, notes = [], []
        for soil_type in ("sand", "clay"):
            g0_errors, g0_notes = check_g0_inputs("PISA", soil_layers,
                                                  soil_type)
            errors += g0_errors
            notes += g0_notes
        errors += missing_parameters(soil_layers, ("dr",))
        clay_errors, clay_notes = check_clay_calibrations("PISA", soil_layers)
        errors += clay_errors
        notes += clay_notes
        if errors:
            raise ValueError("PISA: soil parameters in soil.csv (%s)."
                             % "; ".join(errors))
        PisaModel.__init__(self, pile_sections, soil_layers, mudline_depth,
                           load_height, diameter_reference)
        self.notes.extend(notes)

        sands = [l for l in soil_layers if l.soil_type == "sand"]
        clays = [l for l in soil_layers if l.soil_type == "clay"]
        # The sand distributed moment depends on p(v): iterate if needed
        self.iterative = bool(sands)
        for layer in sands:
            if not 100 * DR_RANGE[0] <= layer.dr <= 100 * DR_RANGE[1]:
                self.warn("PISA: Dr = %g%% outside calibration range "
                          "45-90%%; clamped." % layer.dr)
        if sands and clays:
            self.notes.append("PISA: layered profile, sand and clay "
                              "reaction curves applied layer by layer "
                              "(Burd et al. 2020b)")
        soft = [l for l in clays if clay_calibration_name(l) == "bothkennar"]
        dense = [l for l in sands if l.dr >= VERY_DENSE_DR]
        if soft and dense:
            self.warn("PISA: soft clay (Bothkennar) combined with very dense "
                      "sand (Dr >= %g%%): Burd et al. (2020b) found the "
                      "ultimate response less accurate for such profiles, "
                      "especially for L/D = 2; check with 3D FE analysis."
                      % VERY_DENSE_DR)

    # -- soil reactions ----------------------------------------------------- #
    def distributed_springs(self, point, node_result):
        if point.layer.soil_type == "sand":
            return self.sand_distributed_springs(point, node_result)
        return self.clay_distributed_springs(point)

    def base_springs(self, point, node_result):
        if point.layer.soil_type == "sand":
            return self.sand_base_springs(point)
        return self.clay_base_springs(point)
