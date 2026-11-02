"""API p-y curves with the initial stiffness of Kallehave et al. (2012),
"Optimization of monopiles for offshore wind turbines".

The API sand curve (see api.py)

    p = A * pu * tanh(E_py * y / (A * pu))

is kept, but the initial slope E_py = k * z of API is replaced by

    E_py = k * z0 * (z / z0)^m * (D / D0)^0.5
    z0 = 2.5 m, D0 = 0.61 m, m = 0.6

with k from API. Clay layers use the API Matlock curves unchanged. The
springs do not depend on the pile response (no iteration).

Units inside the curve class: kN, m, kPa.

Soil parameters (soil.csv):
    all layers  gamma_eff  effective unit weight [kN/m^3], required
    sand        phi        friction angle [deg], required -> pu
                k          initial modulus of subgrade reaction [MN/m^3],
                           optional (default from the API chart, see api.py)
    clay        su, eps50, J as for API (api.py)
"""

from .api import ApiPyModel, SandPY, api_sand_k

Z0 = 2.5      # [m]
D0 = 0.61     # [m]
M = 0.6


class KallehaveSandPY(SandPY):
    def __init__(self, layer, z, sigma_v, diameter, loading):
        super(KallehaveSandPY, self).__init__(layer, z, sigma_v, diameter,
                                              loading)
        k = layer.k * 1000.0 if layer.k is not None else api_sand_k(layer.phi)
        self.kz = k * Z0 * (z / Z0) ** M * (diameter / D0) ** 0.5


class KallehavePyModel(ApiPyModel):
    """API p-y springs with the Kallehave sand stiffness, no iteration."""

    name = "KALLEHAVE"
    label = "Kallehave"
    sand_curve = KallehaveSandPY
