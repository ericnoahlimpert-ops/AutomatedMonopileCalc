"""API p-y curves with the initial stiffness of Sorensen et al. (2010),
"Effects of diameter on initial stiffness of p-y curves for large-diameter
piles in sand".

The API sand curve (see api.py)

    p = A * pu * tanh(E_py * y / (A * pu))

is kept, but the initial slope E_py = k * z of API is replaced by

    E_py = a * (z / z_ref)^b * (D / D_ref)^c * phi^d
    a = 50 MPa, b = 0.6, c = 0.5, d = 3.6, z_ref = D_ref = 1 m,
    phi in radians

Clay layers use the API Matlock curves unchanged. The springs do not depend
on the pile response (no iteration).

Units inside the curve class: kN, m, kPa.

Soil parameters (soil.csv):
    all layers  gamma_eff  effective unit weight [kN/m^3], required
    sand        phi        friction angle [deg], required -> pu and E_py
    clay        su, eps50, J as for API (api.py)
"""

import math

from .api import ApiPyModel, SandPY

A = 50.0e3    # [kPa]
B = 0.6
C = 0.5
D = 3.6
Z_REF = 1.0   # [m]
D_REF = 1.0   # [m]


class SorensenSandPY(SandPY):
    def __init__(self, layer, z, sigma_v, diameter, loading):
        super(SorensenSandPY, self).__init__(layer, z, sigma_v, diameter,
                                             loading)
        self.kz = (A * (z / Z_REF) ** B * (diameter / D_REF) ** C
                   * math.radians(layer.phi) ** D)


class SorensenPyModel(ApiPyModel):
    """API p-y springs with the Sorensen sand stiffness, no iteration."""

    name = "SORENSEN"
    label = "Sorensen"
    sand_curve = SorensenSandPY
