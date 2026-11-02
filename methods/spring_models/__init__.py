"""Soil spring models. Register new models in SPRING_MODELS."""

from .api import ApiPyModel
from .base import SpringModel, SpringSet, NodeSprings, write_springs_csv
from .tott_buswell import TottBuswellModel
from .kallehave import KallehavePyModel
from .pisa import PisaSandModel
from .pisa_clay import PisaClayModel
from .pisa_layered import PisaLayeredModel
from .sorensen import SorensenPyModel
from .thieken import ThiekenSandModel

SPRING_MODELS = {
    ApiPyModel.name: ApiPyModel,
    KallehavePyModel.name: KallehavePyModel,
    SorensenPyModel.name: SorensenPyModel,
    PisaSandModel.name: PisaSandModel,
    PisaClayModel.name: PisaClayModel,
    PisaLayeredModel.name: PisaLayeredModel,
    ThiekenSandModel.name: ThiekenSandModel,
    TottBuswellModel.name: TottBuswellModel,
}


def create_spring_model(name, pile_sections, soil_layers, mudline_depth,
                        curve_points=None, grid_start=None, **options):
    """Create a spring model by name. Model specific options (e.g. loading
    for API, load_height for PISA, cpt_profile for TOTT_BUSWELL) are passed
    through; unused ones are ignored by the models. curve_points and
    grid_start set the sampling of the spring curves (None: defaults of
    base.py)."""
    key = name.upper()
    if key not in SPRING_MODELS:
        raise KeyError("Unknown spring model '%s' (available: %s)."
                       % (name, ", ".join(sorted(SPRING_MODELS))))
    if curve_points is not None and (isinstance(curve_points, bool)
                                     or not isinstance(curve_points, int)
                                     or curve_points < 3):
        raise ValueError("CURVE_POINTS must be an integer >= 3, got %r."
                         % (curve_points,))
    if grid_start is not None and not 0.0 < grid_start < 1.0:
        raise ValueError("GRID_START must be between 0 and 1, got %r."
                         % (grid_start,))
    model = SPRING_MODELS[key](pile_sections, soil_layers, mudline_depth,
                               **options)
    if curve_points is not None:
        model.curve_points = curve_points
    if grid_start is not None:
        model.grid_start = grid_start
    return model
