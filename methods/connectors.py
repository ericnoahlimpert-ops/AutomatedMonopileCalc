"""Soil springs as connector elements to ground in Abaqus/CAE.

Every spring node gets a wire from its pile vertex to ground with its own
connector section:
    components 1, 2  (U1, U2)   lateral spring curve (force vs displacement)
    components 4, 5  (UR1, UR2) rotational spring curve (moment vs rotation),
                                only for spring models with rotational springs
Component 3 (axial) and 6 (torsion) are left free. Both lateral and both
rotational directions are uncoupled and use the same curve, mirrored for
negative displacements. Base springs (pile tip) get a separate connector.

Calling create_soil_connectors again with new curves (iterative spring
models) only updates the existing connector sections.

Unlike nonlinear SPRING elements, connectors are fully supported in
Abaqus/CAE (pre-processing and Viewer).

Must be run inside Abaqus/CAE.
"""

from abaqusConstants import CARTESIAN, ROTATION, IMPRINT, OFF, NONLINEAR
from abaqusConstants import NONE as ABQ_NONE
import connectorBehavior

LATERAL_COMPONENTS = (1, 2)
ROTATIONAL_COMPONENTS = (4, 5)


def _behaviors(node_springs):
    behaviors = []
    for curve, components in ((node_springs.lateral, LATERAL_COMPONENTS),
                              (node_springs.rotational,
                               ROTATIONAL_COMPONENTS)):
        if curve is None:
            continue
        table = curve.symmetric_table()
        behaviors.extend(connectorBehavior.ConnectorElasticity(
            components=(c,), behavior=NONLINEAR, table=table)
            for c in components)
    return behaviors


def _connector(model, instance, name, node_springs, depth_to_z, rotational):
    """Create or update one connector to ground. Returns its set name, or
    None if nothing was created."""
    assembly = model.rootAssembly
    behaviors = _behaviors(node_springs)

    if name in model.sections:
        model.sections[name].setValues(behaviorOptions=behaviors)
        return name
    if not behaviors:
        return None

    model.ConnectorSection(name=name, translationalType=CARTESIAN,
                           rotationalType=ROTATION if rotational
                           else ABQ_NONE,
                           behaviorOptions=behaviors)
    vertex = instance.vertices.findAt((0.0, 0.0,
                                       depth_to_z(node_springs.z_pile)))
    wire = assembly.WirePolyLine(points=((vertex, None),),
                                 mergeType=IMPRINT, meshable=OFF)
    edges = assembly.getFeatureEdges(wire.name)
    assembly.features.changeKey(fromName=wire.name, toName="Wire-" + name)
    region = assembly.Set(name=name, edges=edges)
    assembly.SectionAssignment(region=region, sectionName=name)
    return name


def create_soil_connectors(model, instance, spring_set, depth_to_z,
                           rotational=False, set_name="SOIL_SPRINGS"):
    """Create (first call) or update (later calls) the soil connectors.

    model       Abaqus model
    instance    pile instance (every spring node must be a vertex)
    spring_set  SpringSet from a spring model
    depth_to_z  function mapping depth below pile head to global Z
    rotational  True if the spring model has rotational springs; the
                connectors then get rotational components 4 and 5

    Returns the number of connectors.
    """
    names = []
    for node in spring_set.nodes:
        name = _connector(model, instance, "PY-N%d" % node.node_label, node,
                          depth_to_z, rotational)
        if name:
            names.append(name)
    if spring_set.base is not None:
        name = _connector(model, instance,
                          "BASE-N%d" % spring_set.base.node_label,
                          spring_set.base, depth_to_z, True)
        if name:
            names.append(name)

    assembly = model.rootAssembly
    if names and set_name not in assembly.sets:
        assembly.SetByBoolean(name=set_name,
                              sets=[assembly.sets[n] for n in names])
    return len(names)
