"""Creation of the 1D monopile beam model in Abaqus/CAE.

Units: SI (m, N, kg, Pa).
Coordinate system: the pile axis is the global Z axis. The pile head sits at
Z = -z_top of the first section, depths from the CSV are mapped to Z = -z,
so the pile points downwards into the soil.

Mesh: element_count elements over the whole pile (head to tip), i.e. a
target element size of pile length / element_count. Nodes are forced at
every pile section boundary and at the given extra depths (mudline, soil
layer boundaries), so no element spans a change of pile or soil
parameters. The elements are distributed over the intervals between these
depths with lengths as uniform as possible; every interval gets at least
one element, so the count is exceeded only if there are more intervals
than element_count.

Conical sections are modelled as prismatic elements, each using the
diameter and wall thickness at its mid-depth (stepped approximation).

Every element end node is also a geometric vertex of the wire (one element
per edge), so soil connectors can be attached to vertices in the assembly.

Cross-section (pipe profile):
    diameter_reference = "outer": D from pile.csv is the outer diameter,
                                  outer radius r = D/2
    diameter_reference = "mid":   D is the mid-surface diameter, r = D/2 + t/2,
                                  which gives A = pi D t and I = pi D^3 t/8
                                  (1 + (t/D)^2), i.e. the thin-walled values
                                  used with the D of the PISA papers
    shear_factor = None:          Abaqus default transverse shear stiffness
    shear_factor = kappa:         K = kappa * G * A (e.g. 0.5 as in PISA)
PISA springs use the mid-surface diameter (D - t for "outer"), the other
spring models use D as given in pile.csv.

Must be run inside Abaqus/CAE (abaqus cae noGUI=main.py).
"""

import math

from abaqus import mdb
from abaqusConstants import (DEFORMABLE_BODY, THREE_D, IMPRINT, ON,
                             DURING_ANALYSIS, N1_COSINES, MIDDLE_SURFACE,
                             FROM_SECTION, STANDARD, FIXED, VALUE)
import abaqusConstants
import mesh
import regionToolset

TOLERANCE = 1e-4


class PileSegment(object):
    """Prismatic piece of the pile with constant cross-section."""

    def __init__(self, section_index, z_top, z_bot, diameter, thickness,
                 material):
        self.section_index = section_index
        self.z_top = z_top
        self.z_bot = z_bot
        self.diameter = diameter
        self.thickness = thickness
        self.material = material


def fixed_node_depths(sections, extra_depths=()):
    """Sorted depths below the pile head that must be nodes: section
    boundaries plus extra_depths inside the pile (mudline, soil layer
    boundaries); depths closer than TOLERANCE are merged."""
    head, tip = sections[0].z_top, sections[-1].z_bot
    depths = sorted(set([s.z_top for s in sections] + [tip]
                        + [z for z in extra_depths
                           if head + TOLERANCE < z < tip - TOLERANCE]))
    merged = [depths[0]]
    for z in depths[1:]:
        if z - merged[-1] > TOLERANCE:
            merged.append(z)
    merged[-1] = tip
    return merged


def split_into_node_depths(fixed_depths, element_count):
    """Depths of all element end nodes. fixed_depths: sorted depths that
    must be nodes, from pile head to tip. element_count elements are
    distributed over the intervals in between with element lengths as
    uniform as possible (at least one per interval)."""
    lengths = [b - a for a, b in zip(fixed_depths[:-1], fixed_depths[1:])]
    size = sum(lengths) / element_count
    counts = [max(int(math.floor(l / size + TOLERANCE)), 1) for l in lengths]
    while sum(counts) < element_count:      # refine the longest elements
        i = max(range(len(counts)), key=lambda k: lengths[k] / counts[k])
        counts[i] += 1
    while sum(counts) > element_count:      # coarsen the shortest ones
        candidates = [k for k in range(len(counts)) if counts[k] > 1]
        if not candidates:
            break                           # more intervals than elements
        i = min(candidates, key=lambda k: lengths[k] / (counts[k] - 1))
        counts[i] -= 1

    depths = [fixed_depths[0]]
    for z_top, length, n in zip(fixed_depths, lengths, counts):
        depths.extend(z_top + length * (j + 1) / n for j in range(n))
    return depths


def split_into_segments(sections, node_depths):
    """Split the CSV sections into prismatic segments. Prismatic sections
    stay one segment, conical sections get one segment per element with
    the properties at its mid-depth."""
    segments = []
    for i, s in enumerate(sections, start=1):
        if s.is_conical:
            nodes = [z for z in node_depths
                     if s.z_top - TOLERANCE <= z <= s.z_bot + TOLERANCE]
            pieces = zip(nodes[:-1], nodes[1:])
        else:
            pieces = [(s.z_top, s.z_bot)]
        for z_top, z_bot in pieces:
            z_mid = 0.5 * (z_top + z_bot)
            segments.append(PileSegment(i, z_top, z_bot, s.diameter_at(z_mid),
                                        s.thickness_at(z_mid), s.material))
    return segments


# Linear elastic properties: E [Pa], nu [-], density [kg/m^3]
MATERIAL_LIBRARY = {
    "S235": {"E": 210.0e9, "nu": 0.3, "density": 7850.0},
    "S275": {"E": 210.0e9, "nu": 0.3, "density": 7850.0},
    "S355": {"E": 210.0e9, "nu": 0.3, "density": 7850.0},
    "S420": {"E": 210.0e9, "nu": 0.3, "density": 7850.0},
    "S460": {"E": 210.0e9, "nu": 0.3, "density": 7850.0},
    # Pile steel of the PISA calibration analyses (Burd et al. 2020);
    # density not given in the paper, standard steel value assumed
    "PISA_STEEL": {"E": 200.0e9, "nu": 0.3, "density": 7850.0},
}


class PileGeometryBuilder(object):
    """Builds part, materials, sections, mesh and instance for a 1D pile."""

    def __init__(self, sections, model_name="Monopile", part_name="Pile",
                 element_count=60, element_type="B31",
                 diameter_reference="outer", shear_factor=None,
                 fixed_depths=()):
        """fixed_depths: extra depths below the pile head that must be
        nodes (mudline, soil layer boundaries)."""
        if diameter_reference not in ("outer", "mid"):
            raise ValueError("diameter_reference must be 'outer' or 'mid', "
                             "got '%s'." % diameter_reference)
        if shear_factor is not None and shear_factor <= 0.0:
            raise ValueError("shear_factor must be positive or None.")
        if isinstance(element_count, bool) \
                or not isinstance(element_count, int) or element_count < 1:
            raise ValueError("ELEMENT_COUNT must be a positive integer, got "
                             "%r." % (element_count,))
        self.sections = sections
        self.model_name = model_name
        self.part_name = part_name
        self.element_type = element_type
        self.diameter_reference = diameter_reference
        self.shear_factor = shear_factor
        self.node_depths = split_into_node_depths(
            fixed_node_depths(sections, fixed_depths), element_count)
        self.segments = split_into_segments(sections, self.node_depths)
        self.model = None
        self.part = None
        self.instance = None

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def depth_to_z(depth):
        return -float(depth)

    @property
    def z_head(self):
        return self.sections[0].z_top

    @property
    def z_tip(self):
        return self.sections[-1].z_bot

    def _edges_between(self, z_top, z_bot):
        """All part edges lying on the pile axis between two depths."""
        return self.part.edges.getByBoundingBox(
            xMin=-TOLERANCE, yMin=-TOLERANCE,
            zMin=self.depth_to_z(z_bot) - TOLERANCE,
            xMax=TOLERANCE, yMax=TOLERANCE,
            zMax=self.depth_to_z(z_top) + TOLERANCE)

    def outer_radius(self, diameter, thickness):
        """Outer radius of the pipe profile for the diameter reference."""
        if self.diameter_reference == "mid":
            return 0.5 * diameter + 0.5 * thickness
        return 0.5 * diameter

    def _vertex_at(self, depth):
        return self.part.vertices.findAt(((0.0, 0.0, self.depth_to_z(depth)),))

    def nodes(self):
        """(node label, depth below pile head) of every element end node.
        Call after meshing."""
        result = []
        for depth in self.node_depths:
            vertex = self.part.vertices.findAt(
                (0.0, 0.0, self.depth_to_z(depth)))
            result.append((vertex.getNodes()[0].label, depth))
        return result

    # ------------------------------------------------------------------ #
    # Build steps
    # ------------------------------------------------------------------ #
    def build(self):
        self.create_model()
        self.create_part()
        self.create_materials()
        self.create_sections()
        self.create_mesh()
        self.create_instance()
        return self.model

    def create_model(self):
        # Abaqus cannot delete the last model of a database, so the new model
        # is created first and renamed after the old one has been removed.
        temp_name = self.model_name + "-new"
        self.model = mdb.Model(name=temp_name)
        for name in list(mdb.models.keys()):
            if name != temp_name and name in (self.model_name, "Model-1"):
                for job_name in [j for j in mdb.jobs.keys()
                                 if mdb.jobs[j].model == name]:
                    del mdb.jobs[job_name]
                del mdb.models[name]
        mdb.models.changeKey(fromName=temp_name, toName=self.model_name)
        self.model = mdb.models[self.model_name]

    def create_part(self):
        """Wire part along Z with one vertex at every element end node."""
        points = tuple((0.0, 0.0, self.depth_to_z(d))
                       for d in self.node_depths)

        self.part = self.model.Part(name=self.part_name,
                                    dimensionality=THREE_D,
                                    type=DEFORMABLE_BODY)
        self.part.WirePolyLine(points=points, mergeType=IMPRINT, meshable=ON)

        self.part.Set(name="PILE", edges=self._edges_between(self.z_head,
                                                             self.z_tip))
        self.part.Set(name="PILE_HEAD", vertices=self._vertex_at(self.z_head))
        self.part.Set(name="PILE_TIP", vertices=self._vertex_at(self.z_tip))

    def create_materials(self):
        for name in sorted(set(s.material for s in self.sections)):
            key = name.upper()
            if key not in MATERIAL_LIBRARY:
                raise KeyError("Material '%s' is not in MATERIAL_LIBRARY "
                               "(available: %s)."
                               % (name, ", ".join(sorted(MATERIAL_LIBRARY))))
            props = MATERIAL_LIBRARY[key]
            material = self.model.Material(name=name)
            material.Elastic(table=((props["E"], props["nu"]),))
            material.Density(table=((props["density"],),))

    def create_sections(self):
        """One edge set SEC-i per CSV row, plus one pipe profile + beam section
        per prismatic segment (several per row for conical sections)."""
        for i, s in enumerate(self.sections, start=1):
            self.part.Set(name="SEC-%d" % i,
                          edges=self._edges_between(s.z_top, s.z_bot))

        for i, s in enumerate(self.sections, start=1):
            segments = [seg for seg in self.segments if seg.section_index == i]
            for j, seg in enumerate(segments, start=1):
                suffix = "%d" % i if len(segments) == 1 else "%d-%d" % (i, j)
                profile_name = "PipeProfile-" + suffix
                section_name = "PipeSection-" + suffix

                r_outer = self.outer_radius(seg.diameter, seg.thickness)
                self.model.PipeProfile(name=profile_name, r=r_outer,
                                       t=seg.thickness)
                props = MATERIAL_LIBRARY[seg.material.upper()]
                section = self.model.BeamSection(
                    name=section_name, integration=DURING_ANALYSIS,
                    profile=profile_name, material=seg.material,
                    poissonRatio=props["nu"])
                if self.shear_factor is not None:
                    t = seg.thickness
                    area = math.pi * (r_outer ** 2 - (r_outer - t) ** 2)
                    shear_modulus = props["E"] / (2.0 * (1.0 + props["nu"]))
                    k = self.shear_factor * shear_modulus * area
                    section.TransverseShearBeam(scfDefinition=VALUE,
                                                k23=k, k13=k)

                region = regionToolset.Region(
                    edges=self._edges_between(seg.z_top, seg.z_bot))
                self.part.SectionAssignment(region=region,
                                            sectionName=section_name,
                                            offsetType=MIDDLE_SURFACE,
                                            thicknessAssignment=FROM_SECTION)

        # Beam orientation: n1 must not be parallel to the pile axis (Z)
        self.part.assignBeamSectionOrientation(region=self.part.sets["PILE"],
                                               method=N1_COSINES,
                                               n1=(1.0, 0.0, 0.0))

    def create_mesh(self):
        # Edges already have element length: exactly one element per edge
        self.part.seedEdgeByNumber(edges=self.part.edges, number=1,
                                   constraint=FIXED)
        elem_code = getattr(abaqusConstants, self.element_type.upper())
        elem_type = mesh.ElemType(elemCode=elem_code, elemLibrary=STANDARD)
        region = regionToolset.Region(edges=self.part.edges)
        self.part.setElementType(regions=region, elemTypes=(elem_type,))
        self.part.generateMesh()

    def create_instance(self):
        assembly = self.model.rootAssembly
        self.instance = assembly.Instance(name=self.part_name + "-1",
                                          part=self.part, dependent=ON)
