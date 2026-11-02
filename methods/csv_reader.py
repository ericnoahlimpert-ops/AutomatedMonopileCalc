"""Reading and validation of the input CSV files.

The first line of both files is a header.

pile.csv
    z_top, z_bot, diameter, thickness, material

    Columns are found by their header names (case-insensitive, surrounding
    spaces ignored, any order).

    All lengths in metres. z is the depth below the pile head (positive
    downwards), so z_bot > z_top for every section.
    Conical sections: diameter and thickness may be given as "top-bottom",
    e.g. "1-1.2" or "0.008-0.01". Values vary linearly from z_top to z_bot.

soil.csv
    z_top, z_bot, soil_type, gamma_eff [, phi, su, eps50, J, k, dr, g0,
                                          nu, es_ref, lambda_es, g0_ref,
                                          lambda_g0, hb_b, e0, pisa_clay]

    Columns are assigned by number in Input.py (SOIL_COLUMNS); the header
    names do not matter.

    z            depth below mudline [m]
    soil_type    "sand" or "clay"
    gamma_eff    effective (submerged) unit weight [kN/m^3]
    phi          sand: friction angle [deg]
    k            sand: initial modulus of subgrade reaction [MN/m^3]
                 (optional, default from the API chart, see spring_models/api.py)
    su           clay: undrained shear strength [kPa], "top-bottom" allowed
    eps50        clay: strain at 50 % of max. stress [-] (optional, 0.01)
    J            clay: empirical Matlock constant [-] (optional, 0.5)
    dr           sand: relative density [%] (PISA)
    g0           small-strain shear modulus [MPa], "top-bottom" allowed
                 (optional; if given, it replaces the G0 calculation)
    nu           sand: Poisson's ratio [-]                       (Thieken)
    es_ref       sand: oedometric stiffness Es at the reference
                 mean stress 100 kPa [MPa]                        (Thieken)
    lambda_es    sand: stress exponent of Es [-]                  (Thieken)
    g0_ref       sand: small-strain shear modulus G0 at the
                 reference mean stress 100 kPa [MPa]   (G0 for POWER_LAW)
    lambda_g0    sand: stress exponent of G0 [-]       (G0 for POWER_LAW)
                 Es = es_ref (s_m/100 kPa)^lambda_es,
                 G0 = g0_ref (s_m/100 kPa)^lambda_g0.
    hb_b         sand: soil specific constant B of Hardin & Black [-]
                 (Dunkirk sand: 875)                (G0 for HARDIN_BLACK)
    e0           sand: void ratio [-]               (G0 for HARDIN_BLACK)
                 G0 = B 101.3 kPa / (0.3 + 0.7 e0^2) sqrt(s_m/101.3 kPa).
                 Mean stress s_m = s'v (1 + 2 k0)/3 with k0 = 1 - sin phi,
                 or k0 = 0.4 if phi is not given.
                 G0 is chosen per layer from the given columns: g0, else
                 g0_ref + lambda_g0 (power law), else hb_b + e0 (Hardin &
                 Black), see spring_models/base.py.
                 The "Soil parameters" section of each spring model lists
                 the columns it needs.

    pisa_clay    clay: PISA clay calibration, "cowden" (Cowden till,
                 default) or "bothkennar" (Bothkennar clay)
                 (PISA, PISA_CLAY)

    Only z_top, z_bot, soil_type and gamma_eff are required columns.
    Columns that do not apply to a layer or spring model are left empty.

cpt.csv (TOTT_BUSWELL only)
    z, qc [, g0]

    Columns are found by their header names (case-insensitive, any order,
    further columns ignored).

    z            depth below mudline [m]
    qc           cone resistance [MPa]
    g0           small-strain shear modulus [MPa] (optional; replaces the
                 G0 of soil.csv)
"""

import csv
import re

_NUMBER = r"\d+(?:\.\d*)?(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?"
_RANGE = re.compile(r"^\s*(%s)\s*(?:-\s*(%s))?\s*$" % (_NUMBER, _NUMBER))

REQUIRED_COLUMNS = ("z_top", "z_bot", "diameter", "thickness", "material")
SOIL_REQUIRED_COLUMNS = ("z_top", "z_bot", "soil_type", "gamma_eff")
SOIL_OPTIONAL_COLUMNS = ("phi", "su", "eps50", "j", "k", "dr", "g0", "nu",
                         "es_ref", "lambda_es", "g0_ref", "lambda_g0",
                         "hb_b", "e0", "pisa_clay")
SOIL_PARAMETERS = SOIL_REQUIRED_COLUMNS + SOIL_OPTIONAL_COLUMNS
SOIL_TYPES = ("sand", "clay")
DEFAULT_EPS50 = 0.01
DEFAULT_J = 0.5
TOLERANCE = 1e-6


def parse_range(value, name):
    """Parse "1.2" -> (1.2, 1.2) or "1-1.2" -> (1.0, 1.2)."""
    match = _RANGE.match(str(value))
    if not match:
        raise ValueError("invalid %s '%s' (expected a number or 'top-bottom', "
                         "e.g. '1-1.2')" % (name, value))
    top = float(match.group(1))
    bot = float(match.group(2)) if match.group(2) is not None else top
    return top, bot


class PileSection(object):
    """One tubular can of the pile. Diameter and wall thickness vary linearly
    from z_top to z_bot (constant if top and bottom values are equal)."""

    def __init__(self, z_top, z_bot, diameter, thickness, material):
        self.z_top = float(z_top)
        self.z_bot = float(z_bot)
        self.diameter_top, self.diameter_bot = parse_range(diameter,
                                                           "diameter")
        self.thickness_top, self.thickness_bot = parse_range(thickness,
                                                             "thickness")
        self.material = str(material).strip()

    @property
    def length(self):
        return self.z_bot - self.z_top

    @property
    def is_conical(self):
        return (abs(self.diameter_top - self.diameter_bot) > TOLERANCE
                or abs(self.thickness_top - self.thickness_bot) > TOLERANCE)

    def _interpolate(self, top, bot, z):
        xi = (z - self.z_top) / self.length
        return top + xi * (bot - top)

    def diameter_at(self, z):
        return self._interpolate(self.diameter_top, self.diameter_bot, z)

    def thickness_at(self, z):
        return self._interpolate(self.thickness_top, self.thickness_bot, z)

    def __repr__(self):
        return ("PileSection(z_top=%g, z_bot=%g, D=%g-%g, t=%g-%g, "
                "material=%s)"
                % (self.z_top, self.z_bot, self.diameter_top,
                   self.diameter_bot, self.thickness_top, self.thickness_bot,
                   self.material))


class SoilLayer(object):
    """One soil layer for the API p-y formulation. Depths below mudline."""

    def __init__(self, z_top, z_bot, soil_type, gamma_eff, phi=None, su=None,
                 eps50=None, j=None, k=None, dr=None, g0=None, nu=None,
                 es_ref=None, lambda_es=None, g0_ref=None, lambda_g0=None,
                 hb_b=None, e0=None, pisa_clay=None):
        self.z_top = float(z_top)
        self.z_bot = float(z_bot)
        self.soil_type = str(soil_type).strip().lower()
        if self.soil_type not in SOIL_TYPES:
            raise ValueError("soil_type must be one of %s, got '%s'"
                             % (", ".join(SOIL_TYPES), soil_type))
        self.gamma_eff = float(gamma_eff)
        self.phi = _optional_float(phi)
        self.k = _optional_float(k)
        if _is_blank(su):
            self.su_top = self.su_bot = None
        else:
            self.su_top, self.su_bot = parse_range(su, "su")
        self.eps50 = _optional_float(eps50, DEFAULT_EPS50)
        self.j = _optional_float(j, DEFAULT_J)
        self.dr = _optional_float(dr)  # [%]
        if _is_blank(g0):
            self.g0_top = self.g0_bot = None
        else:
            self.g0_top, self.g0_bot = parse_range(g0, "g0")  # [MPa]
        self.nu = _optional_float(nu)
        self.es_ref = _optional_float(es_ref)  # [MPa] at 100 kPa
        self.lambda_es = _optional_float(lambda_es)
        self.g0_ref = _optional_float(g0_ref)  # [MPa] at 100 kPa
        self.lambda_g0 = _optional_float(lambda_g0)
        self.hb_b = _optional_float(hb_b)      # Hardin & Black B [-]
        self.e0 = _optional_float(e0)          # void ratio [-]
        # PISA clay calibration (text, e.g. "cowden"), checked by the model
        self.pisa_clay = None if _is_blank(pisa_clay)             else str(pisa_clay).strip().lower()

    @property
    def length(self):
        return self.z_bot - self.z_top

    def _interpolate(self, top, bot, z):
        xi = (z - self.z_top) / self.length
        return top + xi * (bot - top)

    def su_at(self, z):
        return self._interpolate(self.su_top, self.su_bot, z)

    def g0_at(self, z):
        """Small-strain shear modulus [kPa] at depth z, or None if not given."""
        if self.g0_top is None:
            return None
        return 1000.0 * self._interpolate(self.g0_top, self.g0_bot, z)

    def __repr__(self):
        if self.soil_type == "sand":
            # only the parameters given in soil.csv
            given = [(name, getattr(self, name)) for name in
                     ("phi", "k", "dr", "nu", "es_ref", "lambda_es",
                      "g0_ref", "lambda_g0", "hb_b", "e0")]
            props = ", ".join("%s=%g" % (name, value)
                              for name, value in given if value is not None)
        else:
            props = "su=%g-%g, eps50=%g, J=%g" % (self.su_top, self.su_bot,
                                                 self.eps50, self.j)
            if self.pisa_clay is not None:
                props += ", PISA clay=%s" % self.pisa_clay
        if self.g0_top is not None:
            props += ", G0=%g-%g MPa" % (self.g0_top, self.g0_bot)
        return ("SoilLayer(z_top=%g, z_bot=%g, %s, gamma'=%g, %s)"
                % (self.z_top, self.z_bot, self.soil_type, self.gamma_eff,
                   props))


def _is_blank(value):
    return value is None or str(value).strip() == ""


def _optional_float(value, default=None):
    return default if _is_blank(value) else float(value)


def _normalise(name):
    """Column names are compared case-insensitively, without spaces."""
    return str(name).strip().lower()


def _read_header(path, label):
    with open(path, newline="") as f:
        header = next(csv.reader(f, skipinitialspace=True), None)
    if not header:
        raise ValueError("%s CSV '%s' is empty." % (label, path))
    return [_normalise(name) for name in header]


def _columns_by_name(header, parameters, label, path):
    """{parameter: column index} from the header names (pile.csv)."""
    missing = [p for p in parameters if p not in header]
    if missing:
        raise ValueError("%s CSV '%s' is missing column(s): %s"
                         % (label, path, ", ".join(missing)))
    return dict((p, header.index(p)) for p in parameters)


def _columns_by_number(header, columns, parameters, required, label, path):
    """{parameter: column index (0-based) or None} from column numbers.

    columns: {parameter: column number, counted from 1, or None}, e.g.
    SOIL_COLUMNS from Input.py. Parameters that are missing or None are
    not contained in the CSV."""
    resolved, errors, used = {}, [], {}
    for parameter in parameters:
        number = columns.get(parameter)
        if number is None:
            resolved[parameter] = None
            continue
        if isinstance(number, bool) or not isinstance(number, int)                 or not 1 <= number <= len(header):
            errors.append("%s = %r is not a column number between 1 and %d"
                          % (parameter, number, len(header)))
            continue
        if number in used:
            errors.append("%s and %s both use column %d"
                          % (used[number], parameter, number))
        used[number] = parameter
        resolved[parameter] = number - 1

    errors += ["%s has no column number" % p for p in required
               if p in resolved and resolved[p] is None]
    if errors:
        raise ValueError("%s CSV '%s', column numbers in Input.py: %s. "
                         "Columns in the file: %s."
                         % (label, path, "; ".join(errors),
                            ", ".join("%d = %s" % (i + 1, name)
                                      for i, name in enumerate(header))))
    return resolved


def _read_rows(path, indices):
    """Yield (line number, {parameter: value}) for every non-blank data row.
    indices: {parameter: column index or None}. The first line is the
    header and is skipped."""
    with open(path, newline="") as f:
        reader = csv.reader(f, skipinitialspace=True)
        next(reader, None)
        # Data starts on line 2 (line 1 is the header)
        for line_no, row in enumerate(reader, start=2):
            if not any(value.strip() for value in row):
                continue  # skip blank lines
            yield line_no, dict(
                (parameter, row[i] if i is not None and i < len(row)
                 else None)
                for parameter, i in indices.items())


def read_pile_csv(path):
    """Read the pile CSV and return a validated list of PileSection objects,
    sorted from pile head to pile tip."""
    sections = []
    indices = _columns_by_name(_read_header(path, "Pile"), REQUIRED_COLUMNS,
                               "Pile", path)
    for line_no, row in _read_rows(path, indices):
        try:
            section = PileSection(*(row[c] for c in REQUIRED_COLUMNS))
        except (TypeError, ValueError) as err:
            raise ValueError("Pile CSV '%s', line %d: %s"
                             % (path, line_no, err))
        sections.append(section)

    sections.sort(key=lambda s: s.z_top)
    validate_sections(sections)
    return sections


def read_soil_csv(path, columns):
    """Read the soil CSV and return a validated list of SoilLayer objects,
    sorted from mudline downwards.

    columns: {parameter: column number (from 1) or None}, i.e. SOIL_COLUMNS
    from Input.py. The first line of the file is a header and is skipped."""
    layers = []
    indices = _columns_by_number(_read_header(path, "Soil"), columns,
                                 SOIL_PARAMETERS, SOIL_REQUIRED_COLUMNS,
                                 "Soil", path)
    for line_no, row in _read_rows(path, indices):
        try:
            layer = SoilLayer(**row)
        except (TypeError, ValueError) as err:
            raise ValueError("Soil CSV '%s', line %d: %s"
                             % (path, line_no, err))
        layers.append(layer)

    layers.sort(key=lambda l: l.z_top)
    validate_soil_layers(layers)
    return layers


def unused_soil_columns(path, columns):
    """Columns of the soil CSV that are not assigned to any parameter in
    Input.py, as "number (header)" strings."""
    header = _read_header(path, "Soil")
    used = set(n for p, n in columns.items() if p in SOIL_PARAMETERS)
    return ["%d (%s)" % (i + 1, name) for i, name in enumerate(header)
            if i + 1 not in used and name]


def validate_sections(sections):
    """Check that the sections form one continuous, physically valid pile."""
    if not sections:
        raise ValueError("No pile sections defined.")

    for i, s in enumerate(sections, start=1):
        if s.length <= TOLERANCE:
            raise ValueError("Section %d: z_bot (%g) must be greater than "
                             "z_top (%g)." % (i, s.z_bot, s.z_top))
        for end, d, t in (("top", s.diameter_top, s.thickness_top),
                          ("bottom", s.diameter_bot, s.thickness_bot)):
            if d <= 0.0:
                raise ValueError("Section %d (%s): diameter must be positive."
                                 % (i, end))
            if not 0.0 < t < 0.5 * d:
                raise ValueError("Section %d (%s): thickness must be between "
                                 "0 and D/2." % (i, end))
        if not s.material:
            raise ValueError("Section %d: no material given." % i)

    for i in range(1, len(sections)):
        above, below = sections[i - 1], sections[i]
        if abs(below.z_top - above.z_bot) > TOLERANCE:
            raise ValueError("Sections %d and %d are not continuous "
                             "(z_bot=%g, next z_top=%g)."
                             % (i, i + 1, above.z_bot, below.z_top))


def validate_soil_layers(layers):
    """Check that the layers start at the mudline, are continuous and have
    the parameters their soil type needs."""
    if not layers:
        raise ValueError("No soil layers defined.")
    if abs(layers[0].z_top) > TOLERANCE:
        raise ValueError("The first soil layer must start at the mudline "
                         "(z_top = 0), got %g." % layers[0].z_top)

    for i, l in enumerate(layers, start=1):
        if l.length <= TOLERANCE:
            raise ValueError("Soil layer %d: z_bot (%g) must be greater than "
                             "z_top (%g)." % (i, l.z_bot, l.z_top))
        if l.gamma_eff <= 0.0:
            raise ValueError("Soil layer %d: gamma_eff must be positive." % i)
        if l.soil_type == "sand":
            # phi is required by API and THIEKEN (checked in the models)
            if l.phi is not None and not 0.0 < l.phi < 90.0:
                raise ValueError("Soil layer %d (sand): phi must be between 0 "
                                 "and 90 deg." % i)
            if l.k is not None and l.k <= 0.0:
                raise ValueError("Soil layer %d (sand): k must be positive."
                                 % i)
        else:
            if l.su_top is None or min(l.su_top, l.su_bot) <= 0.0:
                raise ValueError("Soil layer %d (clay): su must be positive."
                                 % i)
            if l.eps50 <= 0.0 or l.j <= 0.0:
                raise ValueError("Soil layer %d (clay): eps50 and J must be "
                                 "positive." % i)
        if l.dr is not None and not 1.0 < l.dr <= 100.0:
            raise ValueError("Soil layer %d: dr is the relative density in %% "
                             "(e.g. 75), got %g." % (i, l.dr))
        if l.g0_top is not None and min(l.g0_top, l.g0_bot) <= 0.0:
            raise ValueError("Soil layer %d: g0 must be positive." % i)
        if l.nu is not None and not 0.0 <= l.nu < 0.5:
            raise ValueError("Soil layer %d: nu must be in [0, 0.5)." % i)
        for name in ("es_ref", "lambda_es", "g0_ref", "lambda_g0", "hb_b",
                     "e0"):
            value = getattr(l, name)
            if value is not None and value <= 0.0:
                raise ValueError("Soil layer %d: %s must be positive."
                                 % (i, name))

    for i in range(1, len(layers)):
        above, below = layers[i - 1], layers[i]
        if abs(below.z_top - above.z_bot) > TOLERANCE:
            raise ValueError("Soil layers %d and %d are not continuous "
                             "(z_bot=%g, next z_top=%g)."
                             % (i, i + 1, above.z_bot, below.z_top))


# --------------------------------------------------------------------------- #
# CPT profile (TOTT_BUSWELL)
# --------------------------------------------------------------------------- #
CPT_REQUIRED_COLUMNS = ("z", "qc")


class CptProfile(object):
    """CPT profile below the mudline: cone resistance qc and, optionally,
    the small-strain shear modulus G0 (e.g. from a seismic CPT). Values are
    interpolated linearly between the data points and kept constant above
    the first and below the last point. Depths in m, stresses in kPa."""

    def __init__(self, depths, qc, g0=None, path=""):
        self.depths = list(depths)
        self.qc = list(qc)
        self.g0 = list(g0) if g0 is not None else None
        self.path = path

    @property
    def has_g0(self):
        return self.g0 is not None

    @property
    def z_max(self):
        return self.depths[-1]

    def _values(self, name):
        values = getattr(self, name)
        if values is None:
            raise ValueError("CPT CSV '%s' has no %s column."
                             % (self.path, name))
        return values

    def value_at(self, z, name="qc"):
        """qc or g0 [kPa] at depth z below mudline."""
        values, depths = self._values(name), self.depths
        if z <= depths[0]:
            return values[0]
        if z >= depths[-1]:
            return values[-1]
        for i in range(1, len(depths)):
            if z <= depths[i]:
                xi = (z - depths[i - 1]) / (depths[i] - depths[i - 1])
                return values[i - 1] + xi * (values[i] - values[i - 1])
        return values[-1]

    def average(self, z0, z1, name="qc"):
        """Mean of qc or g0 [kPa] between the depths z0 < z1 (exact integral
        of the piecewise linear profile)."""
        if z1 - z0 <= TOLERANCE:
            return self.value_at(0.5 * (z0 + z1), name)
        points = [z0] + [z for z in self.depths if z0 < z < z1] + [z1]
        total = 0.0
        for a, b in zip(points[:-1], points[1:]):
            total += 0.5 * (self.value_at(a, name)
                            + self.value_at(b, name)) * (b - a)
        return total / (z1 - z0)


def read_cpt_csv(path):
    """Read a CPT CSV with the header columns z [m below mudline],
    qc [MPa] and, optionally, g0 [MPa] (names case-insensitive, any order,
    further columns ignored). Returns a CptProfile with values in kPa."""
    header = _read_header(path, "CPT")
    names = CPT_REQUIRED_COLUMNS + (("g0",) if "g0" in header else ())
    indices = _columns_by_name(header, names, "CPT", path)
    rows = []
    for line_no, row in _read_rows(path, indices):
        try:
            values = dict((name, float(row[name])) for name in names)
        except (TypeError, ValueError):
            raise ValueError("CPT CSV '%s', line %d: z, qc%s must be numbers."
                             % (path, line_no,
                                ", g0" if "g0" in names else ""))
        if values["z"] < 0.0:
            raise ValueError("CPT CSV '%s', line %d: z must not be negative "
                             "(depth below mudline)." % (path, line_no))
        if values["qc"] <= 0.0 or values.get("g0", 1.0) <= 0.0:
            raise ValueError("CPT CSV '%s', line %d: qc and g0 must be "
                             "positive." % (path, line_no))
        rows.append(values)
    if len(rows) < 2:
        raise ValueError("CPT CSV '%s' needs at least two data rows." % path)
    rows.sort(key=lambda r: r["z"])
    for a, b in zip(rows[:-1], rows[1:]):
        if b["z"] - a["z"] <= TOLERANCE:
            raise ValueError("CPT CSV '%s': depth %g appears twice."
                             % (path, a["z"]))
    return CptProfile([r["z"] for r in rows],
                      [1000.0 * r["qc"] for r in rows],
                      [1000.0 * r["g0"] for r in rows]
                      if "g0" in names else None, path)
