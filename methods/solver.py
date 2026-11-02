"""Incremental loading with iterative soil spring models.

The load is applied in levels, either force controlled (H_k = k * dH) or
displacement controlled (target mudline displacement v_k = k * dv; H is
found for every level). At every level the analysis is repeated until the
spring model has converged (iterative models such as PISA and THIEKEN;
non-iterative models need one run per level under force control):

    for each load level k = 1 .. n:
        repeat:
            springs = spring_model.compute(nodes, state)
            update connectors and load H, run job, read node results
            residual = max |results - state| / max |results|  (v and psi,
                       all pile nodes)
            state = relaxed node results
            displacement control: correct H (secant), state scaled to the
                                  target displacement
        until residual < tolerance (and mudline displacement = target)
        check the break criteria

The residual compares the result of a run with the state its springs were
built from (fixed point), so a well predicted level can converge after a
single run.

Each load level is an Abaqus run loaded from zero to H_k (the soil springs
are nonlinear elastic, so the result does not depend on the load path). The
starting state of a level is predicted from the converged levels before
(force control: v ~ H^a, with a from the last two levels; displacement
control: last deflection line scaled to the new target). If a force
controlled level fails from the predicted state it is retried once from the
cold start of the spring model (e.g. PISA with p = p_u); if that fails too,
the incremental loading stops.

Displacement control: H of a level is corrected by H + K (v_target - v).
K is the secant of the last two runs if it is plausible (positive, within a
factor 4 of the stiffness dH/dv between the last two levels), otherwise
that reference stiffness; on the concave load-displacement curve it is
larger than the current tangent, so the correction does not overshoot even
when the springs change between the runs. H is kept inside the bracket of
loads known to give too small / too large displacements (bisection
otherwise); a failed job (e.g. H beyond the capacity) counts as "too
large". For iterative models the springs change from run to run, so only
the bound from failed jobs is kept (an older bound may exclude the
converged solution). The first H of a level
is extrapolated from the last two levels (first level: HORIZONTAL_LOAD *
level factor). The moment stays M = H * lever arm. Non-iterative spring
models compute their springs once and only correct H.

Break criteria (BreakCriteria, all at the mudline): after every level the
loading stops if the mudline displacement, mudline rotation or the toe
rotation ratio |psi_toe| / |psi_mudline| (vertical tangent of the
deflection line at the toe) reaches its limit. The level is kept, and the
load at the criterion is interpolated linearly between the last two levels
(break_info).

All jobs run in a temporary work directory. Only the output database of the
last converged load level is kept (copied to final_odb); afterwards the
model is reset to the springs and load of that level, so the saved .cae
matches the kept .odb.

Must be run inside Abaqus/CAE.
"""

import math
import shutil
import time

import os

from .analysis import (add_horizontal_load, read_head_results,
                       read_head_results_all_frames, read_node_results,
                       read_node_results_all_frames, run_job,
                       set_output_time_points)
from .connectors import create_soil_connectors
from .progress import duration

FORCE = "force"
DISPLACEMENT = "displacement"


class IterationSettings(object):
    def __init__(self, max_iterations=20, tolerance=1e-3, relaxation=0.5):
        self.max_iterations = max_iterations
        self.tolerance = tolerance      # max. relative fixed-point residual
        self.relaxation = relaxation    # share of the new displacements


class LoadLevelResult(object):
    def __init__(self, increment, load_factor, load, moment, head_u, head_ur,
                 iterations, converged, mudline=(None, None),
                 toe=(None, None)):
        self.increment = increment
        self.load_factor = load_factor  # fraction of the load / displacement
        self.load = load                # [N]
        self.moment = moment            # [Nm] at the pile head
        self.head_u = head_u            # (U1, U2, U3) at pile head [m]
        self.head_ur = head_ur          # (UR1, UR2, UR3) at pile head [rad]
        self.iterations = iterations
        self.converged = converged
        self.mudline = mudline          # (v [m], psi [rad]) at the mudline
        self.toe = toe                  # (v [m], psi [rad]) at the pile toe


class BreakInfo(object):
    """A break criterion that stopped the loading."""

    def __init__(self, name, limit, unit, load, displacement):
        self.name = name                # e.g. "max. rotation"
        self.limit = limit              # in the unit of the criterion
        self.unit = unit
        self.load = load                # interpolated H at the criterion [N]
        self.displacement = displacement  # mudline v at the criterion [m]

    @property
    def limit_text(self):
        return ("%g %s" % (self.limit, self.unit)).strip()

    def __str__(self):
        text = "%s (%s) reached" % (self.name, self.limit_text)
        if self.load is None:
            return text + " already at the first load level"
        return text + " at H = %.4g MN (mudline displacement %.2f mm)" % (
            self.load / 1e6, 1000.0 * self.displacement)


class BreakCriteria(object):
    """Break criteria of the incremental loading, evaluated at the mudline
    after every load level (None = not used):

        max_displacement        |v_mudline| >= limit [m]
        max_rotation            |psi_mudline| >= limit [deg]
        vertical_tangent_ratio  |psi_toe| / |psi_mudline| >= limit [-]: the
                                deflection line no longer has a vertical
                                tangent at the pile toe
    """

    def __init__(self, max_displacement=None, max_rotation=None,
                 vertical_tangent_ratio=None):
        # (name, limit, unit, value function, zero at H = 0)
        self.criteria = []
        if max_displacement is not None:
            self.criteria.append(("max. mudline displacement",
                                  max_displacement, "m",
                                  lambda r: abs(r.mudline[0]), True))
        if max_rotation is not None:
            self.criteria.append(("max. mudline rotation", max_rotation,
                                  "deg",
                                  lambda r: math.degrees(abs(r.mudline[1])),
                                  True))
        if vertical_tangent_ratio is not None:
            self.criteria.append(("vertical tangent at the toe "
                                  "(|psi_toe| / |psi_mudline|)",
                                  vertical_tangent_ratio, "",
                                  lambda r: abs(r.toe[1])
                                  / (abs(r.mudline[1]) or 1e-30), False))
        for name, limit, _, _, _ in self.criteria:
            if limit <= 0.0:
                raise ValueError("Break criterion %s: limit must be "
                                 "positive, got %g." % (name, limit))

    def __bool__(self):
        return bool(self.criteria)

    __nonzero__ = __bool__

    def check(self, results):
        """BreakInfo if the last result reaches a criterion, else None. The
        load at the criterion is interpolated between the last two levels
        (from H = 0 for criteria that are 0 there)."""
        last = results[-1]
        for name, limit, unit, value, zero_at_origin in self.criteria:
            if value(last) < limit:
                continue
            if len(results) > 1:
                previous = results[-2]
                v0, h0, d0 = (value(previous), previous.load,
                              abs(previous.mudline[0]))
            elif zero_at_origin:
                v0, h0, d0 = 0.0, 0.0, 0.0
            else:
                return BreakInfo(name, limit, unit, None, None)
            v1, h1, d1 = value(last), last.load, abs(last.mudline[0])
            w = (limit - v0) / (v1 - v0) if v1 > v0 else 1.0
            w = min(max(w, 0.0), 1.0)
            return BreakInfo(name, limit, unit, h0 + w * (h1 - h0),
                             d0 + w * (d1 - d0))
        return None


def load_factors(increment, max_increments):
    """Load factors k * increment up to 1.0, at most max_increments of
    them (the last one is clipped to 1.0)."""
    if not 0.0 < increment <= 1.0:
        raise ValueError("LOAD_INCREMENT must be in (0, 1], got %g."
                         % increment)
    n = int(math.ceil(1.0 / increment - 1e-9))
    return [round(min((k + 1) * increment, 1.0), 12)
            for k in range(min(n, max_increments))]


class IncrementalSolver(object):
    """Runs the incremental, iterative analysis for one model."""

    def __init__(self, model, instance, nodes, spring_model, depth_to_z,
                 job_name, work_dir, final_odb, num_cpus=1, settings=None,
                 report=None, mudline_depth=0.0, criteria=None):
        """nodes: [(label, depth below pile head)]; report: function(message)
        for the progress log; criteria: BreakCriteria (None: no break)."""
        self.report = report or (lambda message: None)
        self.model = model
        self.instance = instance
        self.nodes = sorted(nodes, key=lambda n: n[1])
        self.spring_model = spring_model
        self.depth_to_z = depth_to_z
        self.job_name = job_name
        self.work_dir = work_dir            # temporary job files
        self.final_odb = final_odb          # odb of the last converged level
        self.num_cpus = num_cpus
        self.settings = settings or IterationSettings()
        self.criteria = criteria or BreakCriteria()

        self.head_label = self.nodes[0][0]
        self.toe_label = self.nodes[-1][0]
        # node at the mudline (always a node, see geometry.py)
        self.mudline_label = min(
            self.nodes, key=lambda n: abs(n[1] - mudline_depth))[0]
        self.state = None                   # relaxed node results
        self.spring_set = None
        self.results = []
        self.history = []                   # (load factor, node results)
        self.break_info = None              # BreakInfo if a criterion hit
        self._last_converged = None         # (spring set, load magnitude)
        self._printed_warnings = set()
        # Aitken relaxation factor, carried over from level to level
        self._omega = self.settings.relaxation
        self._last_residual = None

    # ------------------------------------------------------------------ #
    # Force control
    # ------------------------------------------------------------------ #
    def run(self, load, angle, lever_arm, factors):
        """Apply the load levels factors * load (force control). Returns the
        list of LoadLevelResult; stops at the first load level whose job
        fails or that reaches a break criterion (break_info)."""
        if not self.spring_model.iterative and len(factors) > 1:
            return self._run_single_job(load, angle, lever_arm, factors)
        failed = False
        for increment, factor in enumerate(factors, start=1):
            magnitude = factor * load
            moment = add_horizontal_load(self.model, self.instance,
                                         magnitude, angle, lever_arm)
            print("Increment %d/%d: H = %g N (%.1f %%), M = %g Nm"
                  % (increment, len(factors), magnitude, 100.0 * factor,
                     moment))
            starts = [self._predict(factor)]
            if starts[0] is not None and self.spring_model.iterative:
                starts.append(None)         # retry from the cold start
            result = None
            level_start = time.time()
            for start in starts:
                self.state = start
                if start is None:           # cold start: initial factor
                    self._omega = self.settings.relaxation
                try:
                    result = self._solve_level(increment, factor, magnitude,
                                               moment)
                    break
                except RuntimeError as err:
                    print("  job failed: %s" % err)
                    if start is not None:
                        print("  retrying load level %d from the cold start."
                              % increment)
            if result is None:
                print("  Load level %d failed, stopping incremental loading."
                      % increment)
                self.report("load level %d/%d FAILED (job error), loading "
                            "stopped" % (increment, len(factors)))
                failed = True
                break
            self.results.append(result)
            self._last_converged = (self.spring_set, magnitude)
            self._report_level(result, len(factors),
                               time.time() - level_start)
            if self._check_break(len(factors)):
                break

        if failed and self._last_converged is not None:
            # Reset the model to the last converged level (matches final_odb)
            self.spring_set, magnitude = self._last_converged
            self._apply_springs()
            add_horizontal_load(self.model, self.instance, magnitude, angle,
                                lever_arm)
        return self.results

    def _run_single_job(self, load, angle, lever_arm, factors):
        """Non-iterative spring models: the springs do not depend on the load
        level, so all levels are calculated in one job. The load ramps up to
        the last level over the step, and one output frame is written per
        level (step time = factor / last factor). If a break criterion is
        reached before the last level, the job is repeated up to that level
        so that the kept odb and cae belong to it."""
        full = factors[-1]
        add_horizontal_load(self.model, self.instance, full * load, angle,
                            lever_arm)
        set_output_time_points(self.model, [f / full for f in factors])
        print("All %d load levels in one job (springs independent of the "
              "load level): H up to %g N" % (len(factors), full * load))
        self.spring_set = self.spring_model.compute(self.nodes, None)
        self._apply_springs()

        odb_path = os.path.join(self.work_dir, self.job_name + ".odb")
        self.report("one job for all %d load levels (springs independent of "
                    "the load level) started" % len(factors))
        job_start = time.time()
        try:
            run_job(self.model, self.job_name, self.work_dir, self.num_cpus)
        except RuntimeError as err:
            # Keep the load levels reached before the job stopped
            print("  job failed: %s" % err)
            self.report("job FAILED, keeping the load levels reached")
            if not os.path.exists(odb_path):
                return self.results
        job_time = time.time() - job_start

        node_frames = dict(read_node_results_all_frames(odb_path,
                                                        self.instance.name))
        for step_time, head_u, head_ur in read_head_results_all_frames(
                odb_path, self.instance.name):
            # frame -> nearest load level
            increment, factor = min(enumerate(factors, start=1),
                                    key=lambda f: abs(f[1] / full - step_time))
            if abs(factor / full - step_time) > 1e-6:
                continue                    # not a load level (cutback)
            magnitude = factor * load
            nodes = node_frames[step_time]
            self.results.append(LoadLevelResult(
                increment, factor, magnitude, abs(magnitude * lever_arm),
                head_u, head_ur, 1, True, nodes[self.mudline_label],
                nodes[self.toe_label]))
            print("Increment %d/%d: H = %g N, head U1 = %.5f m"
                  % (increment, len(factors), magnitude, head_u[0]))
            self._report_level(self.results[-1], len(factors), None)
            if self._check_break(len(factors)):
                break
        self.report("job finished after %s" % duration(job_time))
        if len(self.results) < len(factors) and self.break_info is None:
            print("  Load levels reached: %d of %d."
                  % (len(self.results), len(factors)))
        if not self.results:
            return self.results

        last = self.results[-1]
        if self.break_info is not None and last.increment < len(factors):
            # Repeat the job up to the last kept level (odb/cae match it)
            self.report("repeating the job up to load level %d (break "
                        "criterion)" % last.increment)
            add_horizontal_load(self.model, self.instance, last.load, angle,
                                lever_arm)
            set_output_time_points(self.model, [1.0])
            run_job(self.model, self.job_name, self.work_dir, self.num_cpus)
        shutil.copyfile(odb_path, self.final_odb)
        return self.results

    def _predict(self, factor):
        """Starting state for a load level, extrapolated from the converged
        levels with v ~ H^a (a = 1 with only one level, clamped to 1-3)."""
        if not self.history:
            return None
        f1, last = self.history[-1]
        exponent = 1.0
        if len(self.history) > 1:
            f0, before = self.history[-2]
            v0 = abs(before[self.head_label][0])
            v1 = abs(last[self.head_label][0])
            if v0 > 0.0 and v1 > 0.0 and f1 > f0:
                exponent = min(max(math.log(v1 / v0) / math.log(f1 / f0),
                                   1.0), 3.0)
        scale = (factor / f1) ** exponent
        return dict((label, (v * scale, psi * scale))
                    for label, (v, psi) in last.items())

    def _solve_level(self, increment, factor, magnitude, moment):
        settings = self.settings
        converged = False
        self._last_residual = None          # omega from the level before
        for iteration in range(1, settings.max_iterations + 1):
            self.spring_set = self.spring_model.compute(self.nodes,
                                                        self.state)
            self._apply_springs()

            odb_path = run_job(self.model, self.job_name, self.work_dir,
                               self.num_cpus)
            node_results = read_node_results(odb_path, self.instance.name)
            residual = self._residual(node_results)
            print("  iteration %d: head v = %.5f m, psi = %.6f rad, "
                  "residual %s"
                  % ((iteration,) + node_results[self.head_label]
                     + ("-" if residual is None else "%.2e" % residual,)))

            if not self.spring_model.iterative:
                converged = True
                break
            if residual is not None and residual <= settings.tolerance:
                converged = True
                break
            self.state = self._relax(node_results)
            print("    relaxation %.3f" % self._omega)
        if self.spring_model.iterative:
            print("  %s after %d iterations."
                  % ("converged" if converged else "WARNING: NOT converged",
                     iteration))

        self.history.append((factor, node_results))
        shutil.copyfile(odb_path, self.final_odb)
        return self._level_result(increment, factor, magnitude, moment,
                                  odb_path, node_results, iteration,
                                  converged)

    # ------------------------------------------------------------------ #
    # Displacement control
    # ------------------------------------------------------------------ #
    def run_displacement(self, target, angle, lever_arm, factors,
                         load_guess):
        """Load levels with the mudline displacements factors * target [m]
        (displacement control); H of every level is found iteratively.
        load_guess [N]: estimate of H at the full target, only used as the
        starting value of the first level. Returns the list of
        LoadLevelResult; stops at the first level that fails or reaches a
        break criterion (break_info)."""
        for increment, factor in enumerate(factors, start=1):
            level_target = factor * target
            print("Increment %d/%d: mudline displacement %g m (%.1f %%)"
                  % (increment, len(factors), level_target, 100.0 * factor))
            level_start = time.time()
            try:
                result = self._solve_displacement_level(
                    increment, factor, level_target, angle, lever_arm,
                    load_guess)
            except RuntimeError as err:
                print("  Load level %d failed: %s" % (increment, err))
                self.report("load level %d/%d FAILED (%s), loading stopped"
                            % (increment, len(factors), err))
                if self._last_converged is not None:
                    self.spring_set, magnitude = self._last_converged
                    self._apply_springs()
                    add_horizontal_load(self.model, self.instance, magnitude,
                                        angle, lever_arm)
                break
            self.results.append(result)
            self._last_converged = (self.spring_set, result.load)
            self._report_level(result, len(factors),
                               time.time() - level_start)
            if self._check_break(len(factors)):
                break
        return self.results

    def _first_load(self, target, factor, load_guess):
        """Starting H for a level: extrapolated from the last two levels
        (secant through the origin with only one level)."""
        points = [(r.load, abs(r.mudline[0])) for r in self.results]
        if len(points) >= 2:
            (h0, v0), (h1, v1) = points[-2:]
            if v1 > v0:
                return max(h1 + (target - v1) * (h1 - h0) / (v1 - v0), h1)
        if points and points[-1][1] > 0.0:
            h1, v1 = points[-1]
            return h1 * target / v1
        return load_guess * factor

    def _reference_stiffness(self):
        """dH/dv between the last two converged levels (secant through the
        origin with one level; None before the first level)."""
        points = [(r.load, abs(r.mudline[0])) for r in self.results]
        if len(points) >= 2:
            (h0, v0), (h1, v1) = points[-2:]
            if v1 > v0:
                return (h1 - h0) / (v1 - v0)
        if points and points[-1][1] > 0.0:
            return points[-1][0] / points[-1][1]
        return None

    def _next_load(self, points, target, lower, upper, k_ref):
        """Next H of a level: step H + K (target - v).

        Non-iterative models (fixed springs): K = secant of the last two
        runs (through the origin after the first run).
        Iterative models: the springs change between the runs, so a single
        secant can be far off. The step converges as long as 0 < K < 2 x
        the current tangent; a too large K oscillates or diverges, a too
        small K only converges more slowly. K is therefore the smallest
        positive estimate of: secant of the last two runs, stiffness
        between the last two levels (k_ref), secant through the origin of
        the first run of the level.
        The result is kept inside the bracket (lower, upper) of H;
        bisection if the step leaves it."""
        h1, v1 = points[-1]
        secant = None
        if len(points) >= 2:
            h0, v0 = points[-2]
            if v1 != v0 and (h1 - h0) / (v1 - v0) > 0.0:
                secant = (h1 - h0) / (v1 - v0)
        if not self.spring_model.iterative:
            slope = secant or (h1 / v1 if v1 > 0.0 else None)
        else:
            estimates = [k for k in (secant, k_ref) if k is not None]
            if points[0][1] > 0.0:
                estimates.append(points[0][0] / points[0][1])
            slope = min(estimates) if estimates else None
        load = h1 + (target - v1) * slope if slope else 2.0 * h1
        if load <= 0.0:
            load = 0.5 * h1
        if self.spring_model.iterative:
            # Only the last run is a valid bound (current springs); older
            # bounds of successful runs are not kept (see above)
            if v1 < target:
                lower = max(lower, h1)
            elif upper is None or h1 < upper:
                upper = h1
        if upper is not None:
            if not lower < load < upper:
                load = 0.5 * (lower + upper)
        elif load <= lower:
            load = 2.0 * lower if lower > 0.0 else h1
        return load

    def _solve_displacement_level(self, increment, factor, target, angle,
                                  lever_arm, load_guess):
        settings = self.settings
        iterative = self.spring_model.iterative
        load = self._first_load(target, factor, load_guess)
        k_ref = self._reference_stiffness()
        lower, upper = 0.0, None            # H bracket: v too small / large
        points = []                         # (H, v_mudline) of this level
        accepted = None                     # last successful run
        level_odb = os.path.join(self.work_dir, "level.odb")
        converged = False
        self._last_residual = None          # omega from the level before

        # Starting state: last deflection line scaled to the new target
        if iterative and self.history:
            _, last = self.history[-1]
            v_last = last[self.mudline_label][0]
            scale = target / v_last if v_last > 0.0 else 1.0
            self.state = dict((label, (v * scale, psi * scale))
                              for label, (v, psi) in last.items())

        for iteration in range(1, settings.max_iterations + 1):
            if iterative or self.spring_set is None:
                self.spring_set = self.spring_model.compute(self.nodes,
                                                            self.state)
                self._apply_springs()
            moment = add_horizontal_load(self.model, self.instance, load,
                                         angle, lever_arm)
            try:
                odb_path = run_job(self.model, self.job_name, self.work_dir,
                                   self.num_cpus)
            except RuntimeError as err:
                # e.g. H beyond the capacity: too large
                print("  iteration %d: H = %g N, job failed (%s)"
                      % (iteration, load, err))
                upper = load if upper is None else min(upper, load)
                load = 0.5 * (lower + upper)
                continue

            node_results = read_node_results(odb_path, self.instance.name)
            shutil.copyfile(odb_path, level_odb)
            v = node_results[self.mudline_label][0]
            points.append((load, v))
            # A bound from a run is only valid for its springs. Iterative
            # models change the springs every run, so only the bound from a
            # failed job (H beyond the capacity) is kept for them.
            if not iterative:
                if v < target:
                    lower = max(lower, load)
                else:
                    upper = load if upper is None else min(upper, load)
            load_error = abs(v - target) / target
            residual = self._residual(node_results) if iterative else None
            accepted = (load, moment, node_results, self.spring_set)
            print("  iteration %d: H = %g N, mudline v = %.5f m (target "
                  "%.5f), load error %.2e%s"
                  % (iteration, load, v, target, load_error,
                     "" if not iterative else ", spring residual %s"
                     % ("-" if residual is None else "%.2e" % residual)))

            springs_ok = not iterative or (
                residual is not None and residual <= settings.tolerance)
            if load_error <= settings.tolerance and springs_ok:
                converged = True
                break
            load = self._next_load(points, target, lower, upper, k_ref)
            if iterative:
                scale = target / v if v > 0.0 else 1.0
                self.state = self._relax(dict(
                    (label, (a * scale, b * scale))
                    for label, (a, b) in node_results.items()))

        if accepted is None:
            raise RuntimeError("no job of the level ran successfully")
        load, moment, node_results, spring_set = accepted
        print("  %s after %d iterations, H = %g N."
              % ("converged" if converged else "WARNING: NOT converged",
                 iteration, load))
        # Model, odb and history of the last successful run
        if spring_set is not self.spring_set:
            self.spring_set = spring_set
            self._apply_springs()
        add_horizontal_load(self.model, self.instance, load, angle,
                            lever_arm)
        self.history.append((factor, node_results))
        shutil.copyfile(level_odb, self.final_odb)
        return self._level_result(increment, factor, load, moment,
                                  level_odb, node_results, iteration,
                                  converged)

    # ------------------------------------------------------------------ #
    # Common
    # ------------------------------------------------------------------ #
    def _level_result(self, increment, factor, load, moment, odb_path,
                      node_results, iterations, converged):
        head_u, head_ur = read_head_results(odb_path, self.instance.name)
        return LoadLevelResult(increment, factor, load, moment, head_u,
                               head_ur, iterations, converged,
                               node_results[self.mudline_label],
                               node_results[self.toe_label])

    def _check_break(self, n_levels):
        """Check the break criteria after the last level; report and store
        the result. True if the loading has to stop."""
        if not self.criteria:
            return False
        info = self.criteria.check(self.results)
        if info is None:
            return False
        self.break_info = info
        print("  Break criterion: %s" % info)
        last = self.results[-1].increment
        self.report("BREAK after load level %d/%d: %s"
                    % (last, n_levels, info))
        return True

    def _report_level(self, result, n_levels, seconds):
        """Progress message of a finished load level (seconds: its run time,
        None if all levels ran in one job)."""
        v, psi = result.mudline
        message = ("load level %d/%d done: H = %.4g MN, mudline displacement "
                   "%.2f mm, rotation %.3f deg"
                   % (result.increment, n_levels, result.load / 1e6,
                      1000.0 * abs(v), math.degrees(abs(psi))))
        if self.spring_model.iterative or result.iterations > 1:
            message += ", %d iterations%s" % (
                result.iterations,
                "" if result.converged else ", NOT converged")
        if seconds is not None:
            message += " (%s)" % duration(seconds)
        self.report(message)

    def _apply_springs(self):
        create_soil_connectors(
            self.model, self.instance, self.spring_set, self.depth_to_z,
            rotational=self.spring_model.has_rotational_springs)
        for message in self.spring_model.warnings:
            if message not in self._printed_warnings:
                print("  WARNING: %s" % message)
                self._printed_warnings.add(message)

    def _relax(self, new_results):
        """Relaxed node displacements for the next spring update, with the
        Aitken relaxation factor (Irons & Tuck; Kuettler & Wall 2008):

            r_k     = F(x_k) - x_k               (normalised residual)
            omega_k = -omega_k-1 * r_k-1 . (r_k - r_k-1) / |r_k - r_k-1|^2
            x_k+1   = x_k + omega_k r_k

        The first update of a level uses the factor of the level before
        (RELAXATION at the start and after a cold start); omega is limited
        to 0.1 ... 1.0. After the cold start (no state) the result itself is
        the next state."""
        if self.state is None:
            return new_results
        labels = sorted(new_results)
        scales = [max(abs(new_results[l][c]) for l in labels) or 1.0
                  for c in (0, 1)]
        residual = [(new_results[l][c] - self.state[l][c]) / scales[c]
                    for l in labels for c in (0, 1)]
        if self._last_residual is not None:
            delta = [r - r0 for r, r0 in zip(residual, self._last_residual)]
            denominator = sum(d * d for d in delta)
            if denominator > 0.0:
                omega = -self._omega * sum(
                    r0 * d for r0, d in zip(self._last_residual, delta)) \
                    / denominator
                self._omega = min(max(omega, 0.1), 1.0)
        self._last_residual = residual
        w = self._omega
        return dict((label, tuple(o + w * (n - o)
                                  for o, n in zip(self.state[label], new)))
                    for label, new in new_results.items())

    def _residual(self, new_results):
        """Fixed-point residual of the spring iteration: largest difference
        between the computed displacements/rotations and the state the
        springs were built from, relative to the largest value along the
        pile (max over v and psi). None for the cold start (no state)."""
        if self.state is None:
            return None
        residual = 0.0
        for component in (0, 1):
            scale = max(abs(values[component])
                        for values in new_results.values())
            if scale <= 0.0:
                continue
            diff = max(abs(new_results[label][component]
                           - self.state[label][component])
                       for label in new_results)
            residual = max(residual, diff / scale)
        return residual
