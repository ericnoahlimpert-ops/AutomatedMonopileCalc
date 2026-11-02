"""Simple result plots (PNG) with matplotlib (included in Abaqus Python).

Per calculation (folder RESULTS_DIR/<pile>/<soil>/<model>):
    load_displacement.png   H over the mudline displacement and rotation of
                            all load levels; the break criterion (if
                            reached) is marked at its interpolated load
    pile_profile.png        last load level over the depth below mudline:
                            deflection v, shear force Q and bending moment M
                            (Abaqus values, see postprocessing.py)
For several calculations (RESULTS_DIR):
    load_displacement_all.png   load-displacement curves of all of them

The non-interactive Agg backend is used, so the plots also work without
GUI (abaqus cae noGUI=main.py). Units in the plots: MN, MNm, mm, deg.
"""

import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after the backend)

DPI = 150


def _curve(results):
    """(mudline displacement [mm], rotation [deg], H [MN]) of all load
    levels, starting at the origin."""
    v = [0.0] + [1000.0 * abs(r.mudline[0]) for r in results]
    psi = [0.0] + [math.degrees(abs(r.mudline[1])) for r in results]
    h = [0.0] + [r.load / 1e6 for r in results]
    return v, psi, h


def plot_load_displacement(results, path, title, break_info=None):
    """H over mudline displacement and rotation (results: list of
    solver.LoadLevelResult)."""
    v, psi, h = _curve(results)
    fig, (ax_v, ax_psi) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax_v.plot(v, h, "o-", color="C0")
    ax_psi.plot(psi, h, "o-", color="C1")
    if break_info is not None and break_info.load is not None:
        h_break = break_info.load / 1e6
        for ax in (ax_v, ax_psi):
            ax.axhline(h_break, color="C3", linestyle="--", linewidth=1,
                       label="%s: H = %.3g MN" % (break_info.name, h_break))
        ax_v.plot([1000.0 * break_info.displacement], [h_break], "x",
                  color="C3", markersize=10)
        ax_v.legend(loc="lower right", fontsize=8)
    ax_v.set_xlabel("Mudline displacement v [mm]")
    ax_psi.set_xlabel("Mudline rotation [deg]")
    for ax in (ax_v, ax_psi):
        ax.set_ylabel("Horizontal load H [MN]")
        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=0.0)
        ax.grid(True, alpha=0.3)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def plot_pile_profile(rows, path, title, load):
    """Deflection v, shear force Q and moment M over the depth below mudline
    (rows: postprocessing.ProfileRow of the last load level, load [N])."""
    z = [r.z_mudline for r in rows]
    fig, axes = plt.subplots(1, 3, figsize=(11, 6), sharey=True)
    for ax, values, label, color in (
            (axes[0], [1000.0 * r.v for r in rows],
             "Deflection v [mm]", "C0"),
            (axes[1], [r.q / 1e6 for r in rows], "Shear force Q [MN]", "C2"),
            (axes[2], [r.m / 1e6 for r in rows],
             "Bending moment M [MNm]", "C3")):
        ax.plot(values, z, "-", color=color)
        ax.axvline(0.0, color="black", linewidth=0.8)
        if z[0] < 0.0:                      # pile head above the mudline
            ax.axhline(0.0, color="saddlebrown", linewidth=1,
                       linestyle="--")
        ax.set_xlabel(label)
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("Depth below mudline z [m]")
    axes[0].invert_yaxis()
    fig.suptitle("%s - last load level, H = %.4g MN" % (title, load / 1e6))
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def plot_load_displacement_all(curves, path):
    """Load-displacement curves of several calculations in one plot
    (curves: [(label, results)])."""
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for label, results in curves:
        if results:
            v, _, h = _curve(results)
            ax.plot(v, h, "o-", markersize=3, label=label)
    ax.set_xlabel("Mudline displacement v [mm]")
    ax.set_ylabel("Horizontal load H [MN]")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
