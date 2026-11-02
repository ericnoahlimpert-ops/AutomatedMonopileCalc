"""Recreate Burd et al. (2020a), Fig. 13 (a) and (c), with this
implementation overlaid, plus the relative difference to the paper's 1D
model. Run with Abaqus Python (matplotlib) in this folder:

    abaqus python plot_fig13.py

Input (this folder):
    paper_Fig13_<pile>_3D.csv, paper_Fig13_<pile>_1D.csv
        curves of the paper, extracted from the vector graphics of the PDF
        (page 14) and scaled with the axis ticks
    pile_head.csv of this implementation: from
        results/verification/PISA_Verification/<pile>_curve (input folders
        Inputs/Verification/PISA_Verification/<pile>_curve) if that run
        exists, otherwise this_implementation_<pile>_pile_head.csv
Output: Fig13_comparison.png and Fig13_comparison.csv
"""
import csv
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.dirname(HERE)            # results/verification/PISA_Verification
OUT = os.path.join(HERE, "Fig13_comparison.png")
TABLE7 = {"C1": 25.5510, "C4": 174.3406}        # H_ult [MN], finest mesh
H_MAX = {"C1": 30.0, "C4": 200.0}


def read_curve(path):
    with open(path) as f:
        rows = list(csv.reader(f))[1:]
    return [float(r[0]) for r in rows], [float(r[1]) for r in rows]


def read_ours(pile):
    runs = glob.glob(os.path.join(RUNS, "%s_curve" % pile, "*", "*", "*",
                                  "pile_head.csv"))
    path = runs[0] if runs else os.path.join(
        HERE, "this_implementation_%s_pile_head.csv" % pile)
    v, h = [0.0], [0.0]
    with open(path) as f:
        for row in csv.DictReader(f):
            v.append(abs(float(row["v_mudline [m]"])))
            h.append(float(row["H [N]"]) / 1e6)
    return v, h


def interp(xs, ys, x):
    if x <= xs[0]:
        return ys[0] * x / xs[0] if xs[0] > 0 else ys[0]
    for i in range(1, len(xs)):
        if x <= xs[i]:
            t = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
            return ys[i - 1] + t * (ys[i] - ys[i - 1])
    return None


fig, axes = plt.subplots(2, 2, figsize=(11, 8.5),
                         gridspec_kw={"height_ratios": [2.2, 1]})
report = []
for col, (pile, letter) in enumerate((("C1", "a"), ("C4", "c"))):
    v3, h3 = read_curve(os.path.join(HERE, "paper_Fig13_%s_3D.csv" % pile))
    v1, h1 = read_curve(os.path.join(HERE, "paper_Fig13_%s_1D.csv" % pile))
    vo, ho = read_ours(pile)

    ax = axes[0][col]
    ax.plot(v3, h3, "k-", lw=1.6, label="3D finite-element (paper)")
    ax.plot(v1, h1, "k--", lw=1.6, label="1D GDSM (paper)")
    ax.plot(vo, ho, "-o", color="#c0392b", ms=4, lw=1.2,
            label="This implementation (Abaqus)")
    ax.plot([1.0], [TABLE7[pile]], "*", color="#2c7fb8", ms=12,
            label="Table 7: H_ult = %.2f MN" % TABLE7[pile])
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, H_MAX[pile])
    ax.set_xlabel("Ground-level displacement, v_G [m]")
    ax.set_ylabel("Lateral load, H [MN]")
    ax.set_title("(%s) Pile %s, D = 10 m, L = %d m, D_R = 75 %%"
                 % (letter, pile, 20 if pile == "C1" else 60), fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)

    # relative differences at our load levels
    ax = axes[1][col]
    vs, d1, d3, dpaper = [], [], [], []
    for v, h in zip(vo[1:], ho[1:]):
        p1, p3 = interp(v1, h1, v), interp(v3, h3, v)
        if p1 is None or p3 is None:
            continue
        vs.append(v)
        d1.append(100.0 * (h / p1 - 1.0))
        d3.append(100.0 * (h / p3 - 1.0))
        dpaper.append(100.0 * (p1 / p3 - 1.0))
        report.append((pile, v, h, p1, p3))
    ax.axhline(0.0, color="k", lw=0.8)
    ax.plot(vs, d1, "-o", color="#c0392b", ms=4, label="this impl. vs paper 1D")
    ax.plot(vs, d3, "-s", color="#7f7f7f", ms=3, label="this impl. vs 3D FE")
    ax.plot(vs, dpaper, ":", color="k", label="paper 1D vs 3D FE")
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Ground-level displacement, v_G [m]")
    ax.set_ylabel("Difference in H [%]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

fig.suptitle("PISA sand (GDSM), recreation of Burd et al. (2020a), Fig. 13 (a) "
             "and (c), with this implementation", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(OUT, dpi=150)

with open(os.path.splitext(OUT)[0] + ".csv", "w") as f:
    f.write("pile,v_G [m],H this impl. [MN],H paper 1D [MN],H paper 3D [MN],"
            "diff to 1D [%],diff to 3D [%]\n")
    for pile, v, h, p1, p3 in report:
        f.write("%s,%.4f,%.4f,%.4f,%.4f,%.2f,%.2f\n"
                % (pile, v, h, p1, p3, 100 * (h / p1 - 1), 100 * (h / p3 - 1)))
for pile, v, h, p1, p3 in report:
    print("%s v=%.3f  H=%8.3f  paper1D=%8.3f (%+.2f %%)  3D=%8.3f (%+.2f %%)"
          % (pile, v, h, p1, 100 * (h / p1 - 1), p3, 100 * (h / p3 - 1)))
print("written", OUT)
