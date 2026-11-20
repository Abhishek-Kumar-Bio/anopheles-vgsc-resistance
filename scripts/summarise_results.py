#!/usr/bin/env python3
"""Summary statistics and figures from the notebook output in results/tables/.

    python scripts/summarise_results.py

Writes figures/vgsc_trends.png, figures/vgsc_linkage.png and results/summary.md.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / "results" / "tables"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3de", "#fcfcfb"
COL = {"V402L": "#2a78d6", "L995F": "#eb6834"}
SITES = {"BF-09": "Burkina Faso, region BF-09", "BF-02": "Burkina Faso, region BF-02", "ML-2": "Mali, region ML-2"}
plt.rcParams.update({"font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": GRID,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK2})


def wilson(k, n, z=1.96):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def main():
    long = pd.read_csv(T / "focal_alleles_by_cohort.csv")
    sets = pd.read_csv(T / "sample_sets.csv", dtype={"release": str})
    rel = sorted(sets.release.unique(), key=lambda v: tuple(int(x) for x in v.split(".")))
    w = long.pivot_table(index=["cohort", "country", "admin1", "taxon", "year", "n_mosquitoes"],
                         columns="allele", values="frequency").reset_index()

    # ---- figure 1: change over time at the most often sampled An. coluzzii sites ----
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, (site, title) in zip(axes, SITES.items()):
        for allele in ("L995F", "V402L"):
            d = long[(long.admin1 == site) & (long.taxon == "coluzzii") & (long.allele == allele) & (long.year >= 2010)].sort_values("year")
            lo, hi = wilson(d.alt_count.to_numpy(), d.alleles_called.to_numpy())
            ax.fill_between(d.year, lo, hi, color=COL[allele], alpha=0.15, lw=0)
            ax.plot(d.year, d.frequency, "-o", color=COL[allele], lw=2, ms=4.5)
            last = d.iloc[-1]
            ax.annotate(f"{last.frequency:.0%}", (last.year, last.frequency), xytext=(5, 6 if allele == "L995F" else -6),
                        textcoords="offset points", va="center", fontsize=8.5, color=INK)
        ax.set_title(title, loc="left", fontsize=9.5, color=INK)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, .25, .5, .75, 1], ["0", "25%", "50%", "75%", "100%"])
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_xlim(2010, 2024)
        ax.set_xticks([2011, 2014, 2017, 2020, 2023])
    axes[0].set_ylabel("Allele frequency (95% CI)")
    fig.legend(handles=[Line2D([], [], color=COL[a], marker="o", lw=2, label=l) for a, l in
                        (("L995F", "L995F (classic kdr)"), ("V402L", "V402L (with I1527T)"))],
               loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("In An. coluzzii, V402L has risen while L995F has fallen", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.text(0.01, 0.885, "Vgsc allele frequencies by collection year (2010 onwards) at the three most often sampled sites. "
             "Data: MalariaGEN Ag3 (unrestricted-use sample sets).", fontsize=8.5, color=INK2)
    fig.subplots_adjust(top=0.78, bottom=0.2, left=0.07, right=0.97, wspace=0.12)
    (ROOT / "figures").mkdir(exist_ok=True)
    fig.savefig(ROOT / "figures" / "vgsc_trends.png", dpi=200)
    plt.close(fig)

    # ---- figure 2: linkage and complementarity ----
    col = w[w.taxon == "coluzzii"]
    other = w[w.taxon != "coluzzii"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 4.2))
    for ax in (a1, a2):
        ax.plot([0, 1], [0, 1] if ax is a1 else [1, 0], color=INK2, lw=0.8, ls=(0, (3, 2)))
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    a1.scatter(other.V402L, other.I1527T, s=14, color="#b4b2a9", lw=0, label="other species")
    a1.scatter(col.V402L, col.I1527T, s=16, color=COL["V402L"], lw=0, alpha=0.8, label="An. coluzzii")
    r = np.corrcoef(w.V402L, w.I1527T)[0, 1]
    a1.set_xlabel("V402L frequency")
    a1.set_ylabel("I1527T frequency")
    a1.set_title(f"V402L and I1527T travel together (r = {r:.4f})", loc="left", fontsize=9.5, color=INK)
    a1.legend(frameon=False, loc="upper left", fontsize=8)
    a2.scatter(other.V402L, other.L995F, s=14, color="#b4b2a9", lw=0)
    a2.scatter(col.V402L, col.L995F, s=16, color=COL["V402L"], lw=0, alpha=0.8)
    a2.set_xlabel("V402L frequency")
    a2.set_ylabel("L995F frequency")
    a2.set_title("V402L and L995F together never exceed 100%", loc="left", fontsize=9.5, color=INK)
    a2.text(0.62, 0.43, "V402L + L995F = 100%", rotation=-42, fontsize=8, color=INK2)
    fig.text(0.01, 0.01, f"Each dot is one cohort (same region, year and species, at least 10 mosquitoes); n = {len(w)} cohorts.",
             fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(ROOT / "figures" / "vgsc_linkage.png", dpi=200)
    plt.close(fig)

    # ---- summary.md ----
    pres = long[long.frequency > 0]
    L = ["# Results summary", "", "Generated by `scripts/summarise_results.py` from `results/tables/`.", "",
         "## Data", "",
         f"- {sets.sample_count.sum():,} mosquitoes in {len(sets)} unrestricted-use sample sets (Ag3 releases "
         f"{rel[0]} to {rel[-1]}); all sets have `unrestricted_use = True`.",
         f"- {len(w)} cohorts (region x year x species, at least 10 mosquitoes) entered the frequency analysis.", "",
         "## The five focal alleles", "", "| Allele | Cohorts where found | Highest frequency | Highest in |", "|---|---|---|---|"]
    for a in ["L995F", "L995S", "N1570Y", "V402L", "I1527T"]:
        p = pres[pres.allele == a]
        top = p.loc[p.frequency.idxmax()]
        L.append(f"| {a} | {p.cohort.nunique()} of {long.cohort.nunique()} | {top.frequency:.0%} | {top.cohort} |")
    v = w[(w.V402L > 0) | (w.I1527T > 0)]
    L += ["", "## V402L and I1527T", "",
          f"- Correlation of V402L and I1527T frequency across cohorts: r = {np.corrcoef(w.V402L, w.I1527T)[0, 1]:.4f}; "
          f"largest difference within a cohort {(v.V402L - v.I1527T).abs().max():.3f}.",
          f"- Cohorts with V402L but no I1527T: {((w.V402L > 0) & (w.I1527T == 0)).sum()} "
          f"(V402L at most {w.loc[(w.V402L > 0) & (w.I1527T == 0), 'V402L'].max():.1%} in these); with I1527T but no V402L: "
          f"{((w.I1527T > 0) & (w.V402L == 0)).sum()}.",
          f"- Cohorts with V402L at 5% or more: {(w.V402L >= .05).sum()}, of which An. coluzzii: "
          f"{((w.V402L >= .05) & (w.taxon == 'coluzzii')).sum()}.",
          f"- Countries with V402L at 5% or more: {', '.join(sorted(w[w.V402L >= .05].country.unique()))}.",
          f"- Highest V402L + L995F in any cohort: {(w.V402L + w.L995F).max():.3f}.",
          f"- Cohorts where N1570Y exceeds L995F by more than 2 points: {(w.N1570Y > w.L995F + .02).sum()}.", "",
          "## Change over time, An. coluzzii (2010 onwards)", "", "| Site | Allele | First year | Frequency | Last year | Frequency |", "|---|---|---|---|---|---|"]
    for site in SITES:
        for a in ("V402L", "L995F"):
            d = long[(long.admin1 == site) & (long.taxon == "coluzzii") & (long.allele == a) & (long.year >= 2010)].sort_values("year")
            L.append(f"| {site} | {a} | {int(d.year.iloc[0])} | {d.frequency.iloc[0]:.0%} | {int(d.year.iloc[-1])} | {d.frequency.iloc[-1]:.0%} |")
    (ROOT / "results" / "summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
