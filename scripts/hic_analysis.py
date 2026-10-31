#!/usr/bin/env python3
"""A/B compartments and domain boundaries for A. marina Hi-C, and overlap with gene sets.

Usage:
  python hic_analysis.py AM.mcool VIEW.tsv ANNOTATION.gff3 GENESET_DIR OUTDIR [HIGHLIGHT_CHROM]

GENESET_DIR holds one text file per gene set (one gene ID per line, IDs as in the GFF3 'ID=' field),
for example salinome.txt, deg_leaf.txt, deg_root.txt, hkt.txt, aquaporins.txt.
Special files in that folder: background.txt (expressed genes, used as the comparison set for every test)
and highlight.txt (genes to mark on the chromosome figure); neither is tested as a gene set.
Not tested in the Claude workspace. Needs: cooler, cooltools, pandas, numpy, scipy, matplotlib.
"""
import re
import sys
import pathlib

import numpy as np
import pandas as pd
import cooler
from scipy.stats import fisher_exact, mannwhitneyu

mcool, view_path, gff, setdir, outdir = sys.argv[1:6]
highlight_chrom = sys.argv[6] if len(sys.argv) > 6 else None
out = pathlib.Path(outdir)
out.mkdir(parents=True, exist_ok=True)

COMP_RES = 100_000   # compartments
TAD_RES = 50_000     # insulation / boundaries
WINDOWS = [150_000, 250_000]

view = pd.read_csv(view_path, sep="\t", header=None, names=["chrom", "start", "end", "name"])
chroms = set(view.chrom)

# ---- genes from GFF3 ----
rows = []
with open(gff) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or f[0] not in chroms:
            continue
        m = re.search(r"ID=([^;]+)", f[8])
        if m:
            rows.append((m.group(1), f[0], int(f[3]), int(f[4])))
genes = pd.DataFrame(rows, columns=["gene", "chrom", "start", "end"])
genes["mid"] = (genes.start + genes.end) // 2
print(f"{len(genes)} genes on the {len(chroms)} chromosomes")

# ---- compartments (eigenvector oriented by gene density, so A = gene-rich) ----
clr = cooler.Cooler(f"{mcool}::resolutions/{COMP_RES}")
bins = clr.bins()[:][["chrom", "start", "end"]]
bins = bins[bins.chrom.isin(chroms)].copy()
bins["bin"] = bins.start // COMP_RES
genes["bin"] = genes.mid // COMP_RES
dens = genes.groupby(["chrom", "bin"]).size().rename("gene_density").reset_index()
track = bins.merge(dens, on=["chrom", "bin"], how="left").fillna({"gene_density": 0})
track[["chrom", "start", "end", "gene_density"]].to_csv(out / "phasing.tsv", sep="\t", index=False)

import cooltools
track_df = track[["chrom", "start", "end", "gene_density"]].reset_index(drop=True)
_, eig_df = cooltools.eigs_cis(clr, phasing_track=track_df, view_df=view, n_eigs=3)
eig_df.to_csv(out / "comp.cis.vecs.tsv", sep="\t", index=False)
vecs = eig_df.copy()
vecs["bin"] = vecs.start // COMP_RES
vecs["comp"] = np.where(vecs.E1 > 0, "A", np.where(vecs.E1 < 0, "B", "NA"))
genes = genes.merge(vecs[["chrom", "bin", "E1", "comp"]], on=["chrom", "bin"], how="left")
print("Share of genes in A:", round((genes.comp == "A").mean(), 3))

# ---- insulation boundaries ----
clr50 = cooler.Cooler(f"{mcool}::resolutions/{TAD_RES}")
ins = cooltools.insulation(clr50, WINDOWS, view_df=view, threshold="Li", verbose=False)
ins.to_csv(out / "ins.tsv", sep="\t", index=False)
bcol = f"is_boundary_{WINDOWS[0]}"
ins[bcol] = ins[bcol].fillna(False).astype(bool)
bnd = ins[ins[bcol].astype(bool)].copy()
bnd["pos"] = (bnd.start + bnd.end) // 2


def dist_to_boundary(row):
    p = bnd.loc[bnd.chrom == row.chrom, "pos"].values
    return np.abs(p - row.mid).min() if len(p) else np.nan


genes["dist_boundary"] = genes.apply(dist_to_boundary, axis=1)
print(f"{len(bnd)} boundaries at window {WINDOWS[0]}")

# ---- gene-set tests ----
results = []
# Optional background.txt = genes expressed in your RNA-seq (recommended for DEG sets).
bg_file = pathlib.Path(setdir) / "background.txt"
bg_ids = {x.strip() for x in open(bg_file) if x.strip()} if bg_file.exists() else None
print("Background:", f"{len(bg_ids)} genes from background.txt" if bg_ids else "all annotated genes")
for f in sorted(pathlib.Path(setdir).glob("*.txt")):
    if f.stem in {"background", "highlight"}:
        continue
    ids = {x.strip() for x in open(f) if x.strip()}
    called = genes[genes.comp.isin(["A", "B"])]   # genes in masked bins have no compartment call and are not counted as B
    g_in = called[called.gene.isin(ids)]
    g_out = called[~called.gene.isin(ids)]
    if bg_ids is not None:
        # compare like with like: both the set and the comparison genes must be in the expressed background
        g_in = g_in[g_in.gene.isin(bg_ids)]
        g_out = g_out[g_out.gene.isin(bg_ids)]
    if len(g_in) < 3:
        print(f"skip {f.name}: only {len(g_in)} genes found on chromosomes")
        continue
    a_in, a_out = (g_in.comp == "A").sum(), (g_out.comp == "A").sum()
    odds, p_comp = fisher_exact([[a_in, len(g_in) - a_in], [a_out, len(g_out) - a_out]])
    d_in, d_out = g_in.dist_boundary.dropna(), g_out.dist_boundary.dropna()
    p_dist = mannwhitneyu(d_in, d_out).pvalue
    results.append({
        "set": f.stem, "n_genes": len(g_in), "pct_in_A": round(100 * a_in / len(g_in), 1),
        "pct_in_A_background": round(100 * a_out / len(g_out), 1), "odds_ratio_A": round(odds, 2),
        "p_compartment": p_comp, "median_dist_boundary_kb": round(d_in.median() / 1000, 1),
        "background_median_kb": round(d_out.median() / 1000, 1), "p_boundary_distance": p_dist,
    })
res = pd.DataFrame(results)
if len(res):
    # Benjamini-Hochberg across sets for each test
    for col in ["p_compartment", "p_boundary_distance"]:
        order = res[col].rank(method="first").astype(int)
        res[col.replace("p_", "q_")] = np.minimum(1, res[col] * len(res) / order)
res.to_csv(out / "geneset_tests.csv", index=False)
genes.to_csv(out / "genes_with_compartment.csv", index=False)
print(res.to_string(index=False))

# ---- one figure: matrix, E1 and highlighted genes for one chromosome ----
if highlight_chrom:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    c50 = cooler.Cooler(f"{mcool}::resolutions/{TAD_RES}")
    mat = np.log2(c50.matrix(balance=True).fetch(highlight_chrom) + 1e-4)
    size = c50.chromsizes[highlight_chrom]
    e1 = vecs[vecs.chrom == highlight_chrom]
    hl = pathlib.Path(setdir) / "highlight.txt"
    hl_ids = {x.strip() for x in open(hl)} if hl.exists() else set()
    hg = genes[(genes.chrom == highlight_chrom) & genes.gene.isin(hl_ids)]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7.5), gridspec_kw={"height_ratios": [5, 1.2]}, sharex=False)
    ax1.imshow(mat, cmap="Reds", extent=[0, size / 1e6, size / 1e6, 0], vmin=np.nanpercentile(mat, 5), vmax=np.nanpercentile(mat, 99))
    ax1.set_ylabel(f"{highlight_chrom} (Mb)")
    ax2.fill_between(e1.start / 1e6, e1.E1.fillna(0), 0, where=e1.E1 > 0, color="#c0392b", label="A")
    ax2.fill_between(e1.start / 1e6, e1.E1.fillna(0), 0, where=e1.E1 < 0, color="#2c3e50", label="B")
    ax2.plot(hg.mid / 1e6, np.full(len(hg), e1.E1.min() * 0.9), marker="^", ls="none", color="black", ms=5, label="HKT")
    ax2.set_xlim(0, size / 1e6)
    ax2.set_xlabel("Position (Mb)")
    ax2.set_ylabel("E1")
    ax2.legend(frameon=False, loc="upper center", ncol=3, fontsize=8)
    fig.tight_layout()
    fig.canvas.draw()
    p1, p2 = ax1.get_position(), ax2.get_position()
    ax2.set_position([p1.x0, p2.y0, p1.width, p2.height])  # align the E1 track with the matrix
    fig.savefig(out / f"{highlight_chrom}_hic.png", dpi=300)
    print("figure:", out / f"{highlight_chrom}_hic.png")
