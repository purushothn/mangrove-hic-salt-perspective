#!/usr/bin/env python3
"""Density-matched permutation tests for gene sets (A-compartment share, distance to domain boundary).

Usage: python hic_followup.py ANALYSIS_DIR GENESET_DIR [N_PERM] [COUNT_TABLE]
If COUNT_TABLE (genes x libraries, tab separated, first column gene ID) is given, background genes are also matched on
mean expression (CPM), so a difference cannot be explained by detectability or expression level.
Reads ANALYSIS_DIR/genes_with_compartment.csv, phasing.tsv, ins.tsv (made by hic_analysis.py).
For each gene set, background genes are sampled with the same distribution of local gene density
(genes per 100 kb bin, binned into quantile classes), so a difference cannot be explained by gene position
along the chromosome. Boundary distance is computed for the 150 kb and 250 kb insulation windows.
"""
import sys, pathlib
import numpy as np, pandas as pd

adir, gdir = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
COUNTS = sys.argv[4] if len(sys.argv) > 4 else None
rng = np.random.default_rng(1)

g = pd.read_csv(adir / "genes_with_compartment.csv")
ph = pd.read_csv(adir / "phasing.tsv", sep="\t")
ph["bin"] = ph.start // 100_000
g = g.merge(ph[["chrom", "bin", "gene_density"]], on=["chrom", "bin"], how="left")
g = g.dropna(subset=["comp", "gene_density"])
g = g[g.comp.isin(["A", "B"])].copy()
g["isA"] = (g.comp == "A").astype(float)
ins = pd.read_csv(adir / "ins.tsv", sep="\t")
for w in (150000, 250000):
    col = f"is_boundary_{w}"
    if col not in ins.columns:
        continue
    b = ins[ins[col].fillna(False).astype(bool)].copy()
    b["pos"] = (b.start + b.end) // 2
    pos = {c: grp.pos.values for c, grp in b.groupby("chrom")}
    g[f"d{w}"] = [np.abs(pos[c] - m).min() / 1000 if c in pos else np.nan for c, m in zip(g.chrom, g.mid)]
dcols = [c for c in g.columns if c.startswith("d") and c[1:].isdigit()]

def assign_classes(g):
    g = g.copy()
    if COUNTS:
        ct = pd.read_csv(COUNTS, sep="\t", index_col=0)
        ct = ct[[c for c in ct.columns if c not in ("Root_28d", "Leaf_28d")]]
        mcpm = (ct / ct.sum() * 1e6).mean(axis=1).rename("mcpm")
        g = g.merge(mcpm, left_on="gene", right_index=True, how="inner").reset_index(drop=True)
        g["dclass"] = (pd.qcut(g.gene_density.rank(method="first"), 5, labels=False) * 10
                       + pd.qcut(g.mcpm.rank(method="first"), 5, labels=False))
        print("Matching on gene density (5 classes) x mean expression (5 classes)")
    else:
        g["dclass"] = pd.qcut(g.gene_density.rank(method="first"), 10, labels=False)
    return g

bgf = gdir / "background.txt"
bg = {x.strip() for x in open(bgf) if x.strip()} if bgf.exists() else None
if bg is not None:
    g = g[g.gene.isin(bg)]
g = assign_classes(g).reset_index(drop=True)

rows = []
for f in sorted(gdir.glob("*.txt")):
    if f.stem in {"background", "highlight"}:
        continue
    ids = {x.strip() for x in open(f) if x.strip()}
    s = g[g.gene.isin(ids)]
    pool = g[~g.gene.isin(ids)]
    if len(s) < 20:
        print(f"skip {f.stem}: only {len(s)} genes (too few for a permutation test)")
        continue
    need = s.dclass.value_counts()
    by = {k: pool.index[pool.dclass == k].values for k in need.index}
    empty = [k for k in need.index if len(by[k]) == 0]
    if empty:
        n_drop = int(need[empty].sum())
        print(f"{f.stem}: {n_drop} genes dropped (no comparison genes in their strata)")
        need = need.drop(empty)
        s = s[s.dclass.isin(need.index)]
    obs = {"pct_A": 100 * s.isA.mean(), **{c: s[c].median() for c in dcols}}
    sims = {k: [] for k in obs}
    for _ in range(NP):
        idx = np.concatenate([rng.choice(by[k], n, replace=len(by[k]) < n) for k, n in need.items()])
        t = pool.loc[idx]
        sims["pct_A"].append(100 * t.isA.mean())
        for c in dcols:
            sims[c].append(t[c].median())
    r = {"set": f.stem, "n_genes": len(s)}
    for k, v in obs.items():
        v_s = np.array(sims[k])
        r[f"obs_{k}"] = round(v, 2)
        r[f"matched_mean_{k}"] = round(v_s.mean(), 2)
        r[f"p_{k}"] = (1 + np.sum(np.abs(v_s - v_s.mean()) >= abs(v - v_s.mean()))) / (NP + 1)
    rows.append(r)

res = pd.DataFrame(rows)
res.to_csv(adir / ("geneset_tests_density_expr_matched.csv" if COUNTS else "geneset_tests_density_matched.csv"), index=False)
print(res.to_string(index=False))
