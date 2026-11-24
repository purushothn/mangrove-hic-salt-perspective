#!/usr/bin/env python3
"""Robustness tests for the A-compartment enrichment of gene sets (reviewer response).

Usage: python hic_robustness.py ANALYSIS_DIR GENESET_DIR COUNT_TABLE OUTDIR [N_PERM]

Reads ANALYSIS_DIR/genes_with_compartment.csv and phasing.tsv (from hic_analysis.py).
Uses only expressed background genes (GENESET_DIR/background.txt) with an A/B call.
Tests, for each gene set:
  1. stratified permutation: comparison genes drawn within chromosome x local gene density
     (tertiles) x mean expression (tertiles), so chromosome composition is held fixed;
  2. block (circular-shift) permutation: within each chromosome the gene-set labels are rotated
     along the gene order by a random offset, which keeps the clustering of set genes and the
     run structure of compartments; (a) raw A share, no matching; (b) covariate-adjusted
     statistic: each gene's A call minus the mean A share of its chromosome x density x
     expression stratum, so the test keeps spatial clustering and also adjusts for covariates;
  3. Cochran-Mantel-Haenszel test stratified by chromosome, plus the number of chromosomes on
     which the set has a higher A share than other expressed genes.
Leaf-only and root-only subsets (genes responsive in one organ but not the other) are added.
"""
import sys, pathlib
import numpy as np, pandas as pd
from statsmodels.stats.contingency_tables import StratifiedTable

adir, gdir, counts, odir = (pathlib.Path(a) for a in sys.argv[1:5])
NP = int(sys.argv[5]) if len(sys.argv) > 5 else 10000
odir.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(7)

g = pd.read_csv(adir / "genes_with_compartment.csv")
ph = pd.read_csv(adir / "phasing.tsv", sep="\t")
ph["bin"] = ph.start // 100_000
g = g.merge(ph[["chrom", "bin", "gene_density"]], on=["chrom", "bin"], how="left")
g = g[g.comp.isin(["A", "B"])].dropna(subset=["gene_density"]).copy()
bg = {x.strip() for x in open(gdir / "background.txt") if x.strip()}
g = g[g.gene.isin(bg)]
ct = pd.read_csv(counts, sep="\t", index_col=0)
ct = ct[[c for c in ct.columns if c not in ("Root_28d", "Leaf_28d")]]
mcpm = (ct / ct.sum() * 1e6).mean(axis=1).rename("mcpm")
g = g.merge(mcpm, left_on="gene", right_index=True, how="inner")
g = g.sort_values(["chrom", "start"]).reset_index(drop=True)
g["isA"] = (g.comp == "A").astype(float)
g["d3"] = pd.qcut(g.gene_density.rank(method="first"), 3, labels=False)
g["e3"] = pd.qcut(g.mcpm.rank(method="first"), 3, labels=False)
g["stratum"] = g.chrom + "_" + g.d3.astype(str) + g.e3.astype(str)

def read(n):
    return {x.strip() for x in open(gdir / f"{n}.txt") if x.strip()}

sets = {n: read(n) for n in ("salinome", "deg_leaf", "deg_root", "literature_salt")}
sets["leaf_only"] = sets["deg_leaf"] - sets["deg_root"]
sets["root_only"] = sets["deg_root"] - sets["deg_leaf"]
sets["leaf_and_root"] = sets["deg_leaf"] & sets["deg_root"]

def p2(obs, sims):
    sims = np.asarray(sims)
    return (1 + np.sum(np.abs(sims - sims.mean()) >= abs(obs - sims.mean()))) / (len(sims) + 1)

chrom_idx = [grp.index.values for _, grp in g.groupby("chrom")]
isA = g.isA.values
resid = (g.isA - g.groupby("stratum").isA.transform("mean")).values
rows = []
for name, ids in sets.items():
    lab = g.gene.isin(ids).values
    s, pool = g[lab], g[~lab]
    obs = 100 * s.isA.mean()
    # 1. chromosome x density x expression stratified permutation
    need = s.stratum.value_counts()
    by = {k: pool.index[pool.stratum == k].values for k in need.index}
    keep = [k for k in need.index if len(by[k]) > 0]
    n_drop = int(need.drop(keep).sum())
    need = need[keep]
    obs1 = 100 * s[s.stratum.isin(keep)].isA.mean()
    sim1 = [100 * isA[np.concatenate([rng.choice(by[k], n, replace=len(by[k]) < n)
                                      for k, n in need.items()])].mean() for _ in range(NP)]
    # 2. circular-shift block permutation within chromosomes
    sim2 = []
    for _ in range(NP):
        tot = cnt = 0
        for ix in chrom_idx:
            l = np.roll(lab[ix], rng.integers(1, len(ix)))
            tot += isA[ix][l].sum(); cnt += l.sum()
        sim2.append(100 * tot / cnt)
    # 2b. circular-shift permutation of a covariate-adjusted statistic: residual = isA minus the
    #     mean A share of all expressed genes in the same chromosome x density x expression stratum
    obs2b = 100 * resid[lab].mean()
    sim2b = []
    for _ in range(NP):
        tot = cnt = 0
        for ix in chrom_idx:
            l = np.roll(lab[ix], rng.integers(1, len(ix)))
            tot += resid[ix][l].sum(); cnt += l.sum()
        sim2b.append(100 * tot / cnt)
    # 3. CMH stratified by chromosome
    tabs, up, nch = [], 0, 0
    for c, grp in g.groupby("chrom"):
        l = grp.gene.isin(ids)
        t = np.array([[grp.isA[l].sum(), (1 - grp.isA[l]).sum()],
                      [grp.isA[~l].sum(), (1 - grp.isA[~l]).sum()]])
        if l.sum() >= 5:
            nch += 1
            up += grp.isA[l].mean() > grp.isA[~l].mean()
        tabs.append(t + 0.5)
    st = StratifiedTable(tabs)
    rows.append({
        "set": name, "n_genes": int(lab.sum()), "obs_pct_A": round(obs, 2),
        "strat_chrom_dens_expr_expected": round(np.mean(sim1), 2), "strat_p": p2(obs1, sim1),
        "strat_genes_dropped": n_drop,
        "circshift_expected": round(np.mean(sim2), 2), "circshift_p": p2(obs, sim2),
        "adj_excess_pctpts": round(obs2b, 2), "circshift_adj_p": p2(obs2b, sim2b),
        "cmh_OR": round(st.oddsratio_pooled, 3), "cmh_p": st.test_null_odds().pvalue,
        "chrom_higher": f"{up}/{nch}",
    })
res = pd.DataFrame(rows)
res.to_csv(odir / "geneset_tests_robustness.csv", index=False)
pd.set_option("display.width", 250)
print(res.to_string(index=False))
