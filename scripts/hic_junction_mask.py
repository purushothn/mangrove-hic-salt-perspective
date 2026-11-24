#!/usr/bin/env python3
"""Assembly-junction checks for the A/B compartment analysis (reviewer response).

Usage: python hic_junction_mask.py ANALYSIS_DIR GAPS.txt.gz ASSEMBLY_REPORT.txt OUTDIR [FLANK_KB ...]

GAPS: NCBI *_genomic_gaps.txt.gz for GCA_019155195.1 (AmGA_v1.0); every gap is treated as a
possible contig/scaffold junction joined during Hi-C scaffolding (conservative).
1. Writes OUTDIR/mask_<F>kb/ copies of genes_with_compartment.csv and phasing.tsv with genes
   within F kb of any gap removed, ready for hic_robustness.py.
2. Reports whether A/B compartment transitions (100 kb bins) fall near gaps more often than
   expected if transitions were placed at random bins on the same chromosome.
"""
import sys, pathlib, shutil
import numpy as np, pandas as pd

adir, gapf, rep, odir = (pathlib.Path(a) for a in sys.argv[1:5])
flanks = [int(x) for x in sys.argv[5:]] or [100, 250]
rng = np.random.default_rng(3)

cols = ["name", "role", "mol", "loc", "acc", "rel", "refseq", "unit", "length", "ucsc"]
r = pd.read_csv(rep, sep="\t", comment="#", header=None, names=cols)
acc2chr = dict(zip(r.acc, r.name))
gaps = pd.read_csv(gapf, sep="\t", comment="#", header=None,
                   names=["acc", "start", "stop", "len", "type", "evidence"])
gaps["chrom"] = gaps.acc.map(acc2chr)
gaps = gaps[gaps.chrom.str.match(r"chr\d\d$", na=False)]
print(f"{len(gaps)} gaps on chr01-chr31 (median length {gaps.len.median():.0f} bp)")
gpos = {c: ((d.start + d.stop) // 2).values for c, d in gaps.groupby("chrom")}

g = pd.read_csv(adir / "genes_with_compartment.csv")
def dist_to_gap(chrom, pos):
    p = gpos.get(chrom)
    return np.inf if p is None else np.abs(p - pos).min()
g["gap_kb"] = [dist_to_gap(c, m) / 1000 for c, m in zip(g.chrom, g.mid)]

rows = []
for F in flanks:
    keep = g[g.gap_kb > F].drop(columns="gap_kb")
    d = odir / f"mask_{F}kb"
    d.mkdir(parents=True, exist_ok=True)
    keep.to_csv(d / "genes_with_compartment.csv", index=False)
    shutil.copy(adir / "phasing.tsv", d / "phasing.tsv")
    called = g.comp.isin(["A", "B"])
    rows.append({"flank_kb": F, "genes_removed": int((g.gap_kb <= F).sum()),
                 "genes_kept_with_call": int(keep.comp.isin(["A", "B"]).sum()),
                 "pct_called_genes_removed": round(100 * ((g.gap_kb <= F) & called).sum() / called.sum(), 1)})
print(pd.DataFrame(rows).to_string(index=False))

# compartment transitions vs gaps
v = pd.read_csv(adir / "comp.cis.vecs.tsv", sep="\t")
v = v[v.chrom.str.match(r"chr\d\d$")].dropna(subset=["E1"]).sort_values(["chrom", "start"])
obs = exp = n = 0
sims = []
trans_by_chrom = {}
for c, d in v.groupby("chrom"):
    s = np.sign(d.E1.values)
    edge = (d.end.values[:-1])[s[:-1] != s[1:]]          # transition position = bin edge
    trans_by_chrom[c] = (edge, d.end.values[:-1])
for F in flanks:
    near = tot = 0
    sim_near = np.zeros(2000)
    for c, (edge, cand) in trans_by_chrom.items():
        p = gpos.get(c, np.array([-1e12]))
        near += sum(np.abs(p - e).min() <= F * 1000 for e in edge)
        tot += len(edge)
        cand_near = np.array([np.abs(p - e).min() <= F * 1000 for e in cand])
        k = len(edge)
        if k:
            sim_near += np.array([cand_near[rng.choice(len(cand), k, replace=False)].sum() for _ in range(2000)])
    p_val = (1 + np.sum(sim_near >= near)) / 2001
    print(f"transitions within {F} kb of a gap: {near}/{tot} ({100*near/tot:.1f}%); "
          f"expected at random bin edges {sim_near.mean():.1f} ({100*sim_near.mean()/tot:.1f}%); one-sided P = {p_val:.3f}")
