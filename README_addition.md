## Reviewer-response analyses (v1.1)
- `scripts/hic_junction_mask.py`: treats all assembly gaps (NCBI GCA_019155195.1 gap report) as possible scaffolding junctions; tests whether A/B transitions cluster near gaps and writes gene tables with genes within 100 or 250 kb of a gap removed.
  `python scripts/hic_junction_mask.py results GCA_019155195.1_amgenome_v1.0_genomic_gaps.txt.gz GCA_019155195.1_amgenome_v1.0_assembly_report.txt jm 100 250`
- `scripts/hic_robustness.py`: chromosome x density x expression stratified permutation, circular-shift block permutation (raw and covariate-adjusted), Cochran-Mantel-Haenszel test by chromosome, leaf-only and root-only subsets.
  `python scripts/hic_robustness.py results genesets data/Counttable_avisalt.txt results/robustness 10000`
- `results/robustness/`: outputs for all genes and with 100 kb and 250 kb masking around assembly gaps.
- `data/degs/`: leaf (1,759) and root (1,487) DEG tables from Natarajan et al. (2021). Their sum is the published 3,246; the 2,918 unique genes form `genesets/salinome.txt`.
