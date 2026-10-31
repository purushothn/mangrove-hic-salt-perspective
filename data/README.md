# Input data

- `Counttable_avisalt.txt`: pooled RNA-seq count table for the A. marina salt experiment (one column per condition, no replicates), from Natarajan et al. (2021). Used by `hic_followup.py` for expression matching. The Leaf_28d and Root_28d columns are present in the file but were excluded from the analysis (28 d not part of the design; Root_28d is a low-depth library).

Not stored here (public):
- Hi-C reads: NCBI SRA run SRR12126768 (BioProject PRJNA392013); downloaded by `scripts/hic_pipeline.sh`.
- Genome assembly and annotation (AmGA_v1.0; accession JACDXK000000000), Natarajan et al. (2021).
