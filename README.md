# mangrove-hic-salt-perspective

Scripts, gene sets and results supporting the Perspective "Three-dimensional genome organization: the missing regulatory layer in mangrove salt adaptation" (Natarajan, Frontiers in Forests and Global Change, Research Topic "From Genes to Forests").

The analysis reuses the public *Avicennia marina* leaf Hi-C library (SRA run SRR12126768, BioProject PRJNA392013) and the reference assembly and annotation described in Natarajan et al. (2021), Communications Biology 4:851.

## Layout
- `scripts/hic_pipeline.sh`: download, alignment (bwa mem), pairtools processing, cooler matrix and multi-resolution .mcool.
- `scripts/hic_analysis.py`: A/B compartments (100 kb, cooltools eigs_cis phased by gene density), insulation and boundaries (50 kb), gene-set tests and the chromosome figure.
- `scripts/hic_followup.py`: permutation tests matched on local gene density and mean expression (10,000 permutations).
- `genesets/`: gene ID lists (salinome, leaf and root subsets, HKT genes, literature homologs, expressed-gene background).
- `results/`: output tables from the analysis (gene-set tests, compartment and insulation tracks, per-gene compartment calls).
- `figures/`: Figure 1 (chromosome 19 Hi-C matrix and E1 track).
- `data/`: RNA-seq count table used for expression matching.
- `environment.yml`: conda environment (cooltools 0.7.1).

## Run order
```
bash scripts/hic_pipeline.sh GENOME.fa OUTDIR [THREADS] [TMPDIR]
python scripts/hic_analysis.py AM.mcool VIEW.tsv ANNOTATION.gff3 genesets OUTDIR [HIGHLIGHT_CHROM]
python scripts/hic_followup.py OUTDIR genesets [N_PERM] [COUNT_TABLE]
```
Inputs: the raw reads from SRA, the genome FASTA and the GFF3 annotation, and the RNA-seq count table from Natarajan et al. (2021) for expression matching.


## Reviewer-response analyses (v1.1)
- `scripts/hic_junction_mask.py`: treats all assembly gaps (NCBI GCA_019155195.1 gap report) as possible scaffolding junctions; tests whether A/B transitions cluster near gaps and writes gene tables with genes within 100 or 250 kb of a gap removed.
  `python scripts/hic_junction_mask.py results GCA_019155195.1_amgenome_v1.0_genomic_gaps.txt.gz GCA_019155195.1_amgenome_v1.0_assembly_report.txt jm 100 250`
- `scripts/hic_robustness.py`: chromosome x density x expression stratified permutation, circular-shift block permutation (raw and covariate-adjusted), Cochran-Mantel-Haenszel test by chromosome, leaf-only and root-only subsets.
  `python scripts/hic_robustness.py results genesets data/Counttable_avisalt.txt results/robustness 10000`
- `results/robustness/`: outputs for all genes and with 100 kb and 250 kb masking around assembly gaps.
- `data/degs/`: leaf (1,759) and root (1,487) DEG tables from Natarajan et al. (2021). Their sum is the published 3,246; the 2,918 unique genes form `genesets/salinome.txt`.


## Citation
Natarajan P. (2026). Scripts and gene sets for 3D genome organization in *Avicennia marina*. Zenodo. DOI:10.5281/zenodo.23825794

Contact: pnatarajan@umes.edu
