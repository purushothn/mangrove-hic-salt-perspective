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

## Citation
Natarajan P. (2026). Scripts and gene sets for 3D genome organization in *Avicennia marina*. Zenodo. DOI to be added.

Contact: pnatarajan@umes.edu
