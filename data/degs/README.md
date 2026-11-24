# Salt-responsive DEGs (Natarajan et al., 2021)

- `Salt_DEGs_31July2020.xlsx`: DEG workbook from the 2021 study (leaf and root sheets plus up/down and annotation subsets).
- `leaf_DEGs.csv` (1,759 genes) and `root_DEGs.csv` (1,487 genes): the first two sheets, with fold changes at 24 h, 48 h and 14 d
  (log2 fold change of at least 1, FDR below 0.05).

The 3,246 salinity-responsive genes reported in Natarajan et al. (2021) are the sum 1,759 + 1,487. Because 328 genes are
differentially expressed in both organs, the number of unique genes is 2,918. These 2,918 genes are `genesets/salinome.txt`,
and the leaf and root lists are `genesets/deg_leaf.txt` and `genesets/deg_root.txt`.
