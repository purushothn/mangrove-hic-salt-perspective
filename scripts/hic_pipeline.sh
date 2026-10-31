#!/usr/bin/env bash
# Hi-C processing for Avicennia marina (SRR12126768, Dovetail Hi-C, DpnII, leaf).
# Usage: bash hic_pipeline.sh GENOME.fa OUTDIR [THREADS] [TMPDIR]
# GENOME.fa = AmGA_v1.0 chromosome-scale assembly (the 31 chromosomes plus any small scaffolds).
# Needs: sra-tools, pigz, bwa, samtools, pairtools, cooler, cooltools.
# Suggested install: mamba create -n hic -c conda-forge -c bioconda sra-tools pigz bwa samtools pairtools cooler cooltools
set -euo pipefail

GENOME=$1
OUT=$2
T=${3:-16}
TMP=${4:-$OUT/tmp}   # optional 4th argument: temp dir for sorting, e.g. /dev/shm/hic_tmp (RAM disk, much faster than the HDD)
RUN=SRR12126768
NCHR=31   # chromosomes chr01-chr31 in AmGA_v1.0 (the 'un' scaffold is excluded from the analysis view)

mkdir -p "$OUT"/{fastq,pairs,cool,tmp,logs} "$TMP"

# 1. Download reads (about 14 Gb SRA, 133M pairs, 2 x 150 bp)
if [ ! -s "$OUT/fastq/${RUN}_1.fastq.gz" ]; then
  prefetch "$RUN" -O "$OUT/fastq" 2> "$OUT/logs/prefetch.log"
  fasterq-dump "$OUT/fastq/$RUN" -O "$OUT/fastq" -e "$T" --split-files -t "$OUT/tmp" 2>> "$OUT/logs/prefetch.log"
  pigz -p "$T" "$OUT"/fastq/${RUN}_1.fastq "$OUT"/fastq/${RUN}_2.fastq
fi

# 2. Index genome, chromosome sizes, and a view of the chromosome-scale scaffolds
[ -s "$GENOME.bwt" ] || bwa index "$GENOME" 2> "$OUT/logs/bwa_index.log"
samtools faidx "$GENOME"
cut -f1,2 "$GENOME.fai" > "$OUT/cool/chrom.sizes"
awk 'BEGIN{OFS="\t"} $1 ~ /^chr[0-9]+$/ {print $1,0,$2,$1}' "$GENOME.fai" > "$OUT/cool/view.tsv"
[ "$(wc -l < "$OUT/cool/view.tsv")" -eq "$NCHR" ] || { echo "Expected $NCHR chromosomes named chrNN, found $(wc -l < "$OUT/cool/view.tsv")"; exit 1; }
echo "Chromosomes used (chr01-chr31; 'un' and other scaffolds are excluded from the analysis):"
cut -f1,3 "$OUT/cool/view.tsv"

# 3. Align and build valid pairs (standard Open2C workflow; DpnII, so no restriction-site file is needed)
bwa mem -5SP -t "$T" "$GENOME" "$OUT/fastq/${RUN}_1.fastq.gz" "$OUT/fastq/${RUN}_2.fastq.gz" 2> "$OUT/logs/bwa_mem.log" \
 | pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
     --nproc-in "$T" --nproc-out "$T" --chroms-path "$OUT/cool/chrom.sizes" \
 | pairtools sort --nproc "$T" --tmpdir "$TMP" --memory 100G \
 | pairtools dedup --nproc-in "$T" --nproc-out "$T" --mark-dups --output-stats "$OUT/pairs/dedup_stats.txt" \
 | pairtools split --nproc-in "$T" --nproc-out "$T" --output-pairs "$OUT/pairs/valid.pairs.gz"

# Keep unique, unambiguous pairs only (drops duplicates marked DD and partly mapped types)
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' --nproc-in "$T" --nproc-out "$T" -o "$OUT/pairs/filtered.pairs.gz" "$OUT/pairs/valid.pairs.gz"

# 4. Contact matrices at several resolutions, balanced (ICE)
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 "$OUT/cool/chrom.sizes:10000" "$OUT/pairs/filtered.pairs.gz" "$OUT/cool/am.10kb.cool"
cooler zoomify -r 10000,20000,50000,100000,250000 --balance -p "$T" -o "$OUT/cool/am.mcool" "$OUT/cool/am.10kb.cool"

echo "Done. Contact matrix: $OUT/cool/am.mcool ; view file: $OUT/cool/view.tsv ; QC: $OUT/pairs/dedup_stats.txt"
