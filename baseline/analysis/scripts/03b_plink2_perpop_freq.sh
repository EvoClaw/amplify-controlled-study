#!/bin/bash
# Use plink2 to compute per-superpopulation allele frequencies for indels
# Much faster than bcftools per-pop subsetting

VCF_DIR="/home/yanlin/comp/cursor/1000GP/20220422_3202_phased_SNV_INDEL_SV"
DATA_DIR="/home/yanlin/comp/cursor/analysis/data"
OUT_DIR="/home/yanlin/comp/cursor/analysis/data/plink_freq"
mkdir -p "$OUT_DIR"

CLUSTER_FILE="$DATA_DIR/superpop_clusters.txt"

process_chr() {
    CHR=$1
    VCF="$VCF_DIR/1kGP_high_coverage_Illumina.${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    
    if [ -f "$OUT_DIR/${CHR}.afreq.varsort.tsv.gz" ]; then
        echo "SKIP: $CHR"
        return
    fi
    
    # Use plink2 to compute per-cluster frequencies for indels only
    # --snps-only 'just-acgt' would skip indels, so we DON'T use it
    plink2 --vcf "$VCF" \
        --freq \
        --within "$CLUSTER_FILE" \
        --out "$OUT_DIR/${CHR}" \
        --threads 2 \
        --memory 4000 \
        2>/dev/null
    
    echo "DONE: $CHR"
}

for CHR in chr{1..22}; do
    process_chr "$CHR" &
    if (( $(jobs -r | wc -l) >= 4 )); then
        wait -n
    fi
done
wait

echo "All chromosomes done"
