#!/bin/bash
# Extract per-superpopulation indel allele counts from VCF files

VCF_DIR="/home/yanlin/comp/cursor/1000GP/20220422_3202_phased_SNV_INDEL_SV"
DATA_DIR="/home/yanlin/comp/cursor/analysis/data"
OUT_DIR="/home/yanlin/comp/cursor/analysis/data/indels_perpop"
mkdir -p "$OUT_DIR"

SUPERPOPS="AFR AMR EAS EUR SAS"

extract_job() {
    SPOP=$1
    CHR=$2
    VCF="$VCF_DIR/1kGP_high_coverage_Illumina.${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    SAMPLES="$DATA_DIR/samples_${SPOP}.txt"
    OUT="$OUT_DIR/${CHR}_${SPOP}.tsv.gz"
    
    if [ -f "$OUT" ]; then
        return
    fi
    
    bcftools view -v indels -S "$SAMPLES" "$VCF" 2>/dev/null | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' 2>/dev/null | \
        gzip > "$OUT"
    
    echo "DONE: ${CHR}_${SPOP}"
}

for CHR in chr{1..22}; do
    for SPOP in $SUPERPOPS; do
        extract_job "$SPOP" "$CHR" &
        if (( $(jobs -r | wc -l) >= 8 )); then
            wait -n
        fi
    done
done
wait

echo "All per-population extractions done"
