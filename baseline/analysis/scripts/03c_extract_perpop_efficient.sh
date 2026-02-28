#!/bin/bash
# Efficient per-superpopulation indel extraction - one chromosome at a time per pop
VCF_DIR="/home/yanlin/comp/cursor/1000GP/20220422_3202_phased_SNV_INDEL_SV"
DATA_DIR="/home/yanlin/comp/cursor/analysis/data"
OUT_DIR="$DATA_DIR/indels_perpop"
mkdir -p "$OUT_DIR"

SUPERPOPS="AFR AMR EAS EUR SAS"

for CHR in chr{2..22}; do
    for SPOP in $SUPERPOPS; do
        OUT="$OUT_DIR/${CHR}_${SPOP}.tsv.gz"
        [ -s "$OUT" ] && continue
        
        VCF="$VCF_DIR/1kGP_high_coverage_Illumina.${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        SAMPLES="$DATA_DIR/samples_${SPOP}.txt"
        
        bcftools view -v indels -S "$SAMPLES" "$VCF" 2>/dev/null | \
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' 2>/dev/null | \
            gzip > "$OUT" &
        
        if (( $(jobs -r | wc -l) >= 8 )); then
            wait -n
        fi
    done
    wait
    echo "DONE: $CHR (all pops)"
done

echo "All extractions complete"
