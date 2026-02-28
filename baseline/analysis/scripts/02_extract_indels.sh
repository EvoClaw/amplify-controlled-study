#!/bin/bash
# Extract indels from VCF files for all chromosomes
# Output: per-chromosome TSV with indel info and genotype counts per super-population

VCF_DIR="/home/yanlin/comp/cursor/1000GP/20220422_3202_phased_SNV_INDEL_SV"
DATA_DIR="/home/yanlin/comp/cursor/analysis/data"
OUT_DIR="/home/yanlin/comp/cursor/analysis/data/indels"
mkdir -p "$OUT_DIR"

extract_chr() {
    CHR=$1
    VCF="$VCF_DIR/1kGP_high_coverage_Illumina.${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    OUT="$OUT_DIR/${CHR}_indels.tsv.gz"
    
    if [ -f "$OUT" ]; then
        echo "SKIP: $CHR (already exists)"
        return
    fi
    
    # Extract indels only (TYPE="INDEL"), output: CHROM, POS, REF, ALT, AC, AN, AF
    bcftools view -v indels "$VCF" 2>/dev/null | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\n' 2>/dev/null | \
        gzip > "$OUT"
    
    echo "DONE: $CHR"
}

# Run chromosomes in parallel (4 at a time to avoid I/O bottleneck)
for CHR in chr{1..22}; do
    extract_chr "$CHR" &
    # Limit parallelism
    if (( $(jobs -r | wc -l) >= 6 )); then
        wait -n
    fi
done
wait

echo "All chromosomes done"
