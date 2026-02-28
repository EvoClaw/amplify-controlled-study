#!/bin/bash
set -euo pipefail

DATADIR="/home/yanlin/comp/amplify/1000GP/20220422_3202_phased_SNV_INDEL_SV"
OUTDIR="/home/yanlin/comp/amplify/results/full"
SAMPLES="/home/yanlin/comp/amplify/results/unrelated_samples.txt"
NCPU=64

mkdir -p "$OUTDIR"/{snv,indel,sv}

pop_to_superpop() {
    case "$1" in
        YRI|LWK|GWD|MSL|ESN|ACB|ASW) echo "AFR" ;;
        CEU|TSI|FIN|GBR|IBS) echo "EUR" ;;
        CHB|JPT|CHS|CDX|KHV) echo "EAS" ;;
        GIH|PJL|BEB|STU|ITU) echo "SAS" ;;
        MXL|PUR|CLM|PEL) echo "AMR" ;;
        *) echo "UNK" ;;
    esac
}

SAMFILE="/home/yanlin/comp/amplify/1000GP/samples.info"
> "$OUTDIR/sample_pop_map.tsv"
while IFS=$'\t' read -r sid pop sex; do
    sp=$(pop_to_superpop "$pop")
    echo -e "${sid}\t${pop}\t${sp}\t${sex}"
done < "$SAMFILE" > "$OUTDIR/sample_pop_map_all.tsv"

comm -12 <(sort "$SAMPLES") <(cut -f1 "$OUTDIR/sample_pop_map_all.tsv" | sort) > "$OUTDIR/unrelated_in_vcf.txt"
grep -F -f "$OUTDIR/unrelated_in_vcf.txt" "$OUTDIR/sample_pop_map_all.tsv" > "$OUTDIR/sample_pop_map.tsv"

echo "Unrelated samples in VCF: $(wc -l < "$OUTDIR/unrelated_in_vcf.txt")"
echo "Population counts:"
cut -f2 "$OUTDIR/sample_pop_map.tsv" | sort | uniq -c | sort -rn

process_chr() {
    local chr=$1
    local vcf="${DATADIR}/1kGP_high_coverage_Illumina.${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    if [ "$chr" = "chrX" ]; then
        vcf="${DATADIR}/1kGP_high_coverage_Illumina.${chr}.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    fi

    echo "[$(date +%H:%M:%S)] Processing $chr..."

    # SNVs: biallelic, FILTER=PASS, MAF>=0.01, missingness<5%
    bcftools view -S "$OUTDIR/unrelated_in_vcf.txt" -v snps -m2 -M2 \
        --min-af 0.01:minor -i 'F_MISSING<0.05' \
        "$vcf" -Oz -o "$OUTDIR/snv/${chr}.snv.vcf.gz" --threads 4 2>/dev/null
    bcftools index "$OUTDIR/snv/${chr}.snv.vcf.gz" --threads 4 2>/dev/null

    # INDELs: biallelic, FILTER=PASS, MAF>=0.01, missingness<5%
    bcftools view -S "$OUTDIR/unrelated_in_vcf.txt" -v indels -m2 -M2 \
        --min-af 0.01:minor -i 'F_MISSING<0.05' \
        "$vcf" -Oz -o "$OUTDIR/indel/${chr}.indel.vcf.gz" --threads 4 2>/dev/null
    bcftools index "$OUTDIR/indel/${chr}.indel.vcf.gz" --threads 4 2>/dev/null

    # SVs: biallelic, MAF>=0.005, missingness<10%
    bcftools view -S "$OUTDIR/unrelated_in_vcf.txt" -m2 -M2 \
        --min-af 0.005:minor -i 'F_MISSING<0.10 && STRLEN(REF)>=50 | STRLEN(ALT)>=50 | INFO/SVTYPE!="."' \
        "$vcf" -Oz -o "$OUTDIR/sv/${chr}.sv.vcf.gz" --threads 4 2>/dev/null
    bcftools index "$OUTDIR/sv/${chr}.sv.vcf.gz" --threads 4 2>/dev/null

    local nsnv=$(bcftools view -H "$OUTDIR/snv/${chr}.snv.vcf.gz" 2>/dev/null | wc -l)
    local nindel=$(bcftools view -H "$OUTDIR/indel/${chr}.indel.vcf.gz" 2>/dev/null | wc -l)
    local nsv=$(bcftools view -H "$OUTDIR/sv/${chr}.sv.vcf.gz" 2>/dev/null | wc -l)
    echo "[$(date +%H:%M:%S)] $chr: SNV=$nsnv, INDEL=$nindel, SV=$nsv"
}

export -f process_chr pop_to_superpop
export DATADIR OUTDIR SAMPLES

echo "=== Starting chr22 pilot ==="
process_chr chr22
echo "=== Pilot complete ==="
