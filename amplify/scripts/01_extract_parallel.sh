#!/bin/bash
set -euo pipefail

DATADIR="/home/yanlin/comp/amplify/1000GP/20220422_3202_phased_SNV_INDEL_SV"
OUTDIR="/home/yanlin/comp/amplify/results/full"
SAMPINFO="/home/yanlin/comp/amplify/1000GP/samples.info"
LOGDIR="${OUTDIR}/logs"

mkdir -p "${OUTDIR}/snv" "${OUTDIR}/indel" "${OUTDIR}/sv" "${LOGDIR}"

process_chr() {
    local chr=$1
    local vcf="${DATADIR}/1kGP_high_coverage_Illumina.${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

    if [ ! -f "$vcf" ]; then
        echo "SKIP: $chr"
        return
    fi

    local snv_out="${OUTDIR}/snv/${chr}.snv.vcf.gz"
    local indel_out="${OUTDIR}/indel/${chr}.indel.vcf.gz"
    local sv_out="${OUTDIR}/sv/${chr}.sv.vcf.gz"

    if [ ! -f "${snv_out}.csi" ]; then
        bcftools view -v snps -m2 -M2 --min-af 0.01:minor -i 'F_MISSING<0.05' \
            "$vcf" -Oz -o "$snv_out" --threads 4
        bcftools index "$snv_out" --threads 2
    fi

    if [ ! -f "${indel_out}.csi" ]; then
        bcftools view -v indels -m2 -M2 --min-af 0.01:minor -i 'F_MISSING<0.05' \
            "$vcf" -Oz -o "$indel_out" --threads 4
        bcftools index "$indel_out" --threads 2
    fi

    if [ ! -f "${sv_out}.csi" ]; then
        bcftools view -m2 -M2 --min-af 0.005:minor \
            -i 'F_MISSING<0.10 && (STRLEN(REF)>=50 || STRLEN(ALT)>=50)' \
            "$vcf" -Oz -o "$sv_out" --threads 4
        bcftools index "$sv_out" --threads 2
    fi

    echo "DONE: $chr"
}

echo "=== Parallel variant extraction (22 autosomes) ==="
echo "Start: $(date)"

MAXJOBS=11
running=0

for i in $(seq 1 22); do
    chr="chr${i}"
    process_chr "$chr" > "${LOGDIR}/${chr}.log" 2>&1 &
    running=$((running + 1))
    if [ $running -ge $MAXJOBS ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait
echo "=== All chromosomes extracted ==="
echo "$(date)"

for chr in chr{1..22}; do
    cat "${LOGDIR}/${chr}.log" 2>/dev/null
done

echo ""
echo "=== Counting variants per type ==="
for vtype in snv indel sv; do
    total=0
    for f in "${OUTDIR}/${vtype}"/chr*.${vtype}.vcf.gz; do
        [ -f "$f" ] || continue
        n=$(bcftools index -n "$f" 2>/dev/null || echo 0)
        total=$((total + n))
    done
    echo "  ${vtype}: ${total} total variants"
done

echo ""
echo "=== Merging per-type VCFs ==="
for vtype in snv indel sv; do
    merged="${OUTDIR}/${vtype}/all_autosomes.vcf.gz"
    if [ -f "${merged}.csi" ]; then
        echo "  ${vtype}: already merged"
        continue
    fi
    rm -f "$merged" "${merged}.csi"

    ls "${OUTDIR}/${vtype}"/chr*.${vtype}.vcf.gz 2>/dev/null | sort -V > "${OUTDIR}/${vtype}/vcf_list.txt"
    nfiles=$(wc -l < "${OUTDIR}/${vtype}/vcf_list.txt")
    if [ "$nfiles" -eq 0 ]; then
        echo "  ${vtype}: no files to merge"
        continue
    fi

    bcftools concat -f "${OUTDIR}/${vtype}/vcf_list.txt" -Oz -o "$merged" --threads 16
    bcftools index "$merged" --threads 4
    echo "  ${vtype}: merged"
done

echo ""
echo "=== Converting to PLINK format ==="
for vtype in snv indel sv; do
    prefix="${OUTDIR}/${vtype}/all_autosomes"
    merged="${prefix}.vcf.gz"
    [ ! -f "$merged" ] && continue

    if [ ! -f "${prefix}.pgen" ]; then
        plink2 --vcf "$merged" --make-pgen --out "$prefix" \
            --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 20 \
            --allow-extra-chr --threads 32 --memory 64000
    fi
    echo "  ${vtype}: PLINK conversion done"
done

echo ""
echo "=== LD pruning (SNV and INDEL only) ==="
for vtype in snv indel; do
    prefix="${OUTDIR}/${vtype}/all_autosomes"
    [ ! -f "${prefix}.pgen" ] && continue

    if [ ! -f "${prefix}_ldprune.prune.in" ]; then
        plink2 --pfile "$prefix" --indep-pairwise 50 5 0.2 \
            --out "${prefix}_ldprune" --allow-extra-chr --threads 32
    fi
    if [ ! -f "${prefix}_pruned.pgen" ]; then
        plink2 --pfile "$prefix" --extract "${prefix}_ldprune.prune.in" \
            --make-pgen --out "${prefix}_pruned" --allow-extra-chr --threads 32
    fi
    echo "  ${vtype}: LD pruning done"
done

echo ""
echo "=== Running PCA (per variant type) ==="
for vtype in snv indel sv; do
    if [ "$vtype" = "sv" ]; then
        prefix="${OUTDIR}/${vtype}/all_autosomes"
    else
        prefix="${OUTDIR}/${vtype}/all_autosomes_pruned"
        [ ! -f "${prefix}.pgen" ] && prefix="${OUTDIR}/${vtype}/all_autosomes"
    fi
    [ ! -f "${prefix}.pgen" ] && continue

    if [ ! -f "${prefix}_pca.eigenvec" ]; then
        plink2 --pfile "$prefix" --pca 20 --out "${prefix}_pca" \
            --allow-extra-chr --threads 32
    fi
    echo "  PCA ${vtype}: done"
done

echo ""
echo "=== Computing per-population allele frequencies ==="
POPS=$(cut -f2 "$SAMPINFO" | sort -u)

for vtype in snv indel sv; do
    prefix="${OUTDIR}/${vtype}/all_autosomes"
    [ ! -f "${prefix}.pgen" ] && continue

    running=0
    for pop in $POPS; do
        popfile="${OUTDIR}/${vtype}/pop_${pop}.txt"
        freqout="${OUTDIR}/${vtype}/freq_${pop}"
        [ -f "${freqout}.afreq" ] && continue

        grep -w "$pop" "$SAMPINFO" | awk '{print $1"\t"$1}' > "$popfile"
        n=$(wc -l < "$popfile")
        if [ "$n" -ge 5 ]; then
            plink2 --pfile "$prefix" --keep "$popfile" --freq \
                --out "$freqout" --allow-extra-chr --threads 4 \
                > "${LOGDIR}/freq_${vtype}_${pop}.log" 2>&1 &
            running=$((running + 1))
            if [ $running -ge 8 ]; then
                wait -n
                running=$((running - 1))
            fi
        fi
    done
    wait
    echo "  Frequencies ${vtype}: done"
done

echo ""
echo "=== Computing super-population SFS frequencies ==="
declare -A SPOP_POPS
SPOP_POPS[AFR]="YRI LWK GWD MSL ESN ACB ASW"
SPOP_POPS[EUR]="CEU TSI FIN GBR IBS"
SPOP_POPS[EAS]="CHB JPT CHS CDX KHV"
SPOP_POPS[SAS]="GIH PJL BEB STU ITU"
SPOP_POPS[AMR]="MXL PUR CLM PEL"

for vtype in snv indel sv; do
    prefix="${OUTDIR}/${vtype}/all_autosomes"
    [ ! -f "${prefix}.pgen" ] && continue

    for spop in AFR EUR EAS SAS AMR; do
        popfile="${OUTDIR}/${vtype}/spop_${spop}.txt"
        freqout="${OUTDIR}/${vtype}/sfs_${spop}"
        [ -f "${freqout}.afreq" ] && continue

        > "$popfile"
        for p in ${SPOP_POPS[$spop]}; do
            grep -w "$p" "$SAMPINFO" | awk '{print $1"\t"$1}' >> "$popfile"
        done
        plink2 --pfile "$prefix" --keep "$popfile" --freq \
            --out "$freqout" --allow-extra-chr --threads 8 \
            > "${LOGDIR}/sfs_${vtype}_${spop}.log" 2>&1 &
    done
    wait
    echo "  SFS frequencies ${vtype}: done"
done

echo ""
echo "=== All preprocessing complete ==="
echo "Finished: $(date)"
