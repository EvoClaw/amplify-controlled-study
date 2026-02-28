#!/bin/bash
set -euo pipefail

OUTDIR="/home/yanlin/comp/amplify/results/full"
SAMPINFO="/home/yanlin/comp/amplify/1000GP/samples.info"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$LOGDIR"

echo "=== Continuing pipeline: INDEL and SV PLINK conversion ==="
echo "Start: $(date)"

for vtype in indel sv; do
    prefix="${OUTDIR}/${vtype}/all_autosomes"
    merged="${prefix}.vcf.gz"
    [ ! -f "$merged" ] && continue

    if [ ! -f "${prefix}.pgen" ]; then
        plink2 --vcf "$merged" --make-pgen --out "$prefix" \
            --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 300 missing \
            --allow-extra-chr --threads 32 --memory 64000
    fi
    echo "  ${vtype}: PLINK conversion done"
done

echo ""
echo "=== LD pruning (SNV and INDEL) ==="
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
    npruned=$(wc -l < "${prefix}_ldprune.prune.in")
    echo "  ${vtype}: ${npruned} variants after LD pruning"
done

echo ""
echo "=== Running PCA ==="
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
    echo "  Per-pop frequencies ${vtype}: done"
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
echo "=== Pipeline continuation complete ==="
echo "Finished: $(date)"
