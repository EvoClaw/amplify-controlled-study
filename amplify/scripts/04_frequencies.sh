#!/bin/bash
set -uo pipefail

OUTDIR="/home/yanlin/comp/amplify/results/full"
SAMPINFO="/home/yanlin/comp/amplify/1000GP/samples.info"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$LOGDIR"

echo "=== Computing per-population allele frequencies ==="
echo "Start: $(date)"

POPS=$(cut -f2 "$SAMPINFO" | sort -u)

for vtype in snv indel sv; do
    prefix="${OUTDIR}/${vtype}/all_autosomes"
    [ ! -f "${prefix}.pgen" ] && continue

    pids=()
    for pop in $POPS; do
        popfile="${OUTDIR}/${vtype}/pop_${pop}.txt"
        freqout="${OUTDIR}/${vtype}/freq_${pop}"

        grep -w "$pop" "$SAMPINFO" | awk '{print $1}' > "$popfile"
        n=$(wc -l < "$popfile")
        if [ "$n" -ge 5 ]; then
            plink2 --pfile "$prefix" --keep "$popfile" --freq \
                --out "$freqout" --allow-extra-chr --threads 4 \
                > "${LOGDIR}/freq_${vtype}_${pop}.log" 2>&1 &
            pids+=($!)
            if [ ${#pids[@]} -ge 8 ]; then
                wait "${pids[0]}"
                pids=("${pids[@]:1}")
            fi
        fi
    done
    for pid in "${pids[@]}"; do
        wait "$pid" || true
    done
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

    pids=()
    for spop in AFR EUR EAS SAS AMR; do
        popfile="${OUTDIR}/${vtype}/spop_${spop}.txt"
        freqout="${OUTDIR}/${vtype}/sfs_${spop}"

        > "$popfile"
        for p in ${SPOP_POPS[$spop]}; do
            grep -w "$p" "$SAMPINFO" | awk '{print $1}' >> "$popfile"
        done
        plink2 --pfile "$prefix" --keep "$popfile" --freq \
            --out "$freqout" --allow-extra-chr --threads 8 \
            > "${LOGDIR}/sfs_${vtype}_${spop}.log" 2>&1 &
        pids+=($!)
    done
    for pid in "${pids[@]}"; do
        wait "$pid" || true
    done
    echo "  SFS frequencies ${vtype}: done"
done

echo ""
echo "=== Frequency computation complete ==="
echo "Finished: $(date)"

echo ""
echo "=== Verification ==="
for vtype in snv indel sv; do
    n=$(ls "${OUTDIR}/${vtype}"/freq_*.afreq 2>/dev/null | wc -l)
    s=$(ls "${OUTDIR}/${vtype}"/sfs_*.afreq 2>/dev/null | wc -l)
    echo "  ${vtype}: ${n} pop freq files, ${s} SFS freq files"
done
