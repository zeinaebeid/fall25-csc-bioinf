#!/usr/bin/env bash
set -euo pipefail
[ -n "${TRACE:-}" ] && set -x

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
CODE="$ROOT/week1/code"
DATA="$ROOT/week1/data"
OUT="$ROOT/week1/out"
mkdir -p "$OUT"

DATASETS=(data1 data2 data3 data4)

fmt_time() {
  local s="$1"
  printf "%d:%02d:%02d" "$((s/3600))" "$(((s%3600)/60))" "$((s%60))"
}

n50_from_fasta() {
  local fa="$1"
  awk '/^>/ {if (len) print len; len=0; next} {len+=length($0)} END {if (len) print len}' "$fa" \
  | sort -nr \
  | awk '{sum+=$1; a[NR]=$1} END {half=sum/2; s=0; for(i=1;i<=NR;i++){s+=a[i]; if(s>=half){print a[i]; exit}}}'
}

# Print header with fixed widths
printf "%-10s %-10s %-10s %-10s\n" "Dataset" "Language" "Runtime" "N50"
printf "%-10s %-10s %-10s %-10s\n" "--------" "--------" "-------" "----"

pushd "$CODE" >/dev/null

for d in "${DATASETS[@]}"; do
  (
    fasta_dir="$DATA/$d"

    # --- Python run ---
    t0=$(date +%s)
    python3 "$CODE/main.py" "$fasta_dir" >"$OUT/${d}_python.log" 2>&1
    t1=$(date +%s)
    py_rt=$((t1 - t0))
    mv -f "$fasta_dir/contig.fasta" "$OUT/${d}.python.fasta"
    py_n50="$(n50_from_fasta "$OUT/${d}.python.fasta")"
    printf "%-10s %-10s %-10s %-10s\n" "$d" "python" "$(fmt_time "$py_rt")" "${py_n50:-NA}"

    # --- Codon run ---
    t2=$(date +%s)
    if codon run -release "$CODE/main_codon.py" "$fasta_dir" >"$OUT/${d}_codon.log" 2>&1; then
      t3=$(date +%s)
      co_rt=$((t3 - t2))
      mv -f "$fasta_dir/contig.fasta" "$OUT/${d}.codon.fasta"
      co_n50="$(n50_from_fasta "$OUT/${d}.codon.fasta")"
      printf "%-10s %-10s %-10s %-10s\n" "$d" "codon" "$(fmt_time "$co_rt")" "${co_n50:-NA}"
    else
      printf "%-10s %-10s %-10s %-10s\n" "$d" "codon" "NA" "NA"
    fi
  )
  # pause between datasets so memory is reclaimed
  sleep 5
done

popd >/dev/null
