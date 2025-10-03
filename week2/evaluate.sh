#!/usr/bin/env bash
set -euo pipefail

echo "🚀 Week 2 — TRViz (Python + Codon) evaluation"

echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CODE_DIR="code"

echo "📁 Checking files…"
need=(
  "$CODE_DIR/decomposer.py"
  "$CODE_DIR/motif_encoder.py"
  "$CODE_DIR/motif_aligner.py"
  "$CODE_DIR/main.py"
  "$CODE_DIR/test_python.py"
  "$CODE_DIR/test_codon.py"
)
for f in "${need[@]}"; do
  [[ -f "$f" ]] && echo "  ✅ $f" || { echo "  ❌ missing: $f"; exit 1; }
done

export PYTHONPATH="${PYTHONPATH:-}:$CODE_DIR"

# -------------------------
# 1) Python tests (foreground)
# -------------------------
echo ""
echo "🐍 Python tests"
PY_STATUS="SKIP"
if command -v python3 >/dev/null 2>&1; then
  echo "▶ Running Python tests…"
  if TRVIZ_VERBOSE="${TRVIZ_VERBOSE:-0}" python3 "$CODE_DIR/test_python.py"; then
    PY_STATUS="PASS"
    echo "✅ Python tests completed"
  else
    PY_STATUS="FAIL"
    echo "❌ Python tests failed"
  fi
else
  echo "❌ python3 not found."
  PY_STATUS="SKIP"
fi

# -------------------------
# 2) Codon tests (fully silenced backtrace)
#    We run Codon inside an inner bash so any SIGABRT message
#    is emitted by the inner shell and captured in a tempfile.
# -------------------------
echo ""
echo "🧪 Codon tests"
CO_STATUS="SKIP"
if command -v codon >/dev/null 2>&1; then
  echo "▶ Running Codon tests…"
  tmp_out="$(mktemp)"
  set +e
  # Inner shell runs codon; all output captured
  bash -c 'TRVIZ_VERBOSE="${TRVIZ_VERBOSE:-0}" codon run "$1"' bash "$CODE_DIR/test_codon.py" >"$tmp_out" 2>&1
  exit_code=$?
  set -e
  if [[ $exit_code -eq 0 ]]; then
    CO_STATUS="PASS"
    # print any normal output if present
    [[ -s "$tmp_out" ]] && cat "$tmp_out"
    echo "✅ Codon tests completed"
  else
    CO_STATUS="FAIL"
    if [[ "${SHOW_CODON_TRACE:-0}" == "1" ]]; then
      # Show the full captured output only on demand
      cat "$tmp_out"
    else
      echo "🛠️  codon still under construction 🚧 (tests failed)"
      echo "ℹ️  To see Codon’s full output, run: SHOW_CODON_TRACE=1 ./evaluate.sh"
    fi
  fi
  rm -f "$tmp_out"
else
  echo "⚠️  Codon not found; skipping Codon tests."
  CO_STATUS="SKIP"
fi

echo ""
echo "🧹 Cleanup"
rm -f ./*_motif_map.txt || true

echo ""
echo "📊 Summary"
echo "  Python: $PY_STATUS"
echo "  Codon : $CO_STATUS"

echo ""
echo "All done!"
