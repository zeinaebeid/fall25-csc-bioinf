#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/code"

[[ -f "test.python"     ]] || { echo "Missing week4/code/test.python"; exit 1; }
[[ -f "test_codon.py"   ]] || { echo "Missing week4/code/test_codon.py"; exit 1; }

PY_OUT="$(python3 test.python)"
CODON_OUT="$(codon run -release test_codon.py)"

printf "Method            Language    Runtime\n"
printf "--------------------------------------\n"
printf "%s\n" "${PY_OUT}"
printf "%s\n" "${CODON_OUT}"
