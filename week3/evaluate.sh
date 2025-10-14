#!/bin/bash
set -euo pipefail

export PATH="${HOME}/.codon/bin:${PATH}"
export LC_ALL=C
export LANG=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "Language    Runtime"
echo "-------------------"

# Run Python test in module mode from repo root
python_line=$(cd "$REPO_ROOT" && python3 -m week3.code.test_phylo | grep 'python:')

# Run Codon test from code directory
codon_line=$(cd "${REPO_ROOT}/week3/code" && codon run test_phylo_codon.py | grep 'codon:')

# Format and display results
echo "${python_line//:/ }"
echo "${codon_line//:/ }"