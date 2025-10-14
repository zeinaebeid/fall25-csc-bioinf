#!/bin/bash
set -euo pipefail

export PATH="${HOME}/.codon/bin:$PATH"
export LC_ALL=C
export LANG=C

SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." >/dev/null 2>&1 && pwd)"

echo "Language    Runtime"
echo "-------------------"

# Run Python test in module mode from repo root
pushd "$REPO_ROOT" >/dev/null
python_line="$(python3 -m week3.code.test_phylo | grep 'python:')"
popd >/dev/null

# Run Codon test by script path from code dir 
pushd "${REPO_ROOT}/week3/code" >/dev/null
codon_line="$(codon run test_phylo_codon.py | grep 'codon:')"
popd >/dev/null

echo "${python_line//:/      }"
echo "${codon_line//:/      }"