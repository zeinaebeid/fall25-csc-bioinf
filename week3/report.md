# Week 3 — Phylogenetic Trees (UPGMA & Neighbor-Joining) in Python & Codon

**Repo:** `zeinaebeid/fall25-csc-bioinf`

## Environment
- **OS (local):** macOS
- **Python:** `python3 --version` → 3.12.x
- **Codon:** `codon --version` → 0.17.0+
- **Shell:** `zsh`
- **Hardware (local):** 8 GB RAM (laptop)

## What's included
- `tree.py` / `tree_codon.py` — Tree and TreeNode classes (ported from Biotite Cython)
- `upgma.py` / `upgma_codon.py` — UPGMA clustering algorithm
- `nj.py` / `nj_codon.py` — Neighbor-Joining algorithm
- `test_phylo.py` — Python tests (3 tests)
- `test_phylo_codon.py` — Codon tests (3 tests)
- `evaluate.sh` — runs both test, reports timing, and formats as requested
- `__init__.py` — module initialization

## Goals
1. **Port Biotite's phylo module** from Cython to pure Python and Codon
2. **Implement two clustering algorithms**: UPGMA and Neighbor-Joining
3. **Test coverage**: 3 tests each for Python and Codon
4. **Performance comparison**: measure and report execution time in milliseconds


## Results
**Python**: ✅ PASS (all 3 tests)
**Codon**: ✅ PASS (all 3 tests)

Timing comparison:
```
Language    Runtime
-------------------
python      5ms
codon       2ms
```


## Gotchas / Lessons Learned

1. **Module imports are tricky**: Python's `-m` module mode requires relative imports, while Codon needs simple imports. This required maintaining separate import styles for `test_phylo.py` vs `test_phylo_codon.py`.

2. **TreeNode hashings**: Biotite uses `frozenset` for hashing children, which Codon doesn't support natively. Solution: extended the `set` class with a custom `__hash__()` implementation (as suggested in Piazza).

3. **Reading Cython code**: The original Biotite source is in Cython, which mixes Python and C. Understanding the low-level optimizations helped inform the pure Python port.

4. **Testing infrastructure**: Having separate Python and Codon test files made debugging much easier than trying to make one file work for both.

## Time Investment
Approximately **11 hours** total:
- 3 hours: understanding Biotite's Cython source
- 4 hours: initial Python implementation
- 3 hours: Codon porting and debugging type issues (with the help of AI)
- 1 hour: documentation
