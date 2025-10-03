## Week 2 — TRViz (DP) · Python + Codon

**Repo:** `<zeinaebeid/fall25-csc-bioinf>`  

## Environment

- **OS (local):** macOS 
- **Python:** `python3 --version` → 3.12.x  
- **Codon:** `codon --version` → 0.19.x  
- **Shell:** `zsh`  
- **Hardware (local):** 8 GB RAM (laptop)

## What’s included

`decomposer.py` — pure-Python DP decomposer + refine (TRViz-like, no NumPy).

`motif_encoder.py` — maps motifs → symbols (handles “private” ?).

`motif_aligner.py` — simple pad-to-length “center-star”.

`main.py` — pipeline: decompose → refine → encode → align.

`test_python.py` — Python tests (pass).

`test_codon.py` — Codon-safe tests (currently fail one assert).

`evaluate.sh` — runs both, prints a clean summary.

## Goals

1. **Re-implement TRViz core (DP)** in pure Python, Codon-friendly (no NumPy/HMM/Cython).
2. **Modular pipeline**: Decomposer → refine → MotifEncoder → MotifAligner.
3. **Testing**: separate test_python.py and test_codon.py.
4. **Automation**: evaluate.sh runs both and prints a clean summary.

## How to Run

`cd week2`

`chmod +x evaluate.sh`

`export PYTHONPATH=week2/code`

### Run everything
`./evaluate.sh`

### Or separately:
`python3 code/test_python.py`

`codon run code/test_codon.py`   # may fail; backtrace hidden by evaluate.sh

## Results

**Python**: ✅ PASS

**Codon**: ❌ FAIL (one assertion in simplified pipeline test).

- DP core, encoder, and aligner execute under Codon; one expectation differs from Python (tie-break/alignment padding nuance).

- Backtrace is hidden by default in evaluate.sh. Set SHOW_CODON_TRACE=1 to view.

## Gotchas / Lessons Learned

- This was genuinely hard for me to conceptualize and implement end-to-end. After multiple attempts, I accepted that my Codon path wasn’t working on my machine. I prioritized a clean, fully passing Python pipeline and wrapped Codon in `evaluate.sh` with graceful failure.
