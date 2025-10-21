# Week 4 — Sequence Alignment (Codon & Python)

**Repo:** `zeinaebeid/fall25-csc-bioinf`

## Environment
- **OS (local):** macOS
- **Python:** `python3 --version` → 3.12.x
- **Codon:** `codon --version` → 0.17.0+
- **Shell:** `zsh`
- **Hardware (local):** 8 GB RAM (laptop)


## Goal: 

**Implement four sequence alignment algorithms in Python/Codon:**

1. **Global Alignment** (Needleman-Wunsch)
2. **Local Alignment** (Smith-Waterman)
3. **Semi-global Alignment** (Fitting alignment)
4. **Affine Gap Penalty Global Alignment**

> Implementation constants:
> ```python
> match_score = 3
> mismatch_score = -3
> gap_penalty = -2
> gap_open_penalty = -5
> gap_extension_penalty = -1
> ```


## How to run

From the repo root:

`chmod +x week4/evaluate.sh`

`./week4/evaluate.sh`

## How it works

- Reads two DNA sequences from FASTA files, then runs one of the alignment algorithms to compare them and produce a score.
- Global, local, and fitting use a simple DP matrix with a fixed gap cost.
- Affine global uses three matrices to handle gap-open and gap-extend separately.
- Each run and print the method, language, and runtime as below.

FASTA test files were adapted from: https://github.com/lh3/ksw2/tree/master/test

## Output
```
Method           Language     Runtime
--------------------------------------
global-mt_human   python      75254.518ms
global-q1         python      0.010ms
global-q2         python      0.047ms
.
.
.
global-mt_human   codon       14064.074ms
global-q1         codon       0.004ms
......
```

## Gotchas

- **Local vs CI behavior.** All tests passed on my laptop (see screenshot below), but the GitHub CI job consistently canceled right after the Python `affine-q5` row when it tried to run the **MT-human vs MT-orang** case in Python. This is likely due to runtime/resource differences on the CI runner.

- **Fix applied for CI.** I updated `evaluate.sh` to **skip MT-human in Python** and run the MT pair **only under Codon**, which is faster and stable on CI. The q1–q5 sets still run under both Python and Codon, preserving parity checks.



