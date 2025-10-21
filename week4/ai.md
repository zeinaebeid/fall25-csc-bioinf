# AI Usage Report (Week 4)

This document discloses where and how I used AI/translation tools during the Week 4 assignment.

## Tools & Versions
- **LLM:** ChatGPT — *GPT-5 Thinking*  
  - Accessed on: 2025-10-20 – 2025-10-21
- **Reference Data:** FASTA tests adapted from [ksw2 `test/`](https://github.com/lh3/ksw2/tree/master/test)
- **Prompts Provided by Course:** I **also used the instructor’s “Sample AI prompts for converting Python to Codon”** as guidance for phrasing my requests and for checking Codon-specific constraints.

## What AI Helped With

1. **Algorithm support (Python ↔ Codon)**
   - Verified/clarified recurrences for:
     - Global (Needleman–Wunsch), Local (Smith–Waterman)
     - Semi-global (fitting)
     - Affine Global (Gotoh: `lower/middle/upper` matrices)

2. **Codon compatibility fixes**
   - Avoided `Any`, mixed-type lists, and empty-collection ambiguity.
   - Replaced `%`/`.format()` with f-strings (Codon-safe).
   - Removed `sys.implementation`/`builtins` checks; used a Codon-safe heuristic instead.
   - Kept `NEG_INF` as an **int** for affine initialization.

3. **I/O normalization & correctness**
   - Recommended FASTA normalization (uppercase, strip `-`/spaces) to prevent score mismatches.

4. **Test support**
   - Created `evaluate.sh` to:
     - Run **q1–q5** under both Python and Codon,
     - **Skip MT-human** under Python (CI cancellation),
     - Run MT pair under Codon only, with `--` to pass argv safely.


## Example Prompts (summarized)

- "Create an `alignment.py` compatible with Codon and Python that supports all 4 required methods.”
- “GitHub CI cancels on Python MT—modify `evaluate.sh` to skip Python MT and run MT in Codon.”
- **From the professor’s sample prompts:** “Convert this Python DP to Codon, ensuring homogeneous list types and no empty-literal ambiguity,” “Replace unsupported formatting/imports,” “Mark Optionals and class fields explicitly.”

## Limitations

- **Environment differences:** Locally everything passed; CI required skipping Python MT and running MT only under Codon.
- **Session variance:** LLM output depends on prompts/model/session; identical results aren’t guaranteed across sessions.


## Academic Integrity 

- I understand the algorithms (Needleman–Wunsch, Smith–Waterman, fitting, Gotoh) and can explain recurrences, boundary conditions, and Codon-specific constraints.
- AI served as an assistant for Codon-porting guidelines, debugging, scripting, and documentation—not as a substitute for my understanding or implementation.
