# AI Usage Report (Week 1)

This document discloses where and how I used AI/translation tools during the Week 1 assignment.

## Tools & Versions

- **LLM:** ChatGPT — *GPT-5 Thinking*  
  - Accessed on: 2025-09-12 - 2025-09-17
- **Translation:** Google Translate  
  - Used to translate README/notes from Mandarin → English.
- **Reference Docs:** Exaloop Codon documentation (docs.exaloop.io)

---

## What AI Helped With

1. **Understanding the original repo**
   - Summarized the Mandarin README and clarified input files, expected outputs, and run instructions (cross-checked with Google Translate).

2. **Codon conversion support**
   - Identified Codon stdlib differences (avoid `pathlib`, limited `os.path`).
   - Suggested `join()` helper for paths.
   - Proposed iterative DFS to replace recursive depth computation in Codon, avoiding recursion-depth and stack issues.

3. **Code artifacts produced with AI assistance**
   - `week1/code/main_codon.py` (optional `[n_contigs]` arg)
   - `week1/code/utils_codon.py` (FASTA reader + `join`)
   - `week1/code/dbg_codon.py` (iterative, Codon-friendly DBG; frees child cache on finalize, uses `difference_update`)
   - `week1/evaluate.sh` (runs Python & Codon, times runs, computes N50, prints table)

4. **Debugging & tips**
   - Commands to run Python and Codon.
   - Explanations for Codon “Killed: 9” (memory pressure) and how to work around the issue (limit contigs / run on higher-RAM server).


## Examples of Prompts (summarized)

> *“Read docs.exaloop.io to understand how Codon works. Help me with converting `main.py`, `utils.py`, `dbg.py` to Codon-compatible code. Keep outputs identical to Python.”*

> *“My Codon run errors on `os.path.join`/`exists`. Give me a Codon-compatible version and explain whats going on here.”*

> *“The recursive `_get_depth()` doesn't really work. Give me an iterative DFS version that keeps the same ordering and results.”*

> *“Write an `evaluate.sh` that runs Python and Codon across data1..data4, times runs, computes N50 from `contig.fasta`, and prints a table with the following structure (attached a photo of table layout).”*

> *“Why does Codon gets ‘Killed: 9’ on data2/data4 but Python finishes?”*

> *"How do I change the size of text in Markdown?"*


---

## Human-Touch

- I **reviewed and tested** all the suggested code locally:
  - Verified Codon outputs match Python (N50/contigs) on `data1` and `data3`.
  - Noticed Codon issues on `data2` and `data4` on an 8 GB laptop; documented mitigation and server recommendation.
- I adjusted CI details to match my repository layout as outlined in the deliverable 1 requirements on Piazza.

---

## Limitations

- AI suggestions were tailored to Codon 0.19.x and Python 3.12 on my local macOS; other environments might not get the same results with the same prompts as above.
- I used an LLM (ChatGPT — GPT-5 Thinking) that responds based on the **ongoing conversation** and the **information I provided during this session**. Because LLM outputs depend on session history, model version, and prompt phrasing, **the same prompts may not produce identical answers** in a different session, on another device, or for another user. 
- Translation (Google Translate) was used for convenience - simply uploaded the report document to which it output a translated version. 

---

## Statement

I used **ChatGPT (GPT-5 Thinking)** to help translate, plan, and draft Codon-compatible code. I used **Google Translate** to translate the original README from Mandarin. All generated code and text were reviewed and tested by me before including them in the repo.
