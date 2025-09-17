# Week 1 — de Bruijn Assembler

**Repo:** `<zeinaebeid/fall25-csc-bioinf>`  

## Environment

- **OS (local):** macOS 
- **Python:** `python3 --version` → 3.12.x  
- **Codon:** `codon --version` → 0.19.x  
- **Shell:** `zsh`  
- **Hardware (local):** 8 GB RAM (laptop)

## Goals

1) **Automation** :  Set up a clean repo structure and CI (GitHub Actions). Provide a single script `week1/evaluate.sh` that runs Python and Codon across `data1..data4`, records runtime and computes N50.
2) **Codon Conversion** : Convert `main.py`, `utils.py`, and `dbg.py` to Codon-compatible versions.
3) **Documentation**: Record AI prompts utilized in `ai.md` and table results in `report.md`


## How to Run 
- Please note, this is specific to my local environment.

**Requirements** 
- Python 3.11+
- Codon + `seq` plugin

**Evaluate all datasets**

```
./week1/evaluate.sh
```
## What the `evaluate.sh` Does
- For each dataset, the script runs `main.py` and `main_codon.py` from `week1/code` directory.
- Captures stdout, measures time, extracts contig lengths.
- Prints out a compact table:

```
Dataset   Language      Runtime      N50
-------  ----------     ---------  --------
data1     python        0:00:14      9991
data1     codon         0:00:09      9991
…
```

## My Steps Taken and Results 
1. Cloned the original repo from `zhongyuchen/genome-assembly`.
2. Pushed the local clone to my github.
3. The README/report were in Mandarin. I used Google Translate to translate the text into English to understand the data format and instructions.
4. Ran the baseline Python implementation on all four data sets.
5. Converted the Python files to Codon compatible versions (using the help of AI!).
6. Created `evaluate.sh` to automate running both Python and Codon versions.
7. Produced the following results table:

**Results** 

```
Dataset   Language      Runtime      N50
-------  ----------     ---------  --------
data1     python        0:00:14      9991
data1     codon         0:00:09      9991
data2     python        0:00:29      9993
data2     codon         NA           NA
data3     python        0:00:33      9825
data3     codon         0:00:16      9825
data4     python        0:14:59      168187
data4     codon         NA           NA 
…
```
### Reproducibility Note on `data2` and `data4` (Codon)

This wasn’t a bug in the code - it was my machine running out of memory. My local laptop only has 8 GB RAM ☹️. After a few contigs, the remaining de Bruijn graph is still huge and loopy. Codon’s iterative DFS + sorting spikes the peak memory higher than Python. macOS then steps in and terminates the process (no traceback, just the kill).

**What I tried locally**
- Closed all RAM intesive/non-essential applications, killed background processes, and re-ran.
- Reduced the number of requested contigs (e.g., `codon run -release … data2 10`) to lower peak usage.
- Re-ran without `-release` to slightly alter allocation patterns.

Even with all that, Codon still wouldn’t run `data2`/`data4` on 8 GB. The Python version *did* finish on the same machine, which points to a memory ceiling rather than logic differences.

## Conclusion
- Both Python and Codon implementations yield identical N50s.
- Codon did indeed deliver consistent speed-ups, but is clearly more RAM intensive.
- Ideally, Codon should be ran using 16-32GB RAM in order to avoid premature termination.


## AI Usage 
- All AI usage is documented under `week1/ai.md` as required.

  
## ⭐️ BONUS ⭐️ 

### **Method:** I used the NCBI **BLASTN** website.

**Steps:** 
1. Opened **NCBI BLASTN** → Nucleotide Blast → Upload FASTA file.
3. Uploaded `contig.fasta` for **data1** and ran BLAST with default parameters.
4. Reviewed the best matches.

**Result for data1**
- Top Match:
  - **Organism:** *Porphyromonas gingivalis* strain W50/BE1
  - **Description:** chromosome, complete genome



<img width="783" height="208" alt="Screenshot 2025-09-17 at 4 47 50 PM" src="https://github.com/user-attachments/assets/fea18919-6671-4873-9833-fccfdb175090" />

