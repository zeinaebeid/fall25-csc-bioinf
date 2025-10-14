# AI Usage Report (Week 3)

This document discloses where and how I used AI/translation tools during the Week 3 assignment.

## Tools & Versions
- **LLM:** Claude — *Claude Sonnet 4.5*  
  - Accessed on: 2025-10-10 - 2025-10-14
- **Reference Docs:** [Biotite phylo module](https://github.com/biotite-dev/biotite)

## What AI Helped With

1. **Initial code understanding**
   - Analyzing Biotite's Cython source files (`tree.pyx`, `upgma.pyx`, `nj.pyx`)
   - Explaining the differences between UPGMA and Neighbor-Joining algorithms
   - Breaking down the Newick notation parser logic


2. **Codon-specific porting**
   - Implementing custom `set.__hash__()` extension for TreeNode hashing (frozenset alternative)
   - Converting `TreeNode._children` from tuple to list for dynamic sizing
   - Fixing type annotations and explicit typing for Codon's type inference
   - Handling string parsing differences in Newick parser (manual colon finding vs `split()`)
   - Debugging "no module named" errors and import path issues

3. **Script optimization**
   - Wrote `evaluate.sh` to be clean and robust and return proper format


## Example Prompts (summarized)

- "Help me flatten to `week3/code/` and `week3/data/` and update all path references."

- "Fix my imports: Python test fails with 'ModuleNotFoundError: No module named tree' when running with `-m week3.code.test_phylo`"

- "My Codon test fails with 'no module named week3.code.tree_codon'. I'm running from the `week3/code/` directory."

- "Can this bash script be written better?" [provided evaluate.sh]


## Limitations

- AI suggestions were tailored to Codon 0.17.0+ and Python 3.12 on macOS; other environments might require adjustments.

- I used an LLM (Claude Sonnet 4.5) that responds based on the **ongoing conversation** and the **information I provided during this session**. Because LLM outputs depend on session history, model version, and prompt phrasing, **the same prompts may not produce identical answers** in a different session, on another device, or for another user.

- The AI helped with code structure and debugging but did not write the core algorithms (UPGMA, Neighbor-Joining) from scratch—those were ported from Biotite's implementation with AI assistance in understanding the Cython source.

- Some Codon-specific quirks (like the set hashing extension) were provided in course hints; AI helped implement and integrate them but didn't discover them independently.

