# AI Usage Report (Week 2)

This document discloses where and how I used AI/translation tools during the Week 2 assignment.

## Tools & Versions

- **LLM:** ChatGPT — *GPT-5 Thinking*  
  - Accessed on: 2025-09-28 - 2025-10-01
- **Reference Docs:** [TRviz Repo](https://github.com/Jong-hun-Park/trviz) and [Paper](https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad058/7143378)

## What AI Helped With

1. Repo scaffolding & run flow
2. Wrote evaluate.sh that:
   - Exports PYTHONPATH
   - Runs Python tests first, then Codon tests
   - Suppresses Codon backtraces 
3. Refined Decomposer DP routine to mirror the TRViz logic (simplified, no HMM/Cython)
4. Added TRVIZ_VERBOSE env toggle to control debug prints.
5. Codon-compatibility passes (give or take since codon failed)
6. Debugging & Commands
7. Tips on handling PYTHONPATH, relative imports, and suppressing Codon tracebacks


## Example Prompts (summarized)

“Make a Codon-friendly decomposer—no typing.TYPE_CHECKING, no dynamic kwargs iteration; expose explicit params and keep behavior identical to Python where possible.”

“Write evaluate.sh to run Python then Codon, hide Codon backtrace by default, and print a friendly ‘under construction’ line if Codon fails.”

“My imports break with attempted relative import with no known parent package. Adjust paths and give me exact commands.”

“Reduce verbosity—use an env flag so tests don’t flood the terminal.”


## Limitations

- AI suggestions were tailored to Codon 0.19.x and Python 3.12 on my local macOS; other environments might not get the same results with the same prompts as above.
- I used an LLM (ChatGPT — GPT-5 Thinking) that responds based on the **ongoing conversation** and the **information I provided during this session**. Because LLM outputs depend on session history, model version, and prompt phrasing, **the same prompts may not produce identical answers** in a different session, on another device, or for another user. 

