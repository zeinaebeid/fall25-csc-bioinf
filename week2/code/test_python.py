#!/usr/bin/env python3
# Python-only test runner (concise output). Codon will not run this file.

import os, sys
sys.path.append(os.path.dirname(__file__))

from decomposer import Decomposer
from motif_aligner import MotifAligner
from main import TandemRepeatVizWorker

def ok(expr, msg="assertion failed"):
    if not expr:
        raise AssertionError(msg)

def test_decompose_and_refine():
    d = Decomposer("DP")
    seq = "AAAAACAAAAAAAAAAATAAAAAATTAAAA"
    got = d.decompose(seq, ["AAAAAA"])
    ok(got == ["AAAAAC","AAAAAA","AAAAAT","AAAAAA","TTAAAA"])
    refined = d.refine([got])
    ok(refined and refined[0] and isinstance(refined[0], list))

def test_pipeline_small():
    w = TandemRepeatVizWorker()
    ids, aln, dec = w.generate_trplot(
        "demo",
        ["s1","s2"],
        ["ACTACT","ACT"],   # 2 sequences
        ["ACT"],
        skip_alignment=True,
        verbose=False,
    )
    ok(ids == ["s1","s2"])
    ok(len(aln) == 2 and all(isinstance(x,str) for x in aln))
    ok(len(dec) == 2 and all(isinstance(x,list) for x in dec))

def run():
    test_decompose_and_refine()
    test_pipeline_small()

if __name__ == "__main__":
    # Enable short tracing if desired
    os.environ.setdefault("TRVIZ_VERBOSE", "0")
    run()
