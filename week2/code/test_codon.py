#!/usr/bin/env python3
# Codon-safe test runner (no sys.path hacks). Run with:
#   PYTHONPATH=week2/code codon run week2/code/test_codon.py

from decomposer import Decomposer
from motif_aligner import MotifAligner
from main import TandemRepeatVizWorker

def ok(expr, msg="assertion failed"):
    if not expr:
        raise AssertionError(msg)

def test_decompose_basic():
    d = Decomposer("DP")
    parts = d.decompose("AAAATGATGATGCCC", ["ATG"])
    ok(sum(1 for x in parts if x == "ATG") == 3)

def test_refine_majority():
    d = Decomposer("DP")
    refined = d.refine([["ATG","ATG"], ["ATG"]])
    ok(len(refined[1]) == 2)

def test_pipeline_pad_alignment():
    w = TandemRepeatVizWorker()
    ids, aln, dec = w.generate_trplot(
        "demo",
        ["s1","s2"],
        ["ACT","ACTACT"],
        ["ACT"],
        skip_alignment=True,
        verbose=False,
    )
    ok(ids == ["s1","s2"])
    ok(len(aln[0]) == len(aln[1]))
    ok("".join(dec[0]) == "ACT")

def run():
    test_decompose_basic()
    test_refine_majority()
    test_pipeline_pad_alignment()

if __name__ == "__main__":
    run()
