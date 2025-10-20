from align_codon import (
    read_fasta_one,
    global_align_score,
    local_align_score,
    fitting_align_score,
    global_affine_align_score,
)
import time

def p(fname):
    return f"../data/{fname}"

def run_case(method, lang, func, a_path, b_path, label):
    s = read_fasta_one(a_path)
    t = read_fasta_one(b_path)
    t0 = time.time()
    _ = func(s, t)
    ms = int((time.time() - t0) * 1000.0)
    print(f"{method}-{label:<11} {lang:<10} {ms}ms")

def main():
    # MT_human vs MT_orang in Codon
    run_case("global",  "codon", global_align_score,         p("MT-human.fa"), p("MT-orang.fa"), "mt_human")
    run_case("local",   "codon", local_align_score,          p("MT-human.fa"), p("MT-orang.fa"), "mt_human")
    run_case("fitting", "codon", fitting_align_score,        p("MT-human.fa"), p("MT-orang.fa"), "mt_human")
    run_case("affine",  "codon", global_affine_align_score,  p("MT-human.fa"), p("MT-orang.fa"), "mt_human")

    # q1..q5
    for k in range(1, 6):
        q = f"q{k}.fa"
        t = f"t{k}.fa"
        tag = f"q{k}"
        run_case("global",  "codon", global_align_score,        p(q), p(t), tag)
        run_case("local",   "codon", local_align_score,         p(q), p(t), tag)
        run_case("fitting", "codon", fitting_align_score,       p(q), p(t), tag)
        run_case("affine",  "codon", global_affine_align_score, p(q), p(t), tag)

if __name__ == "__main__":
    main()
