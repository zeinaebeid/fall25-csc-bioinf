import time
import sys
from typing import List


match_score = 3
mismatch_score = -3
gap_penalty = -2
gap_open_penalty = -5
gap_extension_penalty = -1
NEG_INF = -1_000_000


def read_fasta(filename: str) -> str:
    """Read FASTA file and return sequence"""
    seq_lines: List[str] = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
    return "".join(seq_lines)


def global_align(seq1: str, seq2: str) -> int:
    """Needleman-Wunsch global alignment"""
    n = len(seq1)
    m = len(seq2)
    
    score: List[List[int]] = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    
    for i in range(1, n + 1):
        score[i][0] = i * gap_penalty
    for j in range(1, m + 1):
        score[0][j] = j * gap_penalty
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diag = score[i - 1][j - 1] + match_score
            else:
                diag = score[i - 1][j - 1] + mismatch_score
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty
            
            best = diag
            if up > best:
                best = up
            if left > best:
                best = left
            score[i][j] = best
    
    return score[n][m]


def local_align(seq1: str, seq2: str) -> int:
    """Smith-Waterman local alignment"""
    n = len(seq1)
    m = len(seq2)
    
    score: List[List[int]] = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    max_score = 0
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diag = score[i - 1][j - 1] + match_score
            else:
                diag = score[i - 1][j - 1] + mismatch_score
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty
            
            best = diag
            if up > best:
                best = up
            if left > best:
                best = left
            if 0 > best:
                best = 0
            score[i][j] = best
            
            if best > max_score:
                max_score = best
    
    return max_score


def fitting_align(seq1: str, seq2: str) -> int:
    """Semi-global (fitting) alignment"""
    n = len(seq1)
    m = len(seq2)
    
    score: List[List[int]] = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    
    for j in range(1, m + 1):
        score[0][j] = j * gap_penalty
    
    for i in range(1, n + 1):
        score[i][0] = 0
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diag = score[i - 1][j - 1] + match_score
            else:
                diag = score[i - 1][j - 1] + mismatch_score
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty
            
            best = diag
            if up > best:
                best = up
            if left > best:
                best = left
            score[i][j] = best
    
    best_score = score[0][m]
    for i in range(1, n + 1):
        if score[i][m] > best_score:
            best_score = score[i][m]
    
    return best_score


def affine_align(seq1: str, seq2: str) -> int:
    """Global alignment with affine gap penalties"""
    n = len(seq1)
    m = len(seq2)
    
    lower: List[List[int]] = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    middle: List[List[int]] = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    upper: List[List[int]] = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    
    for i in range(1, n + 1):
        lower[i][0] = gap_open_penalty + (i - 1) * gap_extension_penalty
        middle[i][0] = gap_open_penalty + (i - 1) * gap_extension_penalty
        upper[i][0] = NEG_INF
    
    for j in range(1, m + 1):
        upper[0][j] = gap_open_penalty + (j - 1) * gap_extension_penalty
        middle[0][j] = gap_open_penalty + (j - 1) * gap_extension_penalty
        lower[0][j] = NEG_INF
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            lower[i][j] = max(
                lower[i - 1][j] + gap_extension_penalty,
                middle[i - 1][j] + gap_open_penalty
            )
            
            upper[i][j] = max(
                upper[i][j - 1] + gap_extension_penalty,
                middle[i][j - 1] + gap_open_penalty
            )
            
            if seq1[i - 1] == seq2[j - 1]:
                score = match_score
            else:
                score = mismatch_score
            
            middle[i][j] = max(
                lower[i][j],
                middle[i - 1][j - 1] + score,
                upper[i][j]
            )
    
    return middle[n][m]


def main():
    if len(sys.argv) != 5:
        print("Usage: python align.py <method> <lang> <file1> <file2>")
        sys.exit(1)
    
    method = sys.argv[1]
    lang = sys.argv[2]
    file1 = sys.argv[3]
    file2 = sys.argv[4]
    
    seq1 = read_fasta(file1)
    seq2 = read_fasta(file2)
    
    start = time.time()
    
    if method.startswith("global"):
        score = global_align(seq1, seq2)
    elif method.startswith("local"):
        score = local_align(seq1, seq2)
    elif method.startswith("fitting"):
        score = fitting_align(seq1, seq2)
    elif method.startswith("affine"):
        score = affine_align(seq1, seq2)
    else:
        print(f"Unknown method: {method}")
        sys.exit(1)
    
    end = time.time()
    runtime_ms = (end - start) * 1000
    
    # Format output
    method_col = method.ljust(18)
    lang_col = lang.ljust(12)
    runtime_col = f"{runtime_ms:.3f}ms"
    
    print(f"{method_col}{lang_col}{runtime_col}")


if __name__ == "__main__":
    main()
