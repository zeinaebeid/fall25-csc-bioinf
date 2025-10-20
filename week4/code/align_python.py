# week4/code/align_python.py
# Linear-gap scoring
MATCH = 3
MISMATCH = -3
GAP = -2

# Affine-gap scoring
GAP_OPEN = -5
GAP_EXT = -1

NEG_INF = -10**12  # big negative int (Codon-safe idea carried over)

def _score(a, b):
    return MATCH if a == b else MISMATCH

def read_fasta_one(path):
    """Read a single-sequence FASTA and return the sequence (uppercase, no gaps)."""
    seq_lines = []
    with open(path, "r") as fh:
        for line in fh:
            if not line:
                continue
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                continue
            seq_lines.append(line.replace(" ", "").replace("-", ""))
    return "".join(seq_lines).upper()

# ---------- Global (Needleman–Wunsch), linear gap ----------
def global_align_score(s, t):
    m = len(s)
    n = len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = i * GAP
    for j in range(1, n + 1):
        dp[0][j] = j * GAP
    for i in range(1, m + 1):
        si = s[i - 1]
        for j in range(1, n + 1):
            tj = t[j - 1]
            diag = dp[i - 1][j - 1] + _score(si, tj)
            up = dp[i - 1][j] + GAP
            left = dp[i][j - 1] + GAP
            dp[i][j] = max(diag, up, left)
    return dp[m][n]

# ---------- Local (Smith–Waterman), linear gap ----------
def local_align_score(s, t):
    m = len(s)
    n = len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    best = 0
    for i in range(1, m + 1):
        si = s[i - 1]
        for j in range(1, n + 1):
            tj = t[j - 1]
            diag = dp[i - 1][j - 1] + _score(si, tj)
            up = dp[i - 1][j] + GAP
            left = dp[i][j - 1] + GAP
            val = max(0, diag, up, left)
            dp[i][j] = val
            if val > best:
                best = val
    return best

# ---------- Semi-global (fitting): align full query to a substring of target ----------
def fitting_align_score(query, target):
    s = query
    t = target
    m = len(s)
    n = len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + GAP      # penalize leading gaps in query
    # first row stays 0 → free target prefix
    for i in range(1, m + 1):
        si = s[i - 1]
        for j in range(1, n + 1):
            tj = t[j - 1]
            diag = dp[i - 1][j - 1] + _score(si, tj)
            up = dp[i - 1][j] + GAP
            left = dp[i][j - 1] + GAP
            dp[i][j] = max(diag, up, left)
    # best over last row (free target suffix)
    best = dp[m][0]
    for j in range(1, n + 1):
        if dp[m][j] > best:
            best = dp[m][j]
    return best

# ---------- Global with affine gaps (Gotoh) ----------
def global_affine_align_score(s, t):
    m = len(s)
    n = len(t)
    M = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    X = [[NEG_INF] * (n + 1) for _ in range(m + 1)]  # gap in t (vertical)
    Y = [[NEG_INF] * (n + 1) for _ in range(m + 1)]  # gap in s (horizontal)

    M[0][0] = 0
    for i in range(1, m + 1):
        X[i][0] = GAP_OPEN + (i - 1) * GAP_EXT
    for j in range(1, n + 1):
        Y[0][j] = GAP_OPEN + (j - 1) * GAP_EXT

    for i in range(1, m + 1):
        si = s[i - 1]
        for j in range(1, n + 1):
            tj = t[j - 1]
            X[i][j] = max(M[i - 1][j] + GAP_OPEN, X[i - 1][j] + GAP_EXT)
            Y[i][j] = max(M[i][j - 1] + GAP_OPEN, Y[i][j - 1] + GAP_EXT)
            M[i][j] = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + _score(si, tj)

    return max(M[m][n], X[m][n], Y[m][n])
