# week4/code/align_codon.codon
# Codon implementations; same logic & signatures as Python ones.

MATCH = 3
MISMATCH = -3
GAP = -2

GAP_OPEN = -5
GAP_EXT = -1

NEG_INF = -10**12

def _score(a, b):
    return MATCH if a == b else MISMATCH

def read_fasta_one(path):
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

def global_align_score(s, t):
    m = len(s)
    n = len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    i = 1
    while i <= m:
        dp[i][0] = i * GAP
        i += 1
    j = 1
    while j <= n:
        dp[0][j] = j * GAP
        j += 1
    i = 1
    while i <= m:
        si = s[i - 1]
        j = 1
        while j <= n:
            tj = t[j - 1]
            diag = dp[i - 1][j - 1] + _score(si, tj)
            up = dp[i - 1][j] + GAP
            left = dp[i][j - 1] + GAP
            best = diag
            if up > best:
                best = up
            if left > best:
                best = left
            dp[i][j] = best
            j += 1
        i += 1
    return dp[m][n]

def local_align_score(s, t):
    m = len(s)
    n = len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    best = 0
    i = 1
    while i <= m:
        si = s[i - 1]
        j = 1
        while j <= n:
            tj = t[j - 1]
            diag = dp[i - 1][j - 1] + _score(si, tj)
            up = dp[i - 1][j] + GAP
            left = dp[i][j - 1] + GAP
            val = diag
            if up > val:
                val = up
            if left > val:
                val = left
            if val < 0:
                val = 0
            dp[i][j] = val
            if val > best:
                best = val
            j += 1
        i += 1
    return best

def fitting_align_score(query, target):
    s = query
    t = target
    m = len(s)
    n = len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    i = 1
    while i <= m:
        dp[i][0] = dp[i - 1][0] + GAP
        i += 1
    i = 1
    while i <= m:
        si = s[i - 1]
        j = 1
        while j <= n:
            tj = t[j - 1]
            diag = dp[i - 1][j - 1] + _score(si, tj)
            up = dp[i - 1][j] + GAP
            left = dp[i][j - 1] + GAP
            best = diag
            if up > best:
                best = up
            if left > best:
                best = left
            dp[i][j] = best
            j += 1
        i += 1
    best = dp[m][0]
    j = 1
    while j <= n:
        if dp[m][j] > best:
            best = dp[m][j]
        j += 1
    return best

def global_affine_align_score(s, t):
    m = len(s)
    n = len(t)
    M = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    X = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    Y = [[NEG_INF] * (n + 1) for _ in range(m + 1)]

    M[0][0] = 0
    i = 1
    while i <= m:
        X[i][0] = GAP_OPEN + (i - 1) * GAP_EXT
        i += 1
    j = 1
    while j <= n:
        Y[0][j] = GAP_OPEN + (j - 1) * GAP_EXT
        j += 1

    i = 1
    while i <= m:
        si = s[i - 1]
        j = 1
        while j <= n:
            tj = t[j - 1]
            xo = M[i - 1][j] + GAP_OPEN
            xe = X[i - 1][j] + GAP_EXT
            Xij = xo if xo >= xe else xe
            yo = M[i][j - 1] + GAP_OPEN
            ye = Y[i][j - 1] + GAP_EXT
            Yij = yo if yo >= ye else ye
            dm = M[i - 1][j - 1]
            dx = X[i - 1][j - 1]
            dy = Y[i - 1][j - 1]
            prev = dm
            if dx > prev:
                prev = dx
            if dy > prev:
                prev = dy
            M[i][j] = prev + _score(si, tj)
            X[i][j] = Xij
            Y[i][j] = Yij
            j += 1
        i += 1

    end = M[m][n]
    if X[m][n] > end:
        end = X[m][n]
    if Y[m][n] > end:
        end = Y[m][n]
    return end
