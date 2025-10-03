# week2/code/decomposer.py
# Codon-safe DP decomposer + column-majority refine()

from typing import List, Tuple

NEG_INF = -10**12  # internal "minus infinity"


# --- small helpers -----------------------------------------------------------
def is_valid_sequence(sequence: str) -> bool:
    """DNA alphabet check (A/C/G/T only)."""
    for ch in sequence:
        if ch not in ("A", "C", "G", "T"):
            return False
    return True


def _env_verbose() -> bool:
    """
    Codon-safe check of TRVIZ_VERBOSE in the environment.
    We do not call os.environ.get(); Codon EnvMap has no .get.
    """
    try:
        import os  # Codon supports this
        if "TRVIZ_VERBOSE" in os.environ:
            v = os.environ["TRVIZ_VERBOSE"]
            return v is not None and v != "" and v != "0"
        return False
    except Exception:
        return False


def _argmax3(a: float, b: float, c: float) -> int:
    """Return index of max among three numbers: 0/1/2."""
    if a >= b and a >= c:
        return 0
    if b >= c:
        return 1
    return 2


# --- main class --------------------------------------------------------------
class Decomposer:
    """
    Minimal, Codon-friendly decomposer:
      - mode kept for API parity, but only DP behavior is implemented here
      - refine(): column-majority padding
      - decompose(): dynamic-programming with backtracking (pure lists)
    """

    def __init__(self, mode: str = "DP"):
        if mode not in ("DP", "DP_CY", "HMM"):
            raise ValueError(str(mode) + " is invalid mode for tandem repeat decomposer.")
        self.mode = "DP"  # force DP semantics in this file

    # ------------------------------------------------------------------ refine
    @staticmethod
    def refine(decomposed_trs: List[List[str]], verbose: bool = False) -> List[List[str]]:
        """
        Column-majority refinement:
        1) Compute max row length.
        2) For each column k, choose the most frequent motif among rows that have it
           (tie-break lexicographically).
        3) Pad shorter rows to max length using that column's majority motif.
        """
        if not decomposed_trs:
            return []

        # 1) max length
        max_len = 0
        for tr in decomposed_trs:
            ln = len(tr)
            if ln > max_len:
                max_len = ln

        # 2) column majorities
        col_major: List[str] = []
        for c in range(max_len):
            counts = {}  # motif -> count
            for tr in decomposed_trs:
                if c < len(tr):
                    m = tr[c]
                    if m in counts:
                        counts[m] = counts[m] + 1
                    else:
                        counts[m] = 1
            best = None
            bestc = -1
            for k in counts:
                v = counts[k]
                if v > bestc or (v == bestc and (best is None or k < best)):
                    best = k
                    bestc = v
            if best is None:
                best = ""
            col_major.append(best)

        # 3) pad shorter rows
        out: List[List[str]] = []
        do_print = verbose or _env_verbose()
        for tr in decomposed_trs:
            new = list(tr)
            cur_len = len(new)
            if cur_len < max_len:
                for c in range(cur_len, max_len):
                    fill = col_major[c]
                    if do_print:
                        print("[refine] fill col", c, "with", fill)
                    new.append(fill)
            out.append(new)

        return out

    # --------------------------------------------------------------- public API
    def decompose(
        self,
        sequence: str,
        motifs_in: List[str],
        match_score: float = 1.0,
        mismatch_score: float = -1.0,
        insertion_score: float = -1.0,
        deletion_score: float = -1.0,
        # IMPORTANT for Codon: use a literal, not float(NEG_INF)
        min_score_threshold: float = -1e12,
        verbose: bool = False,
    ) -> List[str]:
        """
        Decompose sequence into motifs using DP (pure lists; Codon-safe).
        """
        # Basic type checks (Codon-friendly)
        if not isinstance(sequence, str):
            raise TypeError("Sequence must be a string")
        if isinstance(motifs_in, str):
            motifs = [motifs_in]
        else:
            if not isinstance(motifs_in, list):
                raise TypeError("Motifs must be a list of strings")
            motifs = motifs_in

        seq = sequence.upper()
        if not is_valid_sequence(seq):
            raise ValueError("Sequence has invalid characters: " + seq)

        motifs_up: List[str] = []
        for m in motifs:
            mu = str(m).upper()
            if not is_valid_sequence(mu):
                raise ValueError("The motif has invalid characters: " + mu)
            motifs_up.append(mu)

        return self._decompose_dp(
            seq,
            motifs_up,
            float(match_score),
            float(mismatch_score),
            float(insertion_score),
            float(deletion_score),
            float(min_score_threshold),
            bool(verbose),
        )

    # -------------------------------------------------------------------- core
    def _decompose_dp(
        self,
        sequence: str,
        motifs: List[str],
        match_score: float,
        mismatch_score: float,
        insertion_score: float,
        deletion_score: float,
        min_score_threshold: float,
        verbose: bool,
    ) -> List[str]:
        n = len(sequence)
        M = len(motifs)

        # longest motif length
        max_len = 0
        for mot in motifs:
            L = len(mot)
            if L > max_len:
                max_len = L

        # DP tables: s[i][m][j], back[i][m][j] -> (i2,m2,j2)
        s: List[List[List[float]]] = []
        back: List[List[List[Tuple[int, int, int]]]] = []
        for i in range(n + 1):
            row_s: List[List[float]] = []
            row_b: List[List[Tuple[int, int, int]]] = []
            for m_idx in range(M):
                row_s.append([0.0] * (max_len + 1))
                row_b.append([(0, m_idx, 0)] * (max_len + 1))
            s.append(row_s)
            back.append(row_b)

        # boundaries
        for m_idx in range(M):
            Lm = len(motifs[m_idx])
            for i in range(n + 1):
                for j in range(Lm + 1):
                    if i == 0 and j == 0:
                        s[0][m_idx][0] = 0.0
                        back[0][m_idx][0] = (0, m_idx, 0)
                    elif i == 0 and j != 0:
                        s[0][m_idx][j] = s[0][m_idx][j - 1] + insertion_score
                        back[0][m_idx][j] = (0, m_idx, j - 1)
                    elif i != 0 and j == 0:
                        s[i][m_idx][0] = s[i - 1][m_idx][0] + insertion_score
                        back[i][m_idx][0] = (i - 1, m_idx, 0)

        # fill
        for i in range(1, n + 1):
            si = sequence[i - 1]
            for m_idx in range(M):
                mot = motifs[m_idx]
                Lm = len(mot)
                for j in range(1, Lm + 1):
                    mj = mot[j - 1]
                    if j == 1:
                        if i == 1:
                            diag = s[i - 1][m_idx][j - 1] + (match_score if si == mj else mismatch_score)
                            left = s[i - 1][m_idx][j] + insertion_score
                            up = s[i][m_idx][j - 1] + deletion_score
                            choice = _argmax3(diag, left, up)
                            if choice == 0:
                                s[i][m_idx][j] = diag
                                back[i][m_idx][j] = (0, m_idx, 0)
                            elif choice == 1:
                                s[i][m_idx][j] = left
                                back[i][m_idx][j] = (i - 1, m_idx, j)
                            else:
                                s[i][m_idx][j] = up
                                back[i][m_idx][j] = (i, m_idx, 0)
                        else:
                            # Can start a new motif here by coming from the end of ANY motif at i-1
                            best_val = -1e30
                            best_m = 0
                            best_jend = 0
                            for t in range(M):
                                mend = s[i - 1][t][len(motifs[t])]
                                if mend > best_val:
                                    best_val = mend
                                    best_m = t
                                    best_jend = len(motifs[t])

                            motif_end = best_val + (match_score if si == mot[0] else mismatch_score)
                            left = s[i - 1][m_idx][1] + insertion_score
                            up = s[i][m_idx][0] + deletion_score
                            choice = _argmax3(motif_end, left, up)
                            if choice == 0:
                                s[i][m_idx][j] = motif_end
                                back[i][m_idx][j] = (i - 1, best_m, best_jend)
                            elif choice == 1:
                                s[i][m_idx][j] = left
                                back[i][m_idx][j] = (i - 1, m_idx, 1)
                            else:
                                s[i][m_idx][j] = up
                                back[i][m_idx][j] = (i, m_idx, 0)
                    else:
                        diag = s[i - 1][m_idx][j - 1] + (match_score if si == mj else mismatch_score)
                        left = s[i - 1][m_idx][j] + insertion_score
                        up = s[i][m_idx][j - 1] + deletion_score
                        choice = _argmax3(diag, left, up)
                        if choice == 0:
                            s[i][m_idx][j] = diag
                            back[i][m_idx][j] = (i - 1, m_idx, j - 1)
                        elif choice == 1:
                            s[i][m_idx][j] = left
                            back[i][m_idx][j] = (i - 1, m_idx, j)
                        else:
                            s[i][m_idx][j] = up
                            back[i][m_idx][j] = (i, m_idx, j - 1)

        # choose best end
        start_ptr = None  # (i, m, j)
        best_score = min_score_threshold
        for m_idx in range(M):
            Lm = len(motifs[m_idx])
            val = s[n][m_idx][Lm]
            if val > best_score:
                best_score = val
                start_ptr = (n, m_idx, Lm)

        if start_ptr is None:
            raise ValueError("No good match greater than score threshold of " + str(min_score_threshold))

        do_print = verbose or _env_verbose()
        if do_print:
            print("[decompose] best_score:", best_score)

        # backtrack: gather (i,m,j) path
        path_rev: List[Tuple[int, int, int]] = []
        cur = start_ptr
        while True:
            path_rev.append(cur)
            i0, m0, j0 = cur
            if i0 == 0 and j0 == 0:
                break
            cur = back[i0][m0][j0]

        path = list(reversed(path_rev))

        # reconstruct emitted motifs
        decomposed: List[str] = []
        cur_chars: List[str] = []
        prev_i, prev_m, prev_j = -1, -1, -1

        for (i, m, j) in path:
            if prev_m != -1 and prev_j == len(motifs[prev_m]) and j == 1 and cur_chars:
                decomposed.append("".join(cur_chars))
                cur_chars = []
            if i > prev_i and i > 0:
                cur_chars.append(sequence[i - 1])
            prev_i, prev_m, prev_j = i, m, j

        if cur_chars:
            decomposed.append("".join(cur_chars))

        if do_print:
            rec = "".join(decomposed)
            print("[decompose] motifs:", decomposed)
            print("[decompose] rec == seq ? ", rec == sequence)

        return decomposed
