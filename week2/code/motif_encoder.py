from typing import Dict, List, Optional, Tuple
import string

# Alphabet (letters+digits+safe ASCII), reserve '?' for private motifs
PRIVATE = "?"
SKIP = set("()=<>?-")
ALPHABET = list(string.ascii_lowercase + string.ascii_uppercase + string.digits)
ALPHABET.extend([chr(c) for c in range(33, 127) if chr(c) not in SKIP and chr(c) not in ALPHABET])

def _count_motifs(vntrs: List[List[str]]) -> Dict[str, int]:
    c: Dict[str, int] = {}
    for row in vntrs:
        for m in row:
            c[m] = c.get(m, 0) + 1
    return c

class MotifEncoder:
    """
    Assigns symbols to motifs by frequency; groups rare motifs into '?' if needed.
    """

    def __init__(self) -> None:
        self.motif_to_symbol: Dict[str, str] = {}
        self.symbol_to_motif: Dict[str, str] = {}
        self.motif_counter: Dict[str, int] = {}
        self.score_matrix: Dict[str, Dict[str, float]] = {}  # reserved

    @staticmethod
    def _threshold(vntrs: List[List[str]], label_count: Optional[int]) -> int:
        max_labels = len(ALPHABET) - 1  # reserve one for '?'
        if label_count is not None:
            max_labels = max(0, label_count - 1)
        counts = sorted(_count_motifs(vntrs).items(), key=lambda kv: (-kv[1], kv[0]))
        used = 0
        thr = 0
        for _, cnt in counts:
            used += 1
            if used > max_labels:
                thr = cnt
                break
        return thr

    @staticmethod
    def _write_map(path: str, m2s: Dict[str, str], cnt: Dict[str, int]) -> None:
        items = sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))
        with open(path, "w") as f:
            for motif, freq in items:
                f.write(f"{motif}\t{m2s[motif]}\t{freq}\n")

    def encode(
        self,
        decomposed_vntrs: List[List[str]],
        motif_map_file: str,
        label_count: Optional[int] = None,
        auto: bool = True
    ) -> List[str]:
        cnt = _count_motifs(decomposed_vntrs)
        self.motif_counter = dict(cnt)
        thr = self._threshold(decomposed_vntrs, label_count) if (auto or label_count) else 0

        # split into normal/private
        normal = []
        private = []
        for m, c in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0])):
            (private if c <= thr else normal).append(m)

        m2s: Dict[str, str] = {}
        s2m: Dict[str, str] = {}

        # assign normal first
        if len(normal) + (1 if private else 0) > len(ALPHABET):
            raise ValueError("Too many unique motifs for available alphabet.")
        for i, m in enumerate(normal):
            sym = ALPHABET[i]
            m2s[m] = sym
            s2m[sym] = m

        # group private
        if private:
            for m in private:
                m2s[m] = PRIVATE
            # represent PRIVATE by last private motif (arbitrary but deterministic)
            s2m[PRIVATE] = private[-1]

        self.motif_to_symbol = dict(m2s)
        self.symbol_to_motif = dict(s2m)
        self.score_matrix = {}  # not used in this simplified pipeline

        # persist map
        self._write_map(motif_map_file, m2s, cnt)

        # encode each VNTR row to symbol string
        encoded: List[str] = []
        for row in decomposed_vntrs:
            encoded.append("".join(m2s[m] for m in row))
        return encoded
