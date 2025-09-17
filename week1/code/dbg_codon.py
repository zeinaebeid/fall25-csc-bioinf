# iterative, Codon-safe
from typing import Dict, List, Set

def reverse_complement(key: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    out: List[str] = []
    for ch in key[::-1]:
        out.append(comp.get(ch, ch))
    return "".join(out)

class DBG:
    # Graph data structure
    def __init__(self, k: int, data_list: List[List[str]]):
        self.k: int = k
        # node 
        self.children: Dict[int, Set[int]] = {}
        self.counts: Dict[int, int] = {}
        self.kmers: Dict[int, str] = {}
        self.depth: Dict[int, int] = {}
        self.max_child: Dict[int, int] = {}

    
        self.kmer2idx: Dict[str, int] = {}
        self.kmer_count: int = 0

        self._check(data_list)
        self._build(data_list)


    def _check(self, data_list: List[List[str]]) -> None:
        assert len(data_list) > 0
        assert len(data_list[0]) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list: List[List[str]]) -> None:
        for reads in data_list:
            for original in reads:
                rc = reverse_complement(original)
                bound = len(original) - self.k
                if bound <= 0:
                    continue
                for i in range(bound):
                    self._add_arc(original[i:i+self.k], original[i+1:i+1+self.k])
                    self._add_arc(rc[i:i+self.k],        rc[i+1:i+1+self.k])

    def _add_node(self, kmer: str) -> int:
        idx = self.kmer2idx.get(kmer, -1)
        if idx == -1:
            idx = self.kmer_count
            self.kmer2idx[kmer] = idx
            self.kmer_count += 1
            self.children[idx] = set()
            self.counts[idx] = 0
            self.kmers[idx] = kmer
            self.depth[idx] = 0
            self.max_child[idx] = -1
        self.counts[idx] = self.counts[idx] + 1
        return idx

    def _add_arc(self, kmer1: str, kmer2: str) -> None:
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.children[idx1].add(idx2)


    def _sorted_children(self, idx: int) -> List[int]:
        ch = list(self.children.get(idx, set()))
        ch.sort(key=lambda c: (-self.counts[c], c))
        return ch

    def _compute_depth_from(self, start: int,
                            status: Dict[int, int],
                            ch_cache: Dict[int, List[int]]) -> None:
   
        stack: List[int] = []
        iter_idx: Dict[int, int] = {}
        stack.append(start)

        while stack:
            v = stack[-1]
            st = status.get(v, 0)

            if st == 0:
                status[v] = 1
                # cache sorted children to avoid re-sorting
                if v not in ch_cache:
                    ch_cache[v] = self._sorted_children(v)
                iter_idx[v] = 0

            elif st == 1:
                ch = ch_cache[v]
                i = iter_idx[v]
                if i < len(ch):
                    c = ch[i]
                    iter_idx[v] = i + 1
                    cst = status.get(c, 0)
                    if cst == 0:
                        stack.append(c)
                    else:
            
                        pass
                else:

                    best_child = -1
                    best_d = 0
                    for c in ch:
                        cst = status.get(c, 0)
                        if cst == 2:
                            d = self.depth[c]
                        elif cst == 1:
                            d = 1
                        else:
                            d = 0
                        if d > best_d:
                            best_d = d
                            best_child = c
                    self.depth[v] = best_d + 1
                    self.max_child[v] = best_child
                    status[v] = 2
                    stack.pop()

            else: 
                stack.pop()

    def _compute_all_depths(self) -> None:
        status: Dict[int, int] = {}       
        ch_cache: Dict[int, List[int]] = {}

        for v in sorted(self.kmers.keys()):
            if status.get(v, 0) == 0:
                self._compute_depth_from(v, status, ch_cache)

    def _get_longest_path(self) -> List[int]:

        self._compute_all_depths()

        start = -1
        max_d = 0
        for v in self.kmers.keys():
            d = self.depth[v]
            if d > max_d:
                max_d = d
                start = v

        path: List[int] = []
        while start != -1:
            path.append(start)
            start = self.max_child[start]
        return path

    def _delete_path(self, path: List[int]) -> None:
        dead: Set[int] = set(path)

        for idx in path:
            if idx in self.children: del self.children[idx]
            if idx in self.counts:   del self.counts[idx]
            if idx in self.kmers:    del self.kmers[idx]
            if idx in self.depth:    del self.depth[idx]
            if idx in self.max_child:del self.max_child[idx]
        for idx in list(self.children.keys()):
            if self.children[idx]:
                self.children[idx] = self.children[idx].difference(dead)

    def _concat_path(self, path: List[int]) -> str:
        if not path:
            return ""
        s = self.kmers[path[0]]
        for i in range(1, len(path)):
            s += self.kmers[path[i]][-1]
        return s

    def get_longest_contig(self) -> str:
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig