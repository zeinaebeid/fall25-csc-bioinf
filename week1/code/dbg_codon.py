import copy
from typing import List, Dict

def reverse_complement(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    r = list(seq[::-1])
    for i in range(len(r)):
        r[i] = comp.get(r[i], r[i])
    return "".join(r)

class Node:
    __slots__ = ("_children", "_count", "kmer", "visited", "depth", "max_depth_child")

    def __init__(self, kmer: str):
        self._children = set()       # set[int]
        self._count = 0
        self.kmer = kmer
        self.visited = False
        self.depth = 0
        self.max_depth_child = -1    # -1 means “none”

    def add_child(self, idx_child: int):
        self._children.add(idx_child)

    def increase(self):
        self._count += 1

    def reset(self):
        self.visited = False
        self.depth = 0
        self.max_depth_child = -1

    def get_count(self) -> int:
        return self._count

    def get_children(self) -> List[int]:
        return list(self._children)

    def remove_children(self, target_set: set):
        if self._children:
            self._children = self._children.difference(target_set)

class DBG:
    def __init__(self, k: int, data_list: List[List[str]]):
        self.k = k
        self.nodes: Dict[int, Node] = {}
        self.kmer2idx: Dict[str, int] = {}
        self.kmer_count = 0
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list: List[List[str]]):
        assert len(data_list) > 0
        any_ok = any(len(read) >= self.k for group in data_list for read in group)
        assert any_ok, "No reads long enough for k"

    def _build(self, data_list: List[List[str]]):
        for reads in data_list:
            for seq in reads:
                if len(seq) < self.k:
                    continue
                rc = reverse_complement(seq)
                limit = len(seq) - self.k
                for i in range(limit):
                    k1 = seq[i:i+self.k]
                    k2 = seq[i+1:i+1+self.k]
                    if len(k1) == self.k and len(k2) == self.k:
                        self._add_arc(k1, k2)
                    rk1 = rc[i:i+self.k]
                    rk2 = rc[i+1:i+1+self.k]
                    if len(rk1) == self.k and len(rk2) == self.k:
                        self._add_arc(rk1, rk2)

    def _add_node(self, kmer: str) -> int:
        if not kmer or len(kmer) != self.k:
            return -1
        idx = self.kmer2idx.get(kmer, -1)
        if idx == -1:
            idx = self.kmer_count
            self.kmer2idx[kmer] = idx
            self.nodes[idx] = Node(kmer)
            self.kmer_count += 1
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1: str, kmer2: str):
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        if idx1 != -1 and idx2 != -1:
            self.nodes[idx1].add_child(idx2)

    def _get_count(self, idx: int) -> int:
        return self.nodes[idx].get_count()

    def _get_sorted_children(self, idx: int) -> List[int]:
        kids = self.nodes[idx].get_children()
        kids.sort(key=self._get_count, reverse=True)
        return kids

    def _get_depth(self, idx: int) -> int:
        node = self.nodes[idx]
        if not node.visited:
            node.visited = True
            max_d, max_c = 0, -1
            for child in self._get_sorted_children(idx):
                d = self._get_depth(child)
                if d > max_d:
                    max_d, max_c = d, child
            node.depth = max_d + 1
            node.max_depth_child = max_c
        return node.depth

    def _reset(self):
        for idx in list(self.nodes.keys()):
            self.nodes[idx].reset()

    def _get_longest_path(self) -> List[int]:
        max_d, max_idx = 0, -1
        for idx in self.nodes.keys():
            d = self._get_depth(idx)
            if d > max_d:
                max_d, max_idx = d, idx
        path: List[int] = []
        while max_idx != -1:
            path.append(max_idx)
            max_idx = self.nodes[max_idx].max_depth_child
        return path

    def _delete_path(self, path: List[int]):
        path_set = set(path)
        for idx in list(self.nodes.keys()):
            self.nodes[idx].remove_children(path_set)
        for idx in path:
            if idx in self.nodes:
                del self.nodes[idx]

    def _concat_path(self, path: List[int]) -> str:
        if not path:
            return ""
        seq = copy.copy(self.nodes[path[0]].kmer)
        for i in range(1, len(path)):
            seq += self.nodes[path[i]].kmer[-1]
        return seq

    def get_longest_contig(self) -> str:
        if not self.nodes:
            return ""
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        if contig != "":
            self._delete_path(path)
        return contig
