from typing import List, Dict, Set, Optional, Tuple

def reverse_complement(key: str) -> str:
    complement: Dict[str, str] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    key_rev: List[str] = list(key[::-1])
    for i in range(len(key_rev)):
        key_rev[i] = complement.get(key_rev[i], key_rev[i])
    return ''.join(key_rev)

class Node:
    __slots__ = ("_children", "_count", "kmer", "visited", "depth", "max_depth_child")

    def __init__(self, kmer: str):
        self._children: Set[int] = set()
        self._count: int = 0
        self.kmer: str = kmer
        self.visited: bool = False
        self.depth: int = 0
        self.max_depth_child: int = -1

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

    def remove_children(self, target_set: Set[int]):
        if self._children:
            self._children = self._children.difference(target_set)

class DBG:
    def __init__(self, k: int, data_list: List[List[str]]):
        self.k: int = k
        self.nodes: Dict[int, Node] = {}
        self.kmer2idx: Dict[str, int] = {}
        self.kmer_count: int = 0
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list: List[List[str]]):
        assert len(data_list) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list: List[List[str]]):
        for data in data_list:
            for original in data:
                rc: str = reverse_complement(original)
                for i in range(len(original) - self.k):
                    self._add_arc(original[i:i + self.k], original[i + 1:i + 1 + self.k])
                    self._add_arc(rc[i:i + self.k], rc[i + 1:i + 1 + self.k])

    def _add_node(self, kmer: str) -> int:
        idx: int = self.kmer2idx.get(kmer, -1)
        if idx == -1:
            idx = self.kmer_count
            self.kmer2idx[kmer] = idx
            self.nodes[idx] = Node(kmer)
            self.kmer_count += 1
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1: str, kmer2: str):
        idx1: int = self._add_node(kmer1)
        idx2: int = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child_idx: int) -> int:
        return self.nodes[child_idx].get_count()

    def _get_sorted_children(self, idx: int) -> List[int]:
        children: List[int] = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx: int) -> int:
        node: Node = self.nodes[idx]
        if not node.visited:
            node.visited = True
            children: List[int] = self._get_sorted_children(idx)
            max_depth: int = 0
            max_child: int = -1
            for child in children:
                d: int = self._get_depth(child)
                if d > max_depth:
                    max_depth = d
                    max_child = child
            node.depth = max_depth + 1
            node.max_depth_child = max_child
        return node.depth

    def _reset(self):
        for idx in list(self.nodes.keys()):
            self.nodes[idx].reset()

    def _get_longest_path(self) -> List[int]:
        max_depth: int = 0
        max_idx: int = -1
        for idx in self.nodes.keys():
            d: int = self._get_depth(idx)
            if d > max_depth:
                max_depth = d
                max_idx = idx

        path: List[int] = []
        while max_idx != -1:
            path.append(max_idx)
            max_idx = self.nodes[max_idx].max_depth_child
        return path

    def _delete_path(self, path: List[int]):
        for idx in path:
            if idx in self.nodes:
                del self.nodes[idx]
        path_set: Set[int] = set(path)
        for idx in self.nodes.keys():
            self.nodes[idx].remove_children(path_set)

    def _concat_path(self, path: List[int]) -> str:
        if not path:
            return ""
        
        concat: str = self.nodes[path[0]].kmer
        for i in range(1, len(path)):
            concat += self.nodes[path[i]].kmer[-1]
        return concat

    def get_longest_contig(self) -> Optional[str]:
        self._reset()
        path: List[int] = self._get_longest_path()
        contig: str = self._concat_path(path)
        self._delete_path(path)
        return contig