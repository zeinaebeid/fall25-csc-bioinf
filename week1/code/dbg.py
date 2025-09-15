import copy

def reverse_complement(key: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    key_rev = list(key[::-1])
    for i in range(len(key_rev)):
        key_rev[i] = complement.get(key_rev[i], key_rev[i])
    return ''.join(key_rev)

class Node:
    __slots__ = ("_children", "_count", "kmer", "visited", "depth", "max_depth_child")

    def __init__(self, kmer: str):
        # children are indices (ints) into DBG.nodes
        self._children = set()       # set[int]
        self._count = 0              # int
        self.kmer = kmer             # str
        self.visited = False         # bool
        self.depth = 0               # int
        self.max_depth_child = -1    # int sentinel instead of None

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

    def get_children(self):
        # return a list[int] for deterministic sorting
        return list(self._children)

    def remove_children(self, target_set):
        # remove any edges to nodes in target_set
        if self._children:
            self._children = self._children.difference(target_set)

class DBG:
    def __init__(self, k: int, data_list):
        self.k = k
        self.nodes = {}        # dict[int, Node]
        self.kmer2idx = {}     # dict[str, int]
        self.kmer_count = 0
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list):
        assert len(data_list) > 0
        # allow different read lengths as long as they’re >= k
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list):
        for data in data_list:
            for original in data:
                rc = reverse_complement(original)
                # Correct loop bound: last i should allow slices [i:i+k] and [i+1:i+1+k]
                # so i in [0 .. len - k - 1]
                for i in range(len(original) - self.k):
                    self._add_arc(original[i:i + self.k], original[i + 1:i + 1 + self.k])
                    self._add_arc(rc[i:i + self.k],        rc[i + 1:i + 1 + self.k])

    def _add_node(self, kmer: str) -> int:
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
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child_idx: int) -> int:
        return self.nodes[child_idx].get_count()

    def _get_sorted_children(self, idx: int):
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx: int) -> int:
        node = self.nodes[idx]
        if not node.visited:
            node.visited = True
            children = self._get_sorted_children(idx)
            max_depth = 0
            max_child = -1
            for child in children:
                d = self._get_depth(child)
                if d > max_depth:
                    max_depth = d
                    max_child = child
            node.depth = max_depth + 1
            node.max_depth_child = max_child
        return node.depth

    def _reset(self):
        # iterating over list(self.nodes.keys()) guards against size changes elsewhere
        for idx in list(self.nodes.keys()):
            self.nodes[idx].reset()

    def _get_longest_path(self):
        max_depth = 0
        max_idx = -1
        for idx in self.nodes.keys():
            d = self._get_depth(idx)
            if d > max_depth:
                max_depth = d
                max_idx = idx

        path = []
        while max_idx != -1:
            path.append(max_idx)
            max_idx = self.nodes[max_idx].max_depth_child
        return path

    def _delete_path(self, path):
        # delete nodes on the path
        for idx in path:
            if idx in self.nodes:
                del self.nodes[idx]
        # remove edges pointing to any deleted node
        path_set = set(path)
        for idx in self.nodes.keys():
            self.nodes[idx].remove_children(path_set)

    def _concat_path(self, path):
        if len(path) < 1:
            return ""
        concat = copy.copy(self.nodes[path[0]].kmer)
        for i in range(1, len(path)):
            concat += self.nodes[path[i]].kmer[-1]
        return concat

    def get_longest_contig(self) -> str:
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig
