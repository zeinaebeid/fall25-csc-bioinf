# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.


import copy
import numpy as np
from typing import List


@extend
class set:
    def __hash__(self):
        MAX = int.MAX
        MASK = 2 * MAX + 1
        n = len(self)
        h = 1927868237 * (n + 1)
        h &= MASK
        for x in self:
            hx = hash(x)
            h ^= (hx ^ (hx << 16) ^ 89869747) * 3644798167
            h &= MASK
        h = h * 69069 + 907133923
        h &= MASK
        if h > MAX:
            h -= MASK + 1
        if h == -1:
            h = 590923713
        return h


class TreeNode:
    """
    Python port of Biotite's TreeNode (simplified).

    Supports: ** insert what it supports
    """
    _is_root: bool
    _distance: float
    _parent: Optional[TreeNode]
    _index: int
    _children: List[TreeNode]

    def __init__(self, children=None, distances=None, index=None):
        self._is_root = False
        self._distance = 0.0
        self._parent = None
        self._index = -1   
        self._children = []

        # CASE 1: Intermediate (internal) node
        if index is None:
            if children is None or distances is None:
                raise TypeError(
                    "Either a reference index (for leaf) or "
                    "children+distances (for intermediate) must be set."
                )

            if len(children) == 0:
                raise ValueError("Intermediate nodes must have at least one child.")
            if len(children) != len(distances):
                raise ValueError("Number of children must equal number of distances.")

            for item in children:
                if not isinstance(item, TreeNode):
                    raise TypeError(f"Expected 'TreeNode', got '{type(item).__name__}'")
            for item in distances:
                if not isinstance(item, (float, int)):
                    raise TypeError(f"Expected 'float' or 'int', got '{type(item).__name__}'")

            if len(set(id(c) for c in children)) != len(children):
                raise ValueError("Two child nodes cannot be the same object.")

            self._index = -1
            # store as list for field type consistency
            self._children = [i for i in children]

            for child, distance in zip(children, distances):
                child._set_parent(self, float(distance))

        # CASE 2: Leaf node
        elif index < 0:
            raise ValueError("Index cannot be negative.")
        else:
            if children is not None or distances is not None:
                raise TypeError(
                    "Reference index and child nodes are mutually exclusive."
                )
            self._index = index
            self._children = []

  
    def _set_parent(self, parent, distance):
        """Link this node to its parent and store branch length."""
        if self._parent is not None or self._is_root:
            raise ValueError("Node already has a parent.")
        self._parent = parent
        self._distance = float(distance)

    def as_root(self):
        """
        Convert this node into a root node.

        Raises TreeError if node already has a parent.
        """
        if self._parent is not None:
            raise ValueError("Node has a parent; cannot be root.")
        self._is_root = True

    def get_leaves(self):
        """
        Return all leaf nodes under (or including) this node.
        """
        leaf_list = []
        self._get_leaves(leaf_list)
        return leaf_list

    
    def _get_leaves(self: TreeNode, leaf_list: List[TreeNode]) -> None:
        """
        Recursively collect all leaf nodes under the given node.
        """
        if self._index == -1:
            for child in self._children:
                child._get_leaves(leaf_list)
        else:
            leaf_list.append(self)


    
    @property
    def index(self):
        """Return the leaf index if leaf, else None."""
        return None if self._index == -1 else self._index

    @property
    def children(self):
        """Return the tuple of child nodes (None if leaf)."""
        return self._children

    @property
    def parent(self):
        """Return the parent node (None if root)."""
        return self._parent

    @property
    def distance(self):
        """Return distance to parent, or None if root."""
        return None if self._parent is None else self._distance

    def distance_to(self, node, topological=False):
        """
        Return the distance between this node and another node.
        """
        distance = 0.0
        lca = self.lowest_common_ancestor(node)

        if lca is None:
            raise ValueError("The nodes do not have a common ancestor")

        
        current_node = self
        while current_node is not lca:
            distance += 1 if topological else current_node._distance
            current_node = current_node._parent

    
        current_node = node
        while current_node is not lca:
            distance += 1 if topological else current_node._distance
            current_node = current_node._parent

        return distance

    def lowest_common_ancestor(self, node):
      
        self_path = _create_path_to_root(self)
        other_path = _create_path_to_root(node)
        lca = None
        for i in range(-1, -min(len(self_path), len(other_path)) - 1, -1):
            if self_path[i] is other_path[i]:
                lca = self_path[i]
            else:
                break
        return lca


    @staticmethod
    def from_newick(newick, labels=None):
        """
        Create a TreeNode (and all its child nodes) from a Newick notation.

        Parameters
        ----------
        newick : str
            The Newick notation to create the node from.
        labels : list of str, optional
            If the Newick notation contains non-integer labels,
            this list maps labels to reference indices. The index
            corresponds to the labels position in the list.

        Returns
        -------
        (TreeNode, float)
            The parsed TreeNode and its distance to the parent.
            Distance defaults to 0 if not specified.

        Notes
        -----
        - The string must NOT have a terminal semicolon.
        - Labels on intermediate nodes are ignored.
        """
    
        newick = "".join(newick.split())

        subnewick_start_i = -1
        subnewick_stop_i = -1

    
        for i, char in enumerate(newick):
            if char == "(":
                subnewick_start_i = i
                break
            if char == ")":
                raise ValueError("Bracket closed before it was opened")

        for i in range(len(newick) - 1, -1, -1):
            char = newick[i]
            if char == ")":
                subnewick_stop_i = i + 1
                break
            if char == "(":
                raise ValueError("Bracket opened but not closed")


      
        if subnewick_start_i == -1 and subnewick_stop_i == -1:
            label_and_distance: str = newick
            label: str = ""
            distance: float = 0.0
            co: int = -1
    

            for i in range(len(label_and_distance)):
                if label_and_distance[i] == ":":
                    col = i
                    break

            if col != -1:
                label = label_and_distance[:col]
                dist_str: str = label_and_distance[col + 1:]
                distance = float(dist_str)
            else:
                label = label_and_distance

            index: int = -1
            if labels is None:
                index = int(label)
            else:
                index = labels.index(label)

            return TreeNode(index=index), distance



        else:
        
            label: str = ""
            label_and_distance: str = ""
            distance: float = 0.0
            subnewick: str = ""
            comma_pos: List[int] = []
            level: int = 0
            children: List[TreeNode] = []
            distances: List[float] = []
            prev_pos: int = 0
            col: int = -1

      
            if subnewick_stop_i == len(newick):
                label = ""
                distance = 0.0
            else:
                label_and_distance = newick[subnewick_stop_i:]
                label = ""
                distance = 0.0

        
                for i in range(len(label_and_distance)):
                    if label_and_distance[i] == ":":
                        col= i
                        break

                if col != -1:
                    label = label_and_distance[:col]
                    dist_str: str = label_and_distance[col + 1:]
                    distance = float(dist_str)
                else:
                    label = label_and_distance
                    distance = 0.0


            subnewick = newick[subnewick_start_i + 1 : subnewick_stop_i - 1]
            if len(subnewick) == 0:
                raise ValueError("Intermediate node must have at least one child")


            level = 0
            for i in range(len(subnewick)):
                char = subnewick[i]
                if char == "(":
                    level += 1
                elif char == ")":
                    level -= 1
                elif char == "," and level == 0:
                    comma_pos.append(i)
                if level < 0:
                    raise ValueError("Bracket closed before it was opened")


            prev_pos = 0
            for pos in comma_pos + [len(subnewick)]:
                part: str = subnewick[prev_pos:pos]
                child, dist = TreeNode.from_newick(part, labels=labels)
                children.append(child)
                distances.append(dist)
                prev_pos = pos + 1

            return TreeNode(children, distances), distance


    def __eq__(self, other):
        if not isinstance(other, TreeNode):
            return False

        if self._distance != other._distance:
            return False

        if self._index != -1:
            return self._index == other._index


        return set(self._children) == set(other._children)


    def __ne__(self, other):
        """Codon requires explicit inequality for set comparison."""
        return not self.__eq__(other)

    def __hash__(self):
        if self._children is not None:
            children_set = set(self._children)
        else:
            children_set = None
        return hash((self._index, children_set, self._distance))



class Tree:
    _root: TreeNode
    _leaves: List[TreeNode]

    def __init__(self, root):
        root.as_root()
        self._root = root

        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)
        indices = [leaf.index for leaf in leaves_unsorted]

       
        leaves_filled = [root] * leaf_count
        for i, idx in enumerate(indices):
            if idx >= leaf_count or idx < 0:
                raise ValueError("The tree's indices are out of range")
            leaves_filled[idx] = leaves_unsorted[i]
        self._leaves = leaves_filled

    def get_distance(self, index1, index2, topological=False):
        return self._leaves[index1].distance_to(self._leaves[index2], topological)

    @property
    def root(self):
        return self._root

    @property
    def leaves(self):
        return list(self._leaves)

    @staticmethod
    def from_newick(newick, labels=None):
        """
        Create a TreeNode (and subtree) from a Newick notation string.
        """
        newick = newick.strip()
        if len(newick) == 0:
            raise ValueError("Newick string is empty")


        if newick.endswith(";"):
            newick = newick[:-1]

        root, _ = TreeNode.from_newick(newick, labels)
        return Tree(root)

    def __len__(self):
        return len(self._leaves)

    def __eq__(self, item):
        if not isinstance(item, Tree):
            return False
        return self._root == item._root

    def __hash__(self):
        return hash(self._root)




def _create_path_to_root(node: TreeNode) -> List[TreeNode]:
    """
    Return a list of nodes representing the path from the given node
    up to the root.
    """
    path: List[TreeNode] = []
    current_node: Optional[TreeNode] = node
    while current_node is not None:
        path.append(current_node)
        current_node = current_node._parent
    return path