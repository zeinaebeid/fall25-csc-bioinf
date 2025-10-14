# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import copy
import numpy as np


class Tree: 

    def __init__(self, root):
        root.as_root()
        self._root = root


        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)
        indices = [leaf.index for leaf in leaves_unsorted]

        self._leaves = [None] * leaf_count
        for i, idx in enumerate(indices):
            if idx >= leaf_count or idx < 0:
                raise TreeError("The tree's indices are out of range")
            self._leaves[idx] = leaves_unsorted[i]


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





class TreeNode:
    """
    Python port of Biotite's TreeNode (simplified).

    Supports: ** insert what it supports 
    """
    def __init__(self, children=None, distances=None, index=None):
        self._is_root = False
        self._distance = 0.0
        self._parent = None
        self._index = -1   
        self._children = None


        if index is None:
            if children is None or distances is None:
                raise TypeError(
                    "Either a reference index (for leaf) or "
                    "children+distances (for intermediate) must be set."
                )

            if len(children) == 0:
                raise TreeError("Intermediate nodes must have at least one child.")
            if len(children) != len(distances):
                raise ValueError("Number of children must equal number of distances.")

            for item in children:
                if not isinstance(item, TreeNode):
                    raise TypeError(f"Expected 'TreeNode', got '{type(item).__name__}'")
            for item in distances:
                if not isinstance(item, (float, int)):
                    raise TypeError(f"Expected 'float' or 'int', got '{type(item).__name__}'")

            if len(set(id(c) for c in children)) != len(children):
                raise TreeError("Two child nodes cannot be the same object.")

            self._index = -1
            self._children = tuple(children)

            for child, distance in zip(children, distances):
                child._set_parent(self, distance)

      
        elif index < 0:
            raise ValueError("Index cannot be negative.")
        else:
            if children is not None or distances is not None:
                raise TypeError(
                    "Reference index and child nodes are mutually exclusive."
                )
            self._index = index
            self._children = None




    
    def _set_parent(self, parent, distance):
        """Link this node to its parent and store branch length."""
        if self._parent is not None or self._is_root:
            raise TreeError("Node already has a parent.")
        self._parent = parent
        self._distance = float(distance)

    def as_root(self):
        """
        Convert this node into a root node.

        Raises TreeError if node already has a parent.
        """
        if self._parent is not None:
            raise TreeError("Node has a parent; cannot be root.")
        self._is_root = True

    def get_leaves(self):
        """
        Return all leaf nodes under (or including) this node.
        """
        leaf_list = []
        self._get_leaves(leaf_list)
        return leaf_list


    # bring get leaves into the class 
    def _get_leaves(self, leaf_list):
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
            raise TreeError("The nodes do not have a common ancestor")

        # Traverse upward from self to the LCA
        current_node = self
        while current_node is not lca:
            distance += 1 if topological else current_node._distance
            current_node = current_node._parent

        # Traverse upward from node to the LCA
        current_node = node
        while current_node is not lca:
            distance += 1 if topological else current_node._distance
            current_node = current_node._parent

        return distance

    
    def lowest_common_ancestor(self, node):
        """
        Return the lowest common ancestor (LCA) of this node and another.
        """
        # Build both ancestor paths
        self_path = _create_path_to_root(self)
        other_path = _create_path_to_root(node)

        # compare paths from the root downward (reverse order)
        lca = None
        for i in range(-1, -min(len(self_path), len(other_path)) - 1, -1):
            if self_path[i] is other_path[i]:
                lca = self_path[i]
            else:
                break
        return lca


    # reads entire newick string 
    # rebuilds node structure recursively 
    # keep docstring for codon conversion 
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
        # Remove all whitespace
        newick = "".join(newick.split())

        subnewick_start_i = -1
        subnewick_stop_i = -1

        # Find first '(' and matching ')'
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

        # ----------------------------------------
        # CASE 1: No parentheses → leaf node
        # ----------------------------------------
        if subnewick_start_i == -1 and subnewick_stop_i == -1:
            label_and_distance = newick
            try:
                label, distance = label_and_distance.split(":")
                distance = float(distance)
            except ValueError:
                # No colon provided → default distance = 0
                label = label_and_distance
                distance = 0.0

            # Convert label to index
            label = label.strip()
            if label == "":
                raise ValueError("Leaf node label missing")

            try:
                index = int(label)
            except ValueError:
                if labels is None:
                    raise ValueError(f"Label '{label}' cannot be parsed as int and no labels provided")
                index = labels.index(label)

            return TreeNode(index=index), distance

        # ----------------------------------------
        # CASE 2: Internal node (has children)
        # ----------------------------------------
        else:
            # Parse possible label/distance after closing ')'
            if subnewick_stop_i == len(newick):
                label = None
                distance = 0.0
            else:
                label_and_distance = newick[subnewick_stop_i:]
                try:
                    label, distance = label_and_distance.split(":")
                    distance = float(distance)
                except ValueError:
                    label = label_and_distance
                    distance = 0.0
                # Intermediate labels are discarded — only distance matters
                distance = float(distance)

            subnewick = newick[subnewick_start_i + 1 : subnewick_stop_i - 1]
            if len(subnewick) == 0:
                raise ValueError("Intermediate node must have at least one child")

            # Split child sections at top-level commas
            comma_pos = []
            level = 0
            for i, char in enumerate(subnewick):
                if char == "(":
                    level += 1
                elif char == ")":
                    level -= 1
                elif char == "," and level == 0:
                    comma_pos.append(i)
                if level < 0:
                    raise ValueError("Bracket closed before it was opened")

            # Recursive construction of child nodes
            children = []
            distances = []
            prev_pos = 0
            for pos in comma_pos + [len(subnewick)]:
                part = subnewick[prev_pos:pos]
                child, dist = TreeNode.from_newick(part, labels=labels)
                children.append(child)
                distances.append(dist)
                prev_pos = pos + 1

            return TreeNode(children, distances), distance


    def __eq__(self, other):
        """Return True if two TreeNode objects are structurally identical."""
        if not isinstance(other, TreeNode):
            return False

        # Compare branch length (distance to parent)
        if self._distance != other._distance:
            return False

        # Leaf node comparison
        if self._index != -1:
            return self._index == other._index

        # Internal node comparison (compare sets of children)
        return frozenset(self._children) == frozenset(other._children)


    def __hash__(self):
        """
        Return a hash value for this TreeNode.

        The hash is based on:
        - the node's index (for leaves),
        - the set of its children (for internal nodes),
        - and the distance to its parent.

        Order of children does not affect the hash.
        """
        if self._children is not None:
            # Children order doesn’t matter → convert to frozenset
            children_set = frozenset(self._children)
        else:
            children_set = None

        return hash((self._index, children_set, self._distance))







# helper
def _create_path_to_root(node):
    """
    Return a list of nodes representing the path from the given node
    up to the root.
    """
    path = []
    current_node = node
    while current_node is not None:
        path.append(current_node)
        current_node = current_node._parent
    return path




class TreeError(Exception):
    """
    An exception that occurs in context of tree topology.
    """
    pass