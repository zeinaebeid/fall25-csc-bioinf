# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.


import time
from tree_codon import Tree, TreeNode
from nj_codon import neighbor_joining
from upgma_codon import upgma

import numpy as np
from numpy import ndarray, array, zeros, ones, full, float32, float64, int32, uint8, uint32, asarray
from typing import List, Optional



def get_upgma_newick() -> str:
    with open("../data/newick_upgma.txt", "r") as file:
        newick: str = file.read().strip()
    return newick


def get_tree(distances: ndarray) -> Tree:
    return upgma(distances)


def get_distances() -> ndarray:
    rows: List[List[int]] = []
    with open("../data/distances.txt", "r") as fh:
        for ln in fh:
            s = ln.strip()
            if not s:
                continue
            rows.append([int(tok) for tok in s.split()])

    # Build an int32 numpy array with the parsed shape
    n = len(rows)
    m = len(rows[0]) if n else 0
    mat = zeros((n, m), dtype=int32)

    for r in range(n):
        vals = rows[r]
        for c in range(m):
            mat[r, c] = int32(vals[c])

    return mat



def test_distances(tree):
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes

    dist = tree.root.distance_to(tree.leaves[0])
    # Explicitly ensure leaves is typed as a list of TreeNode
    leaves: List[TreeNode] = tree.leaves

    for leaf in leaves:
        assert leaf.distance_to(tree.root) == dist

    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10



def test_upgma(tree, upgma_newick):
    """
    Compare the results of `upgma()` with DendroUPGMA.
    """
    ref_tree = Tree.from_newick(upgma_newick)

    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    n: int = len(tree)
    for i in range(n):
        for j in range(n):

            tree_dist = tree.get_distance(i, j)
            ref_dist = ref_tree.get_distance(i, j)

            # manual approx with tolerance 1e-3
            diff: float = abs(tree_dist - ref_dist)
            if diff >= 1e-3:
                # Codon does not support %-formatting or .format()
                raise AssertionError(
                    f"Distance mismatch at ({i}, {j}): {tree_dist} vs {ref_dist}"
                )

            tree_topo = tree.get_distance(i, j, topological=True)
            ref_topo = ref_tree.get_distance(i, j, topological=True)

            if tree_topo != ref_topo:
                raise AssertionError(
                    f"Topological distance mismatch at ({i}, {j}): {tree_topo} vs {ref_topo}"
                )




def test_neighbor_joining():
    """
    Compare the results of `neighbor_join()` with a known tree.
    """
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ])  # fmt: skip

    ref_tree = Tree(
        TreeNode(
            [
                TreeNode(
                    [
                        TreeNode(
                            [
                                TreeNode(index=0),
                                TreeNode(index=1),
                            ],
                            [1, 4],
                        ),
                        TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                TreeNode(
                    [
                        TreeNode(index=3),
                        TreeNode(index=4),
                    ],
                    [3, 2],
                ),
                TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )

    test_tree = neighbor_joining(dist)

    assert test_tree == ref_tree



def main():

    distances = get_distances()
    newick = get_upgma_newick()
    tree = get_tree(distances)
    
    start_time = time.time()

    test_distances(tree)
    print("Distance Tree Passed\n")
    test_neighbor_joining()
    print("Neighbor Joining Passed\n")
    test_upgma(tree, newick)
    print("Upgma Passed")
    

    end_time = time.time()
    elapsed_ms = (end_time - start_time) * 1000
    
    print(f"codon:  {elapsed_ms:.0f}ms")

if __name__ == "__main__":
    main()