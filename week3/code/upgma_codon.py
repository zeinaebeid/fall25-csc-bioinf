# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.


from typing import List, Optional
from tree_codon import Tree, TreeNode
from numpy import ndarray, array, zeros, ones, full, float32, float64, int32, uint8, uint32, asarray
import numpy as np

MAX_FLOAT: float32 = float32(3.4028235e+38)

def upgma(distances: ndarray) -> Tree:

    # Convert to Codon-compatible NumPy array (ndarray[float, 2])
    distances = asarray(distances).astype(float64, copy=True)

    # validation checks
    if distances.shape[0] != distances.shape[1]:
        raise ValueError("Distance matrix must be square")
    if not np.allclose(distances, distances.T):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances < 0).any():
        raise ValueError("Distances must be non-negative")

    n = distances.shape[0]

    # Nodes list must allow None assignments
    nodes: List[TreeNode] = [TreeNode(index=i) for i in range(n)]
    nodes = array(nodes) 

    # Boolean and numeric arrays
    is_clustered = np.zeros(n, dtype=bool)
    cluster_size = np.ones(n, dtype=np.int32)
    node_heights = np.zeros(n, dtype=np.float64)



   # Copy distance matrix as float array
    D = distances.copy().astype(float64, copy=True)

    while True:
        dist_min: float = float("inf")
        i_min: int = -1
        j_min: int = -1

        # find minimum distance among active clusters
        for i in range(n):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                if D[i, j] < dist_min:
                    dist_min = D[i, j]
                    i_min = i
                    j_min = j

        # Stop when no more pairs to merge
        if i_min == -1 or j_min == -1:
            break

        height: float = dist_min / 2.0

        node_i = nodes[i_min]
        node_j = nodes[j_min]

        if node_i is None or node_j is None:
            continue  # not raise — this keeps flow same as the other Codon version

        nodes[i_min] = TreeNode(
            [node_i, node_j],
            [height - node_heights[i_min], height - node_heights[j_min]]
        )

        node_heights[i_min] = height

        # mark j_min as inactive
        nodes[j_min] = None
        is_clustered[j_min] = True

        # update distances using weighted mean
        total_size = cluster_size[i_min] + cluster_size[j_min]
        for k in range(n):
            if not is_clustered[k] and k != i_min:
                mean = (
                    D[i_min, k] * float(cluster_size[i_min]) +
                    D[j_min, k] * float(cluster_size[j_min])
                ) / float(total_size)
                D[i_min, k] = mean
                D[k, i_min] = mean

        cluster_size[i_min] = total_size

    # return the last non-None node as root
    for root in reversed(nodes):
        if root is not None:
            return Tree(root)

    raise ValueError("Failed to build UPGMA tree: input may be empty")


