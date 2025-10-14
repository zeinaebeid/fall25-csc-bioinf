# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import numpy as np
from .tree import Tree, TreeNode

def neighbor_joining(distances):
    """
    neighbor_joining(distances)

    Perform hierarchical clustering using the Neighbor-Joining algorithm.
    Returns a rooted tree whose leaves' indices correspond to rows/cols
    of the input distance matrix.

    Parameters
    ----------
    distances : np.ndarray, shape (n, n)
        Symmetric, non-negative pairwise distance matrix (zeros on diag).

    Returns
    -------
    Tree
        Rooted tree (root will have 3 children; others binary).

    Raises
    ------
    ValueError
        If matrix is not square/symmetric, has NaN/negative values, or n < 4.
    """

    # validation
    D = np.asarray(distances, dtype=float)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("Distance matrix must be square")
    if not np.allclose(D, D.T):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(D).any():
        raise ValueError("Distance matrix contains NaN values")
    if (D < 0).any():
        raise ValueError("Distances must be non-negative")
    n = D.shape[0]
    if n < 4:
        raise ValueError("At least 4 nodes are required")


    nodes = [TreeNode(index=i) for i in range(n)]     # active node objects
    is_clustered = np.zeros(n, dtype=bool)            # retired rows/cols
    divergence = np.zeros(n, dtype=float)             # r_i = sum_k d(i,k)
    Q = np.zeros_like(D)                              # corrected distances

    # inverts the bool array 
    def remaining_count():
        return int((~is_clustered).sum())

    # --- main loop ---
    while True:
        m = remaining_count()
        # compute divergences r_i over active indices
        for i in range(n):
            if is_clustered[i]:
                continue
            s = 0.0
            row = D[i]
            for k in range(n):
                if not is_clustered[k]:
                    s += row[k]
            divergence[i] = s

        # compute Q-matrix (corrected distances) on active upper triangle
        Q.fill(0.0)
        for i in range(n):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                Q[i, j] = (m - 2) * D[i, j] - divergence[i] - divergence[j]

        # find pair with minimum corrected distance
        dist_min = float("inf")
        i_min = j_min = -1
        for i in range(n):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                q = Q[i, j]
                if q < dist_min:
                    dist_min = q
                    i_min, j_min = i, j

        # done if no active pair found
        if i_min == -1 or j_min == -1:
            break

        # limb lengths from NJ formula
        dij = D[i_min, j_min]
        denom = (m - 2)
        node_dist_i = 0.5 * (dij + (divergence[i_min] - divergence[j_min]) / denom)
        node_dist_j = 0.5 * (dij + (divergence[j_min] - divergence[i_min]) / denom)

        if m > 3:
            # create an internal node with two children
            nodes[i_min] = TreeNode([nodes[i_min], nodes[j_min]],
                                    [node_dist_i, node_dist_j])
            # retire j_min
            nodes[j_min] = None
            is_clustered[j_min] = True

            # update distances from new cluster at i_min to other active nodes
            for k in range(n):
                if not is_clustered[k] and k != i_min:
                    dik = D[i_min, k]
                    djk = D[j_min, k]
                    newd = 0.5 * (dik + djk - dij)
                    D[i_min, k] = D[k, i_min] = newd
        else:
            # final step: combine the remaining third node k into the root
            is_clustered[i_min] = True
            is_clustered[j_min] = True
            # the only remaining active index
            k = int(np.where(~is_clustered)[0][0])
            node_dist_k = 0.5 * (D[i_min, k] + D[j_min, k] - dij)
            root = TreeNode([nodes[i_min], nodes[j_min], nodes[k]],
                            [node_dist_i, node_dist_j, node_dist_k])
            return Tree(root)

    # Fallback (should not happen for valid n>=4 inputs):
    for root in reversed(nodes):
        if root is not None:
            return Tree(root)
    raise ValueError("Failed to construct NJ tree")