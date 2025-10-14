# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.


from tree_codon import Tree, TreeNode, ValueError
import numpy as np
from numpy import ndarray, array, zeros, full, float32, float64, asarray, uint8
from typing import List


MAX_FLOAT: float32 = float32(3.4028235e+38)



def neighbor_joining(distances: ndarray) -> Tree:

    D = asarray(distances).astype(float32, copy=True)

    rows: int = D.shape[0]
    cols: int = D.shape[1]

    # validation checks

    # check matrix shape
    if rows != cols:
        raise ValueError("Distance matrix must be square (NxN).")

    is_symmetric: bool = True
    has_nan: bool = False
    has_inf: bool = False
    has_negative: bool = False

    for i in range(rows):
        for j in range(cols):
            val_ij: float32 = D[i, j]
            val_ji: float32 = D[j, i]

            # Symmetry check (within tolerance)
            if abs(val_ij - val_ji) >= float32(1e-8):
                is_symmetric = False

            # NaN check (NaN != NaN)
            if val_ij != val_ij:
                has_nan = True

            # Infinity check (Codon float32 max)
            if val_ij >= MAX_FLOAT:
                has_inf = True

            # Negative check
            if val_ij < float32(0.0):
                has_negative = True

            if not is_symmetric or has_nan or has_inf or has_negative:
                break
        if not is_symmetric or has_nan or has_inf or has_negative:
            break

    if not is_symmetric:
        raise ValueError("Distance matrix must be symmetric.")
    if has_nan:
        raise ValueError("Distance matrix contains NaN values.")
    if has_inf:
        raise ValueError("Distance matrix contains infinity values.")
    if has_negative:
        raise ValueError("Distances must be non-negative.")
    if rows < 4:
        raise ValueError("At least four nodes are required.")


    # matrix initialization 
    n: int = D.shape[0]

    # Prepare node containers and helper arrays
    nodes_list: List[TreeNode] = []
    # 0 = active, 1 = clustered
    is_clustered = full(n, 0, dtype=uint8)
    distances_mat: ndarray[float32, 2] = D
    divergence: ndarray[float32, 1] = zeros(n, dtype=float32)
    corr_distances: ndarray[float32, 2] = zeros((n, n), dtype=float32)

    # Build initial leaf nodes
    for i in range(n):
        node: TreeNode = TreeNode(index=i)
        nodes_list.append(node)

    nodes = array(nodes_list)


    # NJ main loop
    while True:

        n_rem_nodes: int = 0
        for i in range(n):
            if is_clustered[i] == uint8(0):
                n_rem_nodes += 1
        if n_rem_nodes <= 1:
            break

        # (1) Compute divergence = sum of distances for each active node
        for i in range(n):
            if is_clustered[i] != uint8(0):
                continue
            dist_sum: float32 = float32(0.0)
            for k in range(n):
                if is_clustered[k] != uint8(0):
                    continue
                dist_sum += distances_mat[i, k]
            divergence[i] = dist_sum



        # (2) Compute corrected distance matrix (Q-matrix)
        for i in range(n):
            if is_clustered[i] != uint8(0):
                continue
            for j in range(i):
                if is_clustered[j] != uint8(0):
                    continue
                corr_distances[i, j] = (
                    float32(n_rem_nodes - 2) * distances_mat[i, j]
                    - divergence[i]
                    - divergence[j]
                )

        # (3) Find minimum corrected distance pair
        dist_min: float32 = MAX_FLOAT
        i_min: int = -1
        j_min: int = -1
        for i in range(n):
            if is_clustered[i] != uint8(0):
                continue
            for j in range(i):
                if is_clustered[j] != uint8(0):
                    continue
                d: float32 = corr_distances[i, j]
                if d < dist_min:
                    dist_min = d
                    i_min = i
                    j_min = j

        # Stop if no valid pair found
        if i_min == -1 or j_min == -1:
            break

        # (4) Compute new limb lengths (branch distances)
        inv_term: float32 = float32(0.0)
        if n_rem_nodes > 2:
            inv_term = float32(1.0) / float32(n_rem_nodes - 2)

        node_dist_i: float = float(
            float32(0.5) * (
                distances_mat[i_min, j_min]
                + inv_term * (divergence[i_min] - divergence[j_min])
            )
        )
        node_dist_j: float = float(
            float32(0.5) * (
                distances_mat[i_min, j_min]
                + inv_term * (divergence[j_min] - divergence[i_min])
            )
        )

        if n_rem_nodes > 3:
            # Create a node with two children
            # (Assumes TreeNode accepts (children, branch_lengths) as you used)
            # mypy/Codon note: nodes[j_min] is not None here
            node_i = nodes[i_min]
            node_j = nodes[j_min]
            # Safety hints for Codon’s type inference (both must be non-None at this point)
            if node_i is None or node_j is None:
                raise ValueError("Internal error: unexpected None node")

            children_list: List[TreeNode] = [node_i, node_j]
            distances_list: List[float] = [node_dist_i, node_dist_j]

            # Replace node_i with new internal node
            # problem 
            nodes[i_min] = TreeNode(children_list, distances_list)
            nodes[j_min] = None
            is_clustered[j_min] = uint8(1)



        else:
            # Last join: build root with the remaining third node k
            node_i = nodes[i_min]
            node_j = nodes[j_min]
            is_clustered[i_min] = uint8(1)
            is_clustered[j_min] = uint8(1)

            k: int = -1
            for t in range(n):
                if is_clustered[t] == uint8(0):
                    k = t
                    break
            if k == -1:
                raise ValueError("Could not find third node for root")

            node_k = nodes[k]
            if node_i is None or node_j is None or node_k is None:
                raise ValueError("Node is None when creating root")

            node_dist_k: float = float(
                float32(0.5)
                * (
                    distances_mat[i_min, k]
                    + distances_mat[j_min, k]
                    - distances_mat[i_min, j_min]
                )
            )
            
            root_children: List[TreeNode] = [node_i, node_j, node_k]
            root_dists: List[float] = [node_dist_i, node_dist_j, node_dist_k]
            root: TreeNode = TreeNode(root_children, root_dists)

            return Tree(root)


        # 5) Update distances for the new cluster at i_min
        for k in range(n):
            if is_clustered[k] ==uint8(0) and k != i_min:
                newd: float32 = float32(
                    float32(0.5)
                    * (
                        distances_mat[i_min, k]
                        + distances_mat[j_min, k]
                        - distances_mat[i_min, j_min]
                    )
                )
                distances_mat[i_min, k] = newd
                distances_mat[k, i_min] = newd

    # If we ever break without returning (shouldn’t happen in a valid NJ run),
    # just create a trivial tree of the first non-None node.
    # This is defensive; in practice, the return happens in the loop.
    for t in range(n):
        if is_clustered[t] == uint8(0) and nodes[t] is not None:
            return Tree(nodes[t])
    raise ValueError("Failed to build NJ tree.")
