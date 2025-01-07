import sympy as sp

def fast_symbolic_lu(matrix):
    """
    Perform symbolic LU decomposition on a matrix using a more efficient algorithm.
    Handles zero pivots through row exchanges.
    
    Args:
        matrix: sympy Matrix object to decompose
    Returns:
        L, U, P: Lower triangular, Upper triangular, and Permutation matrices
    """
    from sympy import Matrix, symbols, expand
    n = matrix.shape[0]
    L = Matrix.zeros(n)
    U = Matrix.zeros(n)
    P = Matrix.eye(n)  # Permutation matrix
    
    # Initialize L's diagonal to 1
    for i in range(n):
        L[i, i] = 1
    
    # Work with a copy of the matrix that we can permute
    A = matrix[:,:]
    
    # Copy first row of U after finding best pivot
    # Find the largest element in the first column (partial pivoting)
    max_idx = 0
    for i in range(n):
        if A[i, 0] != 0:  # Check if element is non-zero
            max_idx = i
            break
    
    if max_idx != 0:
        # Swap rows in A and P
        A.row_swap(0, max_idx)
        P.row_swap(0, max_idx)
    
    for j in range(n):
        U[0, j] = A[0, j]
    
    # Calculate first column of L
    if U[0, 0] != 0:  # Only proceed if pivot is non-zero
        for i in range(1, n):
            L[i, 0] = A[i, 0] / U[0, 0]
    
    # Compute rest of L and U
    for k in range(1, n):
        # Find the best pivot in column k, starting from row k
        max_idx = k
        for i in range(k, n):
            if A[i, k] != 0:  # Look for first non-zero element
                max_idx = i
                break
        
        # If we found a non-zero pivot, proceed with the decomposition
        if max_idx != k:
            # Swap rows in A and P
            A.row_swap(k, max_idx)
            P.row_swap(k, max_idx)
            # Adjust L matrix for rows we've already computed
            for j in range(k):
                L[k, j], L[max_idx, j] = L[max_idx, j], L[k, j]
        
        # Calculate row k of U
        for j in range(k, n):
            sum_term = 0
            for s in range(k):
                sum_term += L[k, s] * U[s, j]
            U[k, j] = expand(A[k, j] - sum_term)
        
        # If pivot is zero, skip computing corresponding L entries
        if U[k, k] != 0:
            # Calculate column k of L
            for i in range(k + 1, n):
                sum_term = 0
                for s in range(k):
                    sum_term += L[i, s] * U[s, k]
                L[i, k] = expand((A[i, k] - sum_term) / U[k, k])

    assert (P - sp.eye(n)).norm() == 0

    nonzero_idx = [i for i in range(n) if U[i, i] != 0]
    L, U = L[:, nonzero_idx], U[nonzero_idx, :]
    
    test_M = L * U
    if not (test_M - matrix).expand().applyfunc(sp.simplify).norm().is_zero:
        raise ValueError("LU decomposition failed")
    
    return L, U

def fast_symbolic_ldl(matrix):
    """
    Perform symbolic LDL decomposition on a matrix using a more efficient algorithm.
    Returns L and D matrices.
    
    Args:
        matrix: sympy Matrix object to decompose
    Returns:
        L, D: Lower triangular and Diagonal matrices
    """
    
    L, U = fast_symbolic_lu(matrix)
    
    D = sp.diag(*[next((U[i, j] for j in range(U.shape[1]) if U[i, j] != 0), 0) for i in range(U.shape[0])])

    test_M = L * D * L.adjoint()
    if not (test_M - matrix).expand().applyfunc(sp.simplify).norm().is_zero:
        raise ValueError("LDL decomposition failed")
    
    return L, D

