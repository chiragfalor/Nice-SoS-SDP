import sympy as sp
from typing import Sequence, Tuple

def tp(l1, l2):
    # tensor product of two lists
    return [a*b for a in l1 for b in l2]

def adjoint(mat):
    return mat.T.applyfunc(lambda x: x.adjoint())

def simplify(expr, relations={}):
    cur_expr = expr
    while True:
        new_expr = cur_expr.subs(relations).expand()
        if new_expr == cur_expr:
            return new_expr
        cur_expr = new_expr

def simplify_mat(mat, relations={}):
    output = sp.zeros(mat.shape[0], mat.shape[1])
    cache = {}
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if mat[i, j] in cache:
                output[i, j] = cache[mat[i, j]]
            else:
                output[i, j] = simplify(mat[i, j], relations)
                cache[mat[i, j]] = output[i, j]
    return output
    return sp.Matrix([[simplify(mat[i, j], relations) for j in range(mat.shape[1])] for i in range(mat.shape[0])])

def simplify_basis(basis, relations):
    '''remove duplicates from basis'''
    basis = sp.Matrix(basis)
    basis = simplify_mat(basis, relations)
    new_basis = []
    for i in range(basis.shape[0]):
        O = basis[i]
        if all(O != j for j in new_basis):
            new_basis.append(O)
    new_basis = sp.Matrix(new_basis)
    return new_basis


def basis_disp(mat, basis):
    basis_dag = sp.Matrix([O.adjoint() for O in basis]).T
    return sp.Matrix([sp.Matrix([0, *basis]).T, sp.Matrix([[basis_dag.T, mat]])]).expand()


