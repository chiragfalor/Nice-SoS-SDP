import sympy as sp
from src.py_utils.utils import simplify_basis
from src.py_utils.vars import A0, A1, B0, B1, w, get_relations

degree = 3
relations = get_relations(degree)

def get_B3_basis():
    basis1 = [A0, A1, B0**2, B1**2, B0*B1, B1*B0, A0**2, A1**2, B0, B1]
    # basis1 += [x.adjoint() for x in basis1]
    basis2 = [A0*B0, A0*B1]
    basis2 += [x.adjoint() for x in basis2]
    basis3 = [A1*B0, A1*B1]
    basis3 += [x.adjoint() for x in basis3]
    basis = basis1 + basis2 + basis3 + [1]
    basis = simplify_basis(basis, relations)
    basis = sp.Matrix(basis)
    return basis

basis = get_B3_basis()

game_poly = [
    (1, (A0, B0)),
    (1, (A0**2, B0**2)),
    (1, (A0, B1)),
    (1, (A0**2, B1**2)),
    (1, (A1, B0)),
    (1, (A1**2, B0**2)),
    (w, (A1, B1)),
    (w**2, (A1**2, B1**2)),
]


rounding_coeff = 24*26