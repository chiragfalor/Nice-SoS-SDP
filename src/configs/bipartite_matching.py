import sympy as sp
from src.py_utils.utils import simplify_basis
from src.py_utils.vars import A0, A1, A2, B0, B1, B2, get_relations

degree = 2
relations = get_relations(degree)

basis = [A0, A1, A2, B0, B1, B2, 1]
basis = simplify_basis(basis, relations)

game_poly = [
    (1, (A0, B0)),
    (-1, (A0, B1)),
    (-1, (A0, B2)),
    (-1, (A1, B0)),
    (1, (A1, B1)),
    (-1, (A1, B2)),
    (-1, (A2, B0)),
    (-1, (A2, B1)),
    (1, (A2, B2)),
]

rounding_coeff = 24