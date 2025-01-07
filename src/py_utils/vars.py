import sympy as sp

from typing import Sequence, Tuple, Iterable

A0, A1, A2, B0, B1, B2 = sp.symbols("A0 A1 A2 B0 B1 B2", commutative=False)
w = sp.symbols("w", complex=True)
a, b, c, d = sp.symbols("a b c d", commutative=True, complex=True)
i_to_w = {sp.I: (2*w + 1)/sp.sqrt(3)}
w_to_i = {w: -sp.Rational(1, 2) + sp.I * sp.sqrt(3)/2}
A0c, A1c, B0c, B1c = A0.adjoint(), A1.adjoint(), B0.adjoint(), B1.adjoint()

def get_conj_relations(O, degree, conj_fn=lambda O: O.adjoint()):
    if isinstance(O, Iterable):
        output = {}
        for o in O:
            output.update(get_conj_relations(o, degree, conj_fn))
        return output
    return {
        # conj_fn(O) * O: 1,
        # O * conj_fn(O): 1,
        O**degree: 1,
        conj_fn(O): O**(degree-1),
    }

def get_comm_relations(A, B, conj_fn=lambda O: O.adjoint()):
    if isinstance(A, Sequence):
        output = {}
        for a in A:
            output.update(get_comm_relations(a, B, conj_fn))
        return output
    if isinstance(B, Sequence):
        output = {}
        for b in B:
            output.update(get_comm_relations(A, b, conj_fn))
        return output
    return {
        B*A : A*B,
        conj_fn(B) * A : A * conj_fn(B),
        B * conj_fn(A) : conj_fn(A) * B,
        conj_fn(B) * conj_fn(A) : conj_fn(A) * conj_fn(B),
    }

def get_relations(degree):
    A_set = [Ai**d for Ai in [A0, A1, A2] for d in range(1, degree)]
    B_set = [Bi**d for Bi in [B0, B1, B2] for d in range(1, degree)]
    relations = get_conj_relations(A_set + B_set, degree) | get_comm_relations(A_set, B_set) | {w: -sp.Rational(1, 2) + sp.I*sp.sqrt(3)/2}
    return relations