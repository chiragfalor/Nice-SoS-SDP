from typing import List, Union
import sympy as sp
import sys
sys.path.append("..")
from game_simps import simplify

class SoS_term:
    def __init__(self, coeff, u, b):
        self.coeff = coeff
        self.u = u.expand()
        self.b = b
        self.expr = (self.u * self.b)[0]

    def expand(self):
        return self.expr.adjoint() * self.coeff * self.expr
    
    def _repr_latex_(self):
        return f"${sp.latex(self.coeff)} \\cdot \\left|{sp.latex(self.expr)}\\right|^2$"
    
    def __floordiv__(self, other: Union[sp.Rational, int]):
        new_coeff = self.coeff / other
        u = self.u * sp.sqrt(other)
        return SoS_term(new_coeff, u, self.b)

    
    def __pow__(self, other: Union[sp.Rational, int]):
        return SoS_term(self.coeff * other, self.u / sp.sqrt(other), self.b)
    
    def subs(self, subs):
        return SoS_term(self.coeff.subs(subs), self.u.subs(subs), self.b.subs(subs))
    

    def __mul__(self, other):
        assert (other * other.adjoint()).simplify() == 1
        return SoS_term(self.coeff, (self.u * other).expand(), self.b)
    

    def normalize(self):
        try:
            const_id = next(i for i, x in enumerate(self.b) if x == 1)
        except StopIteration:
            return self
        
        c = self.u[const_id]
        c_norm = (c * c.adjoint()).simplify()
        u = (self.u * c.adjoint() / c_norm).expand()
        coeff = (self.coeff * c_norm).simplify()
        return SoS_term(coeff, u, self.b)
    
    

class SoS_expr:
    def __init__(self, terms: List[SoS_term], relations: dict={}):
        self.terms = terms
        self.relations = relations
        self.simplification = simplify(sum(term.expand() for term in self.terms), self.relations)

    @classmethod
    def from_Ucb(cls, U, c, b, relations={}):
        return cls(
            [
                SoS_term(
                    c[i, i],
                    U[i, :],
                    b
                )
                for i in range(U.shape[0])
            ], 
            relations=relations
        )


    def simplify(self):
        return self.simplification
    
    def subs(self, subs):
        return SoS_expr([term.subs(subs) for term in self.terms], self.relations)
    
    def __getitem__(self, key):
        return self.terms[key]

    def __setitem__(self, key, value):
        self.terms[key] = value
        self.simplification = simplify(sum(term.expand() for term in self.terms), self.relations)
    
    def _repr_latex_(self):
        mat = sp.Matrix([[term.coeff, term.expr] for term in self.terms]).T
        return f"${sp.latex(mat)}$"
    

    def output_latex(self):
        terms_str = "\\\\ \n".join([f"\\lambda_{i+1} = {sp.latex(term.coeff)} &,& S_{i+1} = {sp.latex(term.expr)}" for i, term in enumerate(self.terms)])
        game_poly_str = f"""
\\begin{{eqnarray}}
G = {sp.latex(self.simplification)}\\\\
G = \\sum_{{i=1}}^{{{len(self.terms)}}} \\lambda_i {{S_i^\\dagger S_i}}
\\end{{eqnarray}}
where,
\\begin{{eqnarray}}
{terms_str}
\\end{{eqnarray}}
"""
        return game_poly_str
    
    def assemble(self) -> 'SoS_expr':
        """
        Expand each term's basis to include all basis elements from all terms,
        maintaining consistency by adding zeros for missing elements.
        
        Returns:
            SoS_expr: New expression with standardized basis across all terms
        """
        # Get all unique basis elements while preserving order
        all_bases = [term.b for term in self.terms]
        full_basis = merge_bases(all_bases)
        
        # Expand each term to use the full basis
        new_terms = [expand_term_to_basis(term, full_basis) for term in self.terms]
        
        return SoS_expr(new_terms, self.relations)
    
    def get_Ucb(self):
        sos = self.assemble()
        U = sp.Matrix([t.u for t in sos.terms])
        c = sp.diag(*[t.coeff for t in sos.terms])
        b = sos.assemble()[0].b
        return U, c, b
    

def merge_bases(bases):
    """
    Merge multiple bases into a single comprehensive basis while preserving order
    and removing duplicates.
    
    Args:
        bases (List[List]): List of bases to merge
        
    Returns:
        List: Combined basis with unique elements preserving order
    """
    seen = set()
    merged = []
    for basis in bases:
        for element in list(basis):
            if element not in seen:
                seen.add(element)
                merged.append(element)
    return sp.Matrix(merged)

def expand_term_to_basis(term: SoS_term, full_basis) -> SoS_term:
    """
    Expand a term's u vector to match the full basis by inserting zeros
    for missing basis elements.
    
    Args:
        term (SoS_term): Original term to expand
        full_basis (List): Complete basis to expand to
        
    Returns:
        SoS_term: New term with expanded u vector
    """
    # Create mapping from basis elements to their positions in the original basis
    original_basis_map = {b: i for i, b in enumerate(term.b)}
    
    # Create new u vector with zeros for the full basis
    new_u = sp.zeros(1, len(full_basis))
    
    # Fill in values from original u where they exist in the mapping
    for i, basis_element in enumerate(full_basis):
        # basis_str = basis_element
        if basis_element in original_basis_map:
            new_u[i] = term.u[original_basis_map[basis_element]]
    
    return SoS_term(term.coeff, new_u, full_basis)