import sympy as sp
import numpy as np
import pickle 

def get_SoS(basis, op_degree, game_poly, rounding_coeff=1, dump_folder=None, relations=None, cleaning_eps=1e-4):
    

    from src.py_utils.utils import tp, simplify, simplify_basis, simplify_mat, basis_disp
    from src.py_utils.sos_sdp_maker import get_basis_matrix, get_constraint_tuples, output_julia
    from src.py_utils.vars import w, i_to_w, w_to_i, A0, A1, B0, B1

    if relations is None:
        from src.py_utils.vars import get_relations
        relations = get_relations(op_degree)

    gp = sum([coeff * sp.expand(var[0] * var[1]) for coeff, var in game_poly])
    basis_dag = sp.Matrix([O.adjoint() for O in basis]).T

    eqs = get_constraint_tuples(basis, relations)

    # # Julia intermission
    output_julia(eqs, basis, game_poly, op_degree, file_name='src/shared_data/input.jl')
    import os
    os.system("julia src/julia_optim/optim.jl")


    # Cleaning and rounding
    from src.py_utils.clean_read import read_julia_mat, round_coeffs, clean_matrix

    mat = read_julia_mat("src/shared_data/output.csv")
    rc = rounding_coeff
    mat = clean_matrix(mat.subs({sp.I: (2*w + 1)/np.sqrt(3)}), rc, eps=cleaning_eps).subs(w_to_i)
    print(mat, basis, basis_dag)
    G = (basis_dag * mat * basis)[0].simplify()
    G = round_coeffs(simplify(G, relations), cleaning_eps)
    game_value = simplify(G + gp.subs(w_to_i))

    if dump_folder is not None:
        with open(f'{dump_folder}/SoS_mat.pkl', 'wb') as f:
            pickle.dump(mat, f)

    # # Getting into SoS Form

    from src.SoS_utils.lu_utils import fast_symbolic_ldl
    from src.SoS_utils.SoS import SoS_expr

    l, d = fast_symbolic_ldl(mat)
    sos = SoS_expr.from_Ucb(l.adjoint(), d, basis, relations)

    assert simplify(sos.simplification + gp.subs(w_to_i) - game_value) == 0
    if dump_folder is not None:
        with open(f'{dump_folder}/SoS.pkl', 'wb') as f:
            pickle.dump(sos, f)

    return sos, mat, game_value



if __name__=="__main__":

    from src.configs.CHSH import basis, game_poly, rounding_coeff, degree
    sos, mat, game_value = get_SoS(basis, degree, game_poly, rounding_coeff=rounding_coeff, dump_folder=None)

    print("The game value is: ", game_value)
    print("The SoS is: ")
    print(sos)