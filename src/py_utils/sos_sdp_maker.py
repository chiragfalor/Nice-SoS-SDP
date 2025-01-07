import sympy as sp

from .utils import simplify_mat

# def get_connected_components(nodes, is_connected):
#     nodes = set(nodes)
#     connected = []
#     while nodes:
#         node = nodes.pop()
#         group = {node}
#         new_nodes = {node}
#         while new_nodes:
#             new_nodes = {n for n in nodes if any(is_connected(n, n2) for n2 in new_nodes)}
#             group.update(new_nodes)
#             nodes -= new_nodes
#         connected.append(group)
#     return connected


def get_equivalence_classes(nodes, is_equivalent):
    nodes = set(nodes)
    equivalence = []
    while nodes:
        node = nodes.pop()
        for eq in equivalence:
            eq_node = next(iter(eq))
            if is_equivalent(node, eq_node):
                eq.add(node)
                break
        else:
            equivalence.append({node})
    return equivalence

# def eq_cons(equivalence_idx, mat):
#     equivalence_idx = set(equivalence_idx)
#     assert len(equivalence_idx) > 1
#     first_idx = equivalence_idx.pop()
#     cmats = []
#     for idx in equivalence_idx:
#         cmat = sp.zeros(mat.shape[0], mat.shape[1])
#         cmat[first_idx] = 1
#         cmat[idx] = -1
#         cmats.append(cmat)
#     return cmats

def get_basis_matrix(basis, relations):
    return simplify_mat(sp.Matrix([[O2.adjoint() * O1 for O1 in basis] for O2 in basis]), relations)

def get_constraint_tuples(basis, relations):
    
    mat = get_basis_matrix(basis, relations)
    indices = [(i, j) for i in range(mat.shape[0]) for j in range(mat.shape[1])]

    def is_connected(idx1, idx2):
        return mat[idx1] == mat[idx2]

    connected = get_equivalence_classes(indices, is_connected)

    return connected


def get_julia_game_function(poly):

    # Construct the header and initialization part of the Julia function
    julia_code = f"function get_game_matrix(basis)\n"
    julia_code += "    w = exp(2*pi*im/3)\n"
    julia_code += "    basis_len = length(basis)\n"
    julia_code += "    matrix = zeros(ComplexF64, basis_len, basis_len)\n\n"

    # Construct the matrix assignment lines
    for coeff, variables in poly:
        # Convert sympy expressions to Julia-friendly strings
        variables_str = ", ".join([str(var).replace('**', '^') for var in variables])
        coeff_str = str(coeff).replace('**', '^')
        julia_code += f"    matrix[get_id({variables_str}, basis, op_degree)...] = {coeff_str}\n"

    # Finalize the function
    julia_code += "\n    G = -transpose(matrix)\n"
    julia_code += "\n    return G\n"
    julia_code += "end\n"

    return julia_code


def get_julia_moment_function(poly):
    
    # Construct the header and initialization part of the Julia function
    julia_code = f"function get_moment_matrix(basis)\n"
    julia_code += "    w = exp(2*pi*im/3)\n"
    julia_code += "    game_poly = [\n"

    # Construct the matrix assignment lines
    for coeff, variables in poly:
        # Convert sympy expressions to Julia-friendly strings
        variables_str = ", ".join([str(var).replace('**', '^') for var in variables])
        coeff_str = str(coeff).replace('**', '^')
        julia_code += f"       [{coeff_str}, get_id({variables_str}, basis)],\n"

    # Finalize the function
    julia_code += "    ]\n"
    julia_code += "\n    G = sum(real(c) * Re_M[id[1], id[2]] + imag(c) * Im_M[id[1], id[2]] for (c, id) in game_poly)\n"
    julia_code += "\n    return G\n"
    julia_code += "end\n"

    return julia_code

def output_julia(eq_constraints, basis, poly, op_degree, file_name='shared_data/input.jl'):
    '''
    Output the constraints in Julia code as a list of list of lists.
    '''
    eq_constraints = [list(map(list, eq_constraint)) for eq_constraint in eq_constraints]
    # adjust the indices to 1-based
    for eq_constraint in eq_constraints:
        for id in eq_constraint:
            id[0] += 1
            id[1] += 1

    # remove the identity eq constraints
    eq_constraints = [eqc for eqc in eq_constraints if eqc[0][0] != eqc[0][1]]

    with open(file_name, 'w') as f:
        f.write(f"basis = {str(list(basis)).replace('**', '^')}\n")
        f.write(f"basis_len = length(basis)\n")
        f.write(f"op_degree = {op_degree}\n\n")
        f.write(get_julia_game_function(poly))
        f.write("\n\n")
        f.write(get_julia_moment_function(poly))
        f.write("\n\n")
        f.write(f"eq_constraints = [\n")
        for eqc in eq_constraints:
            f.write(f"{eqc},\n")
        f.write(f"];\n\n")

    return eq_constraints




