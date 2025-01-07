basis = [A0, A1, B0, B1, 1]
basis_len = length(basis)
op_degree = 2

function get_game_matrix(basis)
    w = exp(2*pi*im/3)
    basis_len = length(basis)
    matrix = zeros(ComplexF64, basis_len, basis_len)

    matrix[get_id(A0, B0, basis, op_degree)...] = 1
    matrix[get_id(A0, B1, basis, op_degree)...] = 1
    matrix[get_id(A1, B0, basis, op_degree)...] = 1
    matrix[get_id(A1, B1, basis, op_degree)...] = -1

    G = -transpose(matrix)

    return G
end


function get_moment_matrix(Re_M, Im_M, basis)
    w = exp(2*pi*im/3)
    game_poly = [
       [1, get_id(A0, B0, basis, op_degree)],
       [1, get_id(A0, B1, basis, op_degree)],
       [1, get_id(A1, B0, basis, op_degree)],
       [-1, get_id(A1, B1, basis, op_degree)],
    ]

    G = sum(real(c) * Re_M[id[1], id[2]] + imag(c) * Im_M[id[1], id[2]] for (c, id) in game_poly)

    return G
end


eq_constraints = [
[[5, 1], [1, 5]],
[[4, 5], [5, 4]],
[[4, 2], [2, 4]],
[[1, 3], [3, 1]],
[[2, 1]],
[[3, 5], [5, 3]],
[[1, 4], [4, 1]],
[[1, 2]],
[[2, 3], [3, 2]],
[[4, 3]],
[[5, 2], [2, 5]],
[[3, 4]],
];

