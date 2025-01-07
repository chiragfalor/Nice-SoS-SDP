
setprecision(256)
include("temp/input.jl")
include("B3_utils.jl")

function obj_fn(Re_M, Im_M)
    return get_B3_objective(Re_M, Im_M, basis)
end

G = get_B3_matrix(basis)
cross_constraints = find_A0_A1_indices(basis)


include("utils.jl")


M_mat, game_value = game_SoS(G, eq_constraints, cross_constraints)
moment_mat, game_value = moment_SoS(length(basis), eq_constraints, obj_fn)

using CSV, DataFrames
df = DataFrame(moment_mat, :auto)
CSV.write("temp/moment.csv", df, float_format=:compact)

df = DataFrame(M_mat, :auto)
CSV.write("temp/M_mat.csv", df, float_format=:compact)


moment_mat, M_mat, game_value