
setprecision(256)
include("mat_utils.jl")
include("../shared_data/input.jl")

function obj_fn(Re_M, Im_M)
    return get_moment_matrix(Re_M, Im_M, basis)
end

G = get_game_matrix(basis)
cross_constraints = find_A0_A1_indices(basis)

include("SoS_sdp.jl")

M_mat, game_value = game_SoS(G, eq_constraints, cross_constraints)
moment_mat, game_value = moment_SoS(length(basis), eq_constraints, obj_fn)

# Save the output to a shared folder
shared_folder_path = joinpath(@__DIR__, "../shared_data")
isdir(shared_folder_path) || mkdir(shared_folder_path)

using CSV, DataFrames

df = DataFrame(moment_mat, :auto)
output_file_path = joinpath(shared_folder_path, "moment.csv")
CSV.write(output_file_path, df, float_format=:compact)

df = DataFrame(M_mat, :auto)
output_file_path = joinpath(shared_folder_path, "M_mat.csv")
CSV.write(output_file_path, df, float_format=:compact)

moment_mat, M_mat, game_value