
setprecision(256)
include("mat_utils.jl")
include("../shared_data/input.jl")

G = get_game_matrix(basis)
cross_constraints = find_A0_A1_indices(basis)

include("SoS_sdp.jl")
M_mat, game_value = game_SoS(G, eq_constraints, cross_constraints)

using CSV, DataFrames
df = DataFrame(M_mat, :auto)

# Save the output to a shared folder
shared_folder_path = joinpath(@__DIR__, "../shared_data")
isdir(shared_folder_path) || mkdir(shared_folder_path)

output_file_path = joinpath(shared_folder_path, "output.csv")
CSV.write(output_file_path, df, float_format=:compact)

M_mat, game_value