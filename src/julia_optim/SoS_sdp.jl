using JuMP
using SCS  # or using Mosek if you prefer
using LinearAlgebra
using PrettyTables
using Printf

# Helper functions for complex SDP
function create_hermitian_variable(model, n)
    # Real part is symmetric
    Re = @variable(model, [1:n, 1:n], Symmetric)
    
    # Imaginary part is skew-symmetric
    Im = @variable(model, [1:n, 1:n])
    for i in 1:n
        @constraint(model, Im[i,i] == 0)  # Diagonal elements must be real
        for j in (i+1):n
            @constraint(model, Im[j,i] == -Im[i,j])  # Skew-symmetric
        end
    end
    
    return Re, Im
end

function get_complex_solution(Re_M, Im_M)
    n = size(Re_M, 1)
    complex_sol = zeros(ComplexF64, n, n)
    for i in 1:n
        for j in 1:n
            complex_sol[i,j] = value(Re_M[i,j]) + im*value(Im_M[i,j])
        end
    end
    return complex_sol
end


# primal
function moment_SoS(n, eq_constraints, get_objective)
    eps = 1e-10
    model = Model(optimizer_with_attributes(SCS.Optimizer, "eps_abs" => eps, "eps_rel" => eps, "eps_infeas" => eps))
    #set_optimizer_attribute(model, "verbose", 0)  # Optional: reduce output

    Re_M, Im_M = create_hermitian_variable(model, n)


    M_block = [Re_M -Im_M; Im_M Re_M]
    @constraint(model, M_block in PSDCone())

    # Objective: maximize game objective
    G = get_objective(Re_M, Im_M)
    @objective(model, Max, G)

    constraints_dict = Dict()

    constraints_dict["identity"] = [@constraint(model, Re_M[idx, idx] == 1) for idx in 1:n]

    # Handle equality constraints for complex matrices
    for (i, eq) in enumerate(eq_constraints)
        if length(eq) > 1
            eq_str = ""
            for (a, b) in eq
                eq_str *= " + $a,$b"
            end
            eq_str = eq_str[4:end]

            id1 = eq[1]
            constraints_dict[eq_str * "_real"] = [
                @constraint(model, Re_M[id1[1], id1[2]] == Re_M[idx[1], idx[2]]) for idx in eq[2:end]
            ]
            constraints_dict[eq_str * "_imag"] = [
                @constraint(model, Im_M[id1[1], id1[2]] == Im_M[idx[1], idx[2]]) for idx in eq[2:end]
            ]
        end
    end

    optimize!(model)

    optim_M = get_complex_solution(Re_M, Im_M)

    # Get the optimal objective value
    optimal_objective = objective_value(model)

    return optim_M, optimal_objective
end



# Optimizer
# dual
function game_SoS(G, eq_constraints, cross_constraints)
    eps = 1e-10
    model = Model(optimizer_with_attributes(SCS.Optimizer, "eps_abs" => eps, "eps_rel" => eps, "eps_infeas" => eps))
    #set_optimizer_attribute(model, "verbose", 0)  # Optional: reduce output

    n = size(G)[1]
    Re_M, Im_M = create_hermitian_variable(model, n)


    M_block = [Re_M -Im_M; Im_M Re_M]
    @constraint(model, M_block in PSDCone())

    # Objective: minimize trace (which is real for Hermitian matrices)
    @objective(model, Min, tr(Re_M)) #+ 1e-6*sum(Im_M[i, j]^2 for i in 1:n, j in 1:n))

    constraints_dict = Dict()

    # Handle equality constraints for complex matrices
    for (i, eq) in enumerate(eq_constraints)
        eq_str = ""
        for (a, b) in eq
            eq_str *= " + $a,$b"
        end
        eq_str = eq_str[4:end]
        
        # Split real and imaginary parts of the constraint
        constraints_dict[eq_str * "_real"] = @constraint(model, 
            sum(Re_M[a, b] for (a, b) in eq) == real(sum(G[a, b] for (a, b) in eq)))
        
        constraints_dict[eq_str * "_imag"] = @constraint(model, 
            sum(Im_M[a, b] for (a, b) in eq) == imag(sum(G[a, b] for (a, b) in eq)))
    end

    for (i, eq) in enumerate(cross_constraints)
        cross = basis[eq[1]] * basis[eq[2]]
        
        # Split into real and imaginary parts
        constraints_dict[string(cross) * "_real"] = @constraint(model, 
            Re_M[eq[1], eq[2]] == real(G[eq[1], eq[2]]))
        
        constraints_dict[string(cross) * "_imag"] = @constraint(model, 
            Im_M[eq[1], eq[2]] == imag(G[eq[1], eq[2]]))
    end

    optimize!(model)

    optim_M = get_complex_solution(Re_M, Im_M)

    # Get the optimal objective value
    optimal_objective = objective_value(model)

    return optim_M, optimal_objective
end


function RR(M)
    # Make matrix is Hermitian
    H = Hermitian((M + M') / 2)
    return cholesky(H, check=false).U
end

# Formatting functions
function complex_formatter(x, i, j)
    if !isa(x, Complex)  # if x is not a complex number
        return x
    end
    r, im = round(real(x), digits=3), round(imag(x), digits=3)
    if r == 0 && im == 0
        return "0"
    elseif r == 0
        return @sprintf("%.6fi", im)
    elseif im == 0
        return @sprintf("%.6f", r)
    else
        return @sprintf("%.4f%+.4fi", r, im)
    end
end

function mat_disp(mat, basis, complex_format=true)
    mat = ["" basis'; basis mat]
    if complex_format
        PrettyTables.pretty_table(mat, formatters=complex_formatter, alignment=:c, show_header=false)
    else
        PrettyTables.pretty_table(mat, show_header=false)
    end
end;
