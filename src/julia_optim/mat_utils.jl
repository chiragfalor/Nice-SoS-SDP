using Symbolics
@variables A0, A1, A2, B0, B1, B2;

function conj(Bx, degree)
    # Handle base cases first
    if isequal(Bx, B0) || isequal(Bx, B1) || isequal(Bx, B2)
        return Bx^(degree - 1)
    end
    
    # For expressions of form Bx^k, we need to:
    # 1. Extract the base (B0 or B1)
    # 2. Extract the exponent k
    # 3. Return base^(degree - k)
    
    # Try to match pattern of Bx^k
    if isa(Bx, Expr) && Bx.head == :call && Bx.args[1] == :^ &&
       (isequal(Bx.args[2], B0) || isequal(Bx.args[2], B1) || isequal(Bx.args[2], B2))
        base = Bx.args[2]
        k = Bx.args[3]
        return base^(degree - k)
    end
    
    return nothing
end

function get_id(Ax, By, basis, degree)
    Ax_id = findfirst(isequal(Ax), basis)
    By_id = findfirst(isequal(conj(By, degree)), basis)
    if !any(isnothing, [Ax_id, By_id])
        return Ax_id, By_id
    else
        AxBy_id = findfirst(isequal(Ax*By), basis)
        I_id = findfirst(isequal(1), basis)
        if !any(isnothing, [AxBy_id, I_id])
            return AxBy_id, I_id
        else
            return nothing
        end
    end
end


function find_A0_A1_indices(basis)
    # Find indices containing A0
    A0_indices = findall(term -> occursin("A0", string(term)), basis)
    
    # Find indices containing A1
    A1_indices = findall(term -> occursin("A1", string(term)), basis)
    
    # Return as list of [i,j] pairs
    return [[i, j] for i in A0_indices for j in A1_indices]
end