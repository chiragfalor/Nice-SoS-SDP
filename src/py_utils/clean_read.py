import sympy as sp

EPS = 1e-6

def is_int(number, eps=EPS):
    return abs(number - round(number)) < eps

def clean_qtr_int(real, eps=EPS):
    if is_int(4*real, eps):
        return sp.Rational(round(4*real), 4)
    else:
        return real
    

def clean_complex(number, eps=EPS):
    if abs(number) < eps:
        return 0
    elif abs(sp.im(number)) < eps:
        return clean_qtr_int(sp.re(number), eps)
    else:
        re_part = clean_qtr_int(sp.re(number), eps)
        im_part = clean_qtr_int(sp.im(number), eps)
        return re_part + im_part*sp.I
    
    
def round_coeffs(expr, eps=EPS):
    expr = expr.expand()
    rounded_expr = 0
    for term in expr.as_ordered_terms():
        coeff, var = term.as_coeff_Mul()
        rounded_coeff = clean_complex(coeff, eps)
        rounded_expr += rounded_coeff * var
    return rounded_expr.expand()

def clean_matrix(matrix, mult_factor=1, eps=EPS):
    round_c = lambda expr: round_coeffs(expr, eps)
    return (matrix * mult_factor).applyfunc(round_c) / mult_factor

def read_julia_mat(filename, eps=EPS):
    with open(filename, "r") as file:
        # Read each line, replace "im" with "j", and strip whitespace
        lines = [line.replace("im", "j").strip() for line in file]

    lines = lines[1:] # remove the header of the matrix

    lines = [line.replace(" ", "").split(",") for line in lines]
    data = [[complex(number) for number in line] for line in lines]
    data = sp.Matrix(data)
    data = clean_matrix(data, eps=eps)
    return data
