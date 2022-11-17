
# ---------------------------------------------------------------------------------
# A python module that defines a function mp_normal_modes which takes three inputs:
#     1. N: The number of segments in the mp
#     2. m (optional): a list of strings of length N with the masses of the joints. Defaults to (m_1, m_2... m_N)
#     3. l (optional): a list of strings of length N with the lengths of the segments. Defaults to (l_1, l_2... l_N)
# ---------------------------------------------------------------------------------
# Created by Michal Horansky, 2022
# ---------------------------------------------------------------------------------


# Libraries
import numpy as np
import sympy as sp
from sympy.assumptions import global_assumptions
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt

# Preset matrices generation

def zero_matrix(d_y, d_x=0):
    res = []
    if d_x == 0:
        d_x = d_y
    for i in range(d_y):
        res.append([])
        for j in range(d_x):
            res[i].append(0)
    return(res)

def identity_matrix(d):
    res = []
    for i in range(d):
        res.append([])
        for j in range(d):
            if i == j:
                res[i].append(1.0)
            else:
                res[i].append(0.0)
    return(res)

# Meta functions

def sanitize_input_list(inp, ch, N):
    if inp == 0:
        inp = []
        for i in range(N):
            inp.append(sp.symbols(f"{ch}_{i+1}"))
    else:
        for i in range(N):
            inp[i] = sp.parse_expr(inp[i])
    return(inp)

def assume_real_positive(symbols_list):
    for my_symbol in symbols_list:
        global_assumptions.add(sp.Q.real(my_symbol))
        global_assumptions.add(sp.Q.positive(my_symbol))

def unwrap_eigenvectors(raw_output, refine_expr = True):
    # returns a tuple (eigenvectors, eigenvalues)
    result = []
    def unwrap_local_output(output_matrix):
        res_eigenvectors = []
        res_eigenvalues = []
        for i in range(len(output_matrix)):
            if refine_expr:
                res_eigenvalues.append(sp.refine(output_matrix[i][0]))
                cur_ev = []
                for j in range(len(output_matrix[i][1])):
                    cur_ev.append(sp.refine(output_matrix[i][1][j]))
                res_eigenvectors.append(cur_ev)
                #res_eigenvectors.append([sp.refine(output_matrix[i][1][0]), sp.refine(output_matrix[i][1][1]), sp.refine(output_matrix[i][1][2])])
            else:
                res_eigenvalues.append(output_matrix[i][0])
                res_eigenvectors.append(output_matrix[i][1])
        return((res_eigenvectors, res_eigenvalues))
    for i in range(len(raw_output)):
        eigenvalue, multiplicity, eigenspace = raw_output[i]
        for eigenvector in eigenspace:
            result.append([eigenvalue, sum(eigenvector.tolist(), [])])
    return(unwrap_local_output(result))

# Physics

def mu(m, i, N):
    res = 0
    for j in range(i, N):
        res += m[j]
    return(res)

def create_k_matrix(N, m, l, g):
    # create an empty matrix
    k = zero_matrix(N)
    # populate the edge cases i=j=N and j=N-1,i=N
    k[N-1][N-1] = g / l[N-1]
    if N > 1:
        k[N-2][N-1] = -(g/l[N-1]) * sp.sqrt(m[N-1]/m[N-2])
        k[N-1][N-2] = k[N-2][N-1]
        # populate the main diagonal
        for i in range(N-1):
            k[i][i] = g * (1 / l[i] + (mu(m, i+1, N) / m[i]) * (1/l[i]+1/l[i+1]))
        if N > 2:
            # populate the off-one diagonal
            for i in range(N-2):
                k[i+1][i] = -g * mu(m, i+1, N) / (l[i+1]*sp.sqrt(m[i]*m[i+1]))
                k[i][i+1] = k[i+1][i]
    return(k)


# -------------------------------------------------------------------------------------
# ----------------------------- OUTPUT FUNCTIONS --------------------------------------
# -------------------------------------------------------------------------------------

def print_matrix(matrix, margin=2):
    separator = ' | '
    d_y = len(matrix)
    d_x = len(matrix[0])
    max_width_in_column = [0]*d_x
    offset_matrix = zero_matrix(d_y, d_x)
    for i in range(d_y):
        for j in range(d_x):
            cur_str_len = len(str(matrix[i][j]))
            offset_matrix[i][j] = -cur_str_len
            if cur_str_len > max_width_in_column[j]:
                max_width_in_column[j] = cur_str_len
    breakline_len = sum(max_width_in_column) + (d_x-1)*len(separator)
    breakline = '-'*breakline_len
    for i in range(d_y):
        cur_line = ''
        for j in range(d_x-1):
            left_offset = int(np.floor((offset_matrix[i][j]+max_width_in_column[j])/2))
            cur_line += ' '*left_offset + str(matrix[i][j]) + ' '*((offset_matrix[i][j]+max_width_in_column[j]) - left_offset) + separator
        left_offset = int(np.floor((offset_matrix[i][d_x-1]+max_width_in_column[d_x-1])/2))
        cur_line += ' '*left_offset + str(matrix[i][d_x-1]) + ' '*((offset_matrix[i][d_x-1]+max_width_in_column[d_x-1]) - left_offset)
        print(' '*margin + cur_line)
        if i < d_y-1:
            print(' '*margin + breakline)

def print_vectors(vectors, descriptor = ('', ' = '), margin=2, alignment='left', right_strings = ''):
    d_y = len(vectors)
    d_x = len(vectors[0])
    desc_1, desc_2 = descriptor
    separator = ', '
    line_end = ')'
    if type(right_strings) == str:
        right_strings = [right_strings] * d_y
    max_width_in_column = [0]*d_x
    offset_matrix = zero_matrix(d_y, d_x)
    for i in range(d_y):
        for j in range(d_x):
            cur_str_len = len(str(vectors[i][j]))
            offset_matrix[i][j] = -cur_str_len
            if cur_str_len > max_width_in_column[j]:
                max_width_in_column[j] = cur_str_len
    for i in range(d_y):
        cur_line = desc_1 + str(i+1) + desc_2 + '('
        for j in range(d_x-1):
            if alignment == 'center':
                left_offset = int(np.floor((offset_matrix[i][j]+max_width_in_column[j])/2))
                cur_line += ' '*left_offset + str(vectors[i][j]) + ' '*((offset_matrix[i][j]+max_width_in_column[j]) - left_offset) + separator
            elif alignment == 'left':
                cur_line += str(vectors[i][j]) + ' '*(offset_matrix[i][j]+max_width_in_column[j]) + separator
        if alignment == 'center':
            left_offset = int(np.floor((offset_matrix[i][d_x-1]+max_width_in_column[d_x-1])/2))
            cur_line += ' '*left_offset + str(vectors[i][d_x-1]) + ' '*((offset_matrix[i][d_x-1]+max_width_in_column[d_x-1]) - left_offset) + line_end + right_strings[i]
        elif alignment == 'left':
            cur_line += str(vectors[i][d_x-1]) + ' '*(offset_matrix[i][d_x-1]+max_width_in_column[d_x-1]) + line_end + right_strings[i]
        print(' '*margin + cur_line)

# ---------------------------------------------------------------------------------
# The main function
# ---------------------------------------------------------------------------------


def mp_normal_modes(N, m=0, l=0):
    
    # First, sanitize inputs and create a gravitational acceleration symbol
    
    m = sanitize_input_list(m, "m", N)
    l = sanitize_input_list(l, "l", N)
    
    g = sp.symbols("g")
    
    # Add assumptions of positivity and realness
    assume_real_positive(m)
    assume_real_positive(l)
    assume_real_positive([g])
    
    # Now, create the k matrix
    
    k = create_k_matrix(N, m, l, g)
    print_matrix(k)
    
    k_M = Matrix(k.copy())
    k_M_eigenvectors_raw = k_M.eigenvects(simplify=True)
    modes, f_sq = unwrap_eigenvectors(k_M_eigenvectors_raw, refine_expr = True)
    res_freq = []
    for i in range(len(f_sq)):
        res_freq.append(sp.simplify(sp.sqrt(f_sq[i])))
    print("Resonant frequencies:")
    for i in range(len(res_freq)):
        print(f"  omega_{i+1} = {res_freq[i]}")
    

N = int(input("Number of segments in the multipendulum: "))
print("all segments equal:")
mp_normal_modes(N, ['m']*N, ['l']*N)

print("general case:")
mp_normal_modes(N)


