

import numpy as np


def inner_product(a, b):
    res = 0
    for i in range(len(a)):
        res += a[i] * b[i]
    return(res)

def scalar_product(a, k):
    res = a.copy()
    for i in range(len(a)):
        res[i] *= k
    return(res)

def scalar_sum(a, k):
    res = []
    for i in range(len(a)):
        res.append(a[i] + k)
    return(res)

def vector_sum(a, b):
    res = []
    for i in range(len(a)):
        res.append(a[i] + b[i])
    return(res)

def magnitude(v):
    return(np.sqrt(inner_product(v, v)))


def vector_product(a, b):
    return([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])

def proj(u, v):
    return(scalar_product(u, inner_product(v, u)/inner_product(u, u)))


# Matrix functions

def matrix_operator(M, v):
    return(np.matmul(M, v))
    # NxN times Nx1
    res = np.zeros(len(v))
    for i in range(len(M)):
        for j in range(len(M)):
            res[i] += M[i][j] * v[j]
    return(res)

def rotation_matrix(n, c_t, s_t):
    
    # normalised axis vector, cos of theta, sin of theta
    # note: this is RIGHT-HAND rotation
    
    n_1 = n[0]
    n_2 = n[1]
    n_3 = n[2]
    
    rot_11 = c_t + n_1 * n_1 * (1.0 - c_t)
    rot_21 = n_1 * n_2 * (1.0 - c_t) + n_3 * s_t
    rot_31 = n_1 * n_3 * (1.0 - c_t) - n_2 * s_t
    
    rot_12 = n_1 * n_2 * (1.0 - c_t) - n_3 * s_t
    rot_22 = c_t + n_2 * n_2 * (1.0 - c_t)
    rot_32 = n_2 * n_3 * (1.0 - c_t) + n_1 * s_t
    
    rot_13 = n_1 * n_3 * (1.0 - c_t) + n_2 * s_t
    rot_23 = n_2 * n_3 * (1.0 - c_t) - n_1 * s_t
    rot_33 = c_t + n_3 * n_3 * (1.0 - c_t)
    
    result = [[rot_11, rot_21, rot_31], [rot_12, rot_22, rot_32], [rot_13, rot_23, rot_33]]
    return(result)


# Advanced decomposition functions

def gram_schmidt(vectors):
    new_vectors = []
    for i in range(len(vectors)):
        cur_u = vectors[i].copy()
        for j in range(i):
            cur_u = vector_sum(cur_u, scalar_product(proj(new_vectors[j], vectors[i]), -1))
        new_vectors.append(list(cur_u))
    return(new_vectors)

def perpendicularize(vector, unperp_vectors):
    perp_vectors = gram_schmidt(unperp_vectors)
    res = vector.copy()
    for i in range(len(perp_vectors)):
        res = vector_sum(res, scalar_product(proj(perp_vectors[i], res), -1))
    print("test")
    for i in range(len(perp_vectors)):
        print(f"res dot unperp[{i+1}] =", inner_product(res, unperp_vectors[i]))
    return(res)


def decompose_last_column_of_singular_matrix(s_m):
    # takes a square matrix M[row][item], where the last column can be written as a linear combination of the previous ones
    # and finds the coefficients in this linear combination
    
    N = len(s_m)
    # note: we can ignore the last row, since it doesn't convey any information (we reduce the matrix to (N-1)x(N-1))
    r_m = np.zeros((N-1, N-1))
    dep_vec = np.zeros(N-1)
    for i in range(N-1):
        for j in range(N-1):
            r_m[i][j] = s_m[i][j]
        dep_vec[i] = s_m[i][N-1]
    r_m_inv = np.linalg.inv(r_m)
    return(matrix_operator(r_m_inv, dep_vec))






