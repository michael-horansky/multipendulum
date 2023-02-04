

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



