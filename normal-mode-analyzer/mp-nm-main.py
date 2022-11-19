

import numpy as np

k = [[1, 2, 3], [3, 2, 1], [2, 1, 3]]

eigenvalues, eigenvectors = np.linalg.eig(k)
# save the results in a property normal_modes, which is a list where each element is a two-element list,
# first element is the natural frequency and the second element is a list of length N which is the associated eigenvector
v = []
for i in range(len(eigenvalues)):
    v.append([eigenvalues[i], eigenvectors[:,i]])

print(v)
v.sort(key=lambda x : x[0])

print(eigenvectors)
print(v)
