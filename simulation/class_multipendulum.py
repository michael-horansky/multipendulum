
import numpy as np

import time

def tree_access(tree, index_list, make_copy = True):
    if make_copy:
        tree_c = tree.copy()
    else:
        tree_c = tree
    for i in index_list:
        tree_c = tree_c[i]
    return(tree_c)

def minor(square_matrix, my_x, my_y):
    res_matrix = []
    for i in range(len(square_matrix)):
        if i + 1 == my_x:
            continue
        res_matrix.append(square_matrix[i][:my_y - 1] + square_matrix[i][my_y:])
    return(res_matrix)


class dict_tree:
    
    def __init__(self):
        
        self.t = {}
    
    def set_val(self, index_list, my_val):
        cur_node = self.t
        for i in range(len(index_list[:-1])):
            if index_list[i] in cur_node:
                cur_node = cur_node[index_list[i]]
            else:
                cur_node[index_list[i]] = {}
                cur_node = cur_node[index_list[i]]
        cur_node[index_list[-1]] = my_val
    
    def add_val(self, index_list, my_val):
        cur_node = self.t
        for i in range(len(index_list[:-1])):
            if index_list[i] in cur_node:
                cur_node = cur_node[index_list[i]]
            else:
                cur_node[index_list[i]] = {}
                cur_node = cur_node[index_list[i]]
        if index_list[-1] in cur_node:
            cur_node[index_list[-1]] += my_val
        else:
            cur_node[index_list[-1]] = my_val
    
    def get_val(self, index_list):
        cur_node = self.t
        for i in index_list:
            cur_node = cur_node[i]
        return(cur_node)




class physical_system:
    
    def get_acceleration(self):
        # returns q double-dot
        return(0)



class multipendulum(physical_system):
    
    # ---------------- constructors, destructors, descriptors ----------------------
    
    def __init__(self, my_l, my_m, my_g):
        
        self.l = np.array(my_l)
        self.m = np.array(my_m)
        self.m_r = np.sqrt(self.m / 20.0)
        
        self.g = my_g
        
        self.N = len(my_l)
        
        self.theta     = np.zeros(self.N)
        self.theta_dot = np.zeros(self.N)
        
        # optimization functors
        #self.M_functor

    
    # --------------- physical methods ------------------------
    
    def random_state(self, interval = 1):
        self.theta     = np.random.rand(self.N) * interval
        self.theta_dot = np.random.rand(self.N) * interval
    
    def reset_state(self):
        self.theta     = np.zeros(self.N)
        self.theta_dot = np.zeros(self.N)
    
    def set_state(self, new_theta, new_theta_dot):
        self.theta     = np.array(new_theta    )
        self.theta_dot = np.array(new_theta_dot)
    
    def get_mu(self, a, b):
        mu = 0
        for i in range(max(a, b), self.N):
            mu += self.m[i]
        return(mu)
    
    def get_M(self):
        M = []
        for a in range(self.N):
            M.append([])
            for b in range(self.N):
                M[a].append(np.cos( self.theta[a] - self.theta[b] ) * self.get_mu(a, b))
        return(M)
    
    def get_M_coefs(self):
        # like get_M, but the cos terms are stored separately (so that M[a][b] = [mu coef, a, b]), and we exploit the symmetricity of M
        M = []
        for a in range(self.N):
            M.append([])
            for b in range(self.N):
                M[a].append([self.get_mu(a, b), min(a, b), max(a,b)])
        return(M)
    
    def get_det_M(self):
        my_M = self.get_M()
        return(np.linalg.det(my_M))
    
    def get_S(self, external_force):
        S = []
        for a in range(self.N):
            S.append(0.0)
            cur_sum_i = 0.0
            for i in range(a, self.N):
                cur_sum_j = 0.0
                for j in range(i+1):
                    cur_sum_j += self.l[j] * self.theta_dot[j] * self.theta_dot[j] * np.sin(self.theta[a] - self.theta[j])
                cur_sum_i += self.m[i] * (self.g * np.sin(self.theta[a]) + cur_sum_j)
            S[a] -= cur_sum_i
        S[0] += external_force
        return(S)
    
    def get_modified_M(self, n, external_force):
        M = self.get_M()
        S = self.get_S(external_force)
        for a in range(self.N):
            M[a][n] = S[a]
        return(M)
    
    def get_det_modified_M(self, n, external_force):
        my_M = self.get_modified_M(n, external_force)
        return(np.linalg.det(my_M))
    
    def get_acceleration(self, external_force):
        det_M = self.get_det_M()
        acceleration = []
        for n in range(self.N):
            acceleration.append( self.get_det_modified_M(n, external_force) / (self.l[n] * det_M) )
        return(acceleration)
    
    
    # ----------------- optimization methods -----------------
    
    # !!!!!!!!!!!!!!! NOTE: The benchmark test concluded that the polynome preparation method is significantly
    # slower than direct determinant computation. We will not use it as it is right now.
    
    def index_flat(self, i, j):
        # assigns a unique index to every element in the upper triangular exponent matrix. i, j refers to the coefficient of cos(theta_i - theta_j)
        return(int(self.N * i - i * (i+1) / 2 + j))
    def index_square(self, x):
        # inverse func to index_flat
        i,j = 0,0
        cur_x = x
        cur_width = self.N
        while(cur_width <= cur_x):
            i+=1
            cur_x -= cur_width
            cur_width -= 1
        j = cur_x+i
        return(i, j)
    
    def get_C_M_size(self):
        
        self.C_M_size = 0
        def recursive_count(subtree):
            if type(subtree)!=dict:
                self.C_M_size += 1
            else:
                for key, item in subtree.items():
                    recursive_count(item)
        recursive_count(self.C_M.t)
        return(self.C_M_size)
        
    
    def print_C_M(self):
        
        def final_print(index_list, value):
            output_str = "Coef of"
            for x in range(len(index_list)):
                if index_list[x] > 0:
                    i, j = self.index_square(x)
                    output_str += " cos^" + str(index_list[x]) + "(t_" + str(i+1) + "-t_" + str(j+1)+")"
            print(output_str + " = " + str(value))
        
        def recursive_print(subtree, cur_index_list = []):
            if type(subtree)!=dict:
                final_print(cur_index_list, subtree)
            else:
                for key, item in subtree.items():
                    recursive_print(item, cur_index_list + [key])
                
        recursive_print(self.C_M.t)
    
    def get_det_M_opt(self):
        # works the same way as print_C_M but actually works out the thing
        self.det_M_opt_val = 0
        def recursive_eval(subtree, cur_val=1.0, cur_x=0):
            if type(subtree)!=dict:
                self.det_M_opt_val += cur_val * subtree
            else:
                for key, item in subtree.items():
                    i, j = self.index_square(cur_x)
                    if key == 0:
                        cos_factor = 1.0
                    else:
                        #cos_factor = np.power(np.cos(self.theta[i] - self.theta[j]), key)
                        cos_factor = (np.cos(self.theta[i] - self.theta[j])) ** key
                    recursive_eval(item, cur_val * cos_factor, cur_x + 1)
        recursive_eval(self.C_M.t)
        return(self.det_M_opt_val)
            
    
    def set_M_functor(self):
        # a polynomial of order N
        # we have N^2 variables: x_ij=cos(theta_i - theta_j)
        # the matrix is given by M_ab = mu_ab * x_ij
        # we shall use Laplace expansion to find the polynomial coefficients
        # the coefficients are saved in a matrix of N x N x ... x N = N^(2N) in the form C[i_1][j_1][i_2][j_2]...[i_N][j_N]
        # but it's actually more effective to use a tree dictionary, as the number of elements is nCr(2N-1, N)
        
        self.C_M = dict_tree()
        
        # tail-recursive, carries list of exponent (matrix NxN in the form C[i][j])
        
        # the matrix of exponents must be an upper triangular matrix, as we dont want to discern between cos(A-B) and cos(B-A)
        
        
        def polynomial_determinant(square_matrix, exponent_matrix, cur_coef = 1.0):
            
            my_exponent_matrix = exponent_matrix.copy()
            # Matrix is in the form of list of rows: M[x][y]
            # Laplace expansion along the 1st row
            if len(square_matrix) == 1:
                # terminate and store the coefficient
                
                my_exponent_matrix[square_matrix[0][0][1]][square_matrix[0][0][2]] += 1
                cur_coef *= square_matrix[0][0][0]
                
                #return(square_matrix[0][0])
                
                cur_index_list = []
                
                #output_str = "Coef of"
                for i in range(self.N):
                    for j in range(self.N):
                        if i > j:
                            continue
                        cur_val = my_exponent_matrix[i][j]
                        cur_index_list.append(my_exponent_matrix[i][j])
                        #if cur_val > 0:
                        #    output_str += " cos^" + str(cur_val) + "(t_" + str(i+1) + "-t_" + str(j+1)+")"
                #print(output_str + " = " + str(cur_coef))
                self.C_M.add_val(cur_index_list, cur_coef)
                        
            #my_sum = 0.0
            for i in range(len(square_matrix)):
                #my_sum += np.power(-1, i) * square_matrix[i][0] * determinant(minor(square_matrix, i + 1, 1))
                my_exponent_matrix[square_matrix[i][0][1]][square_matrix[i][0][2]] += 1
                polynomial_determinant(minor(square_matrix, i + 1, 1), my_exponent_matrix, cur_coef * np.power(-1, i) * square_matrix[i][0][0])
                my_exponent_matrix[square_matrix[i][0][1]][square_matrix[i][0][2]] -= 1
            #return(my_sum)
        
        #self.M_functor_coefs = polynomial_determinant(self.get_M_coefs, np.zeros( (self.N, self.N) ))
        polynomial_determinant(self.get_M_coefs(), np.zeros( (self.N, self.N) ))

