
import numpy as np

import time

def list_to_str(my_list, unit_str = '', separation_str = ', '):
    if len(my_list) == 0:
        return('')
    if unit_str != '':
        unit_str = " " + unit_str
    output = ""
    for i in range(len(my_list) - 1):
        output += str(my_list[i]) + unit_str + separation_str
    output += str(my_list[-1]) + unit_str
    return(output)

def elementwise_casting(my_list, my_type):
    res_list = []
    for i in range(len(my_list)):
        res_list.append(my_type(my_list[i]))
    return(res_list)

def empty_list_list(length):
    result = []
    for i in range(length):
        result.append([])
    return(result)

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

def nonzero_sign(x):
    if type(x) == list or type(x) == np.ndarray:
        result = np.zeros(len(x))
        for i in range(len(x)):
            result[i] = nonzero_sign(x[i])
        return(result)
    if x >= 0.0:
        return(1)
    else:
        return(-1)

def base_angle(x):
    if type(x) == list or type(x) == np.ndarray:
        result = np.zeros(len(x))
        for i in range(len(x)):
            result[i] = base_angle(x[i])
        return(result)
    res = x
    while res <= - np.pi:
        res += 2.0 * np.pi
    while res > np.pi:
        res -= 2.0 * np.pi
    return(res)

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

def max_item(my_list):
    if type(my_list) == dict:
        assigned = False
        max_val = max_key = 0
        for key, item in df.items():
            if assigned == False:
                assigned = True
                max_key = key
                max_val = item
            if item > max_val:
                max_val = item
                max_key = key
        return(max_key, max_val)
    if type(my_list) == list or type(my_list) == np.ndarray:
        max_val = my_list[0]
        max_i   = 0
        for i in range(len(my_list)):
            if my_list[i] > max_val:
                max_i = i
                max_val = my_list[i]
        return(max_i, max_val)
    #my_list = list(my


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
        
        self.phi       = np.zeros(self.N)
        
        # optimization functors
        #self.M_functor
    
    
    def parameter_description(self):
        # first line is gravity and number of segments
        # second line is length of segments
        # third line is segment masses
        line1 = "g = %0.2f, N = %i" % (self.g, self.N)
        line2 = "segments:"
        line3 = "masses:  "
        for i in range(self.N):
            line2 += " l_%i = %0.2f," % (i, self.l[i])
            line3 += " m_%i = %0.2f," % (i, self.m[i])
        return([line1, line2[:-1], line3[:-1]])

    
    # --------------- manipulation methods ------------------------
    
    def random_state(self, interval = 1):
        self.theta     = np.random.rand(self.N) * interval
        self.theta_dot = np.random.rand(self.N) * interval
    
    def reset_state(self):
        self.theta     = np.zeros(self.N)
        self.theta_dot = np.zeros(self.N)
    
    def set_state(self, new_theta, new_theta_dot):
        self.theta     = np.array(new_theta    )
        self.theta_dot = np.array(new_theta_dot)
    
    # ---------------------------------------------------------
    # --------------- physical methods ------------------------
    # ---------------------------------------------------------
    
    # ---------- calculating simple properties ----------------
    
    def get_total_energy(self):
        # we ignore the driving force potential, as it is technically not a potential
        U_g = 0.0
        T = 0.0
        for i in range(self.N):
            U_g_ps = 0.0
            T_ps_1 = 0.0
            T_ps_2 = 0.0
            for j in range(0, i+1):
                U_g_ps += self.l[j] * np.cos(self.theta[j])
                T_ps_1 += self.theta_dot[j] * self.l[j] * np.cos(self.theta[j])
                T_ps_2 += self.theta_dot[j] * self.l[j] * np.sin(self.theta[j])
            U_g -= self.m[i] * U_g_ps
            T   += 0.5 * self.m[i] * (T_ps_1 * T_ps_1 + T_ps_2 * T_ps_2)
        U_g *= self.g
        return(T, U_g, T+U_g)

    
    def get_mu(self, a, b = 0):
        mu = 0
        for i in range(max(a, b), self.N):
            mu += self.m[i]
        return(mu)
    
    def get_phi(self):
        self.phi[0] = self.theta[0]
        for i in range(1, self.N):
            self.phi[i] = self.theta[i] - self.theta[i-1]
    
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
    
    def get_angular_momentum(self):
        L = 0.0
        for n in range(self.N):
            A = 0.0
            B = 0.0
            C = 0.0
            D = 0.0
            for i in range(n):
                A += self.l[i] * np.sin(self.theta[i])
                B += self.l[i] * np.cos(self.theta[i])
                C += self.theta_dot[i] * self.l[i] * np.sin(self.theta[i])
                D += self.theta_dot[i] * self.l[i] * np.cos(self.theta[i])
            L += A * C + B * D
        return(L)
    
    # --------------- theoretical normal modes analysis ------------------------
    
    def modes_q_space_into_theta_space(self, modes):
        if type(modes) == list:
            new_modes = []
            for old_mode in modes:
                new_modes.append(self.modes_q_space_into_theta_space(old_mode))
            return(new_modes)
        mode = modes
        new_mode = []
        new_mode.append(mode[0] / (self.l[0] * np.sqrt(self.m[0])))
        for i in range(1, self.N):
            new_mode.append((mode[i] / np.sqrt(self.m[i]) - mode[i-1] / np.sqrt(self.m[i-1])) / self.l[i])
        return(new_mode)
    
    
    def get_normal_modes(self):
        
        # create the k matrix
        k = zero_matrix(self.N)
        # populate the edge cases i=j=N and j=N-1,i=N
        k[self.N-1][self.N-1] = self.g / self.l[self.N-1]
        if self.N > 1:
            k[self.N-2][self.N-1] = -(self.g/self.l[self.N-1]) * np.sqrt(self.m[self.N-1]/self.m[self.N-2])
            k[self.N-1][self.N-2] = k[self.N-2][self.N-1]
            # populate the main diagonal
            for i in range(self.N-1):
                k[i][i] = self.g * (1 / self.l[i] + (self.get_mu(i+1) / self.m[i]) * (1/self.l[i]+1/self.l[i+1]))
            if self.N > 2:
                # populate the off-one diagonal
                for i in range(self.N-2):
                    k[i+1][i] = -self.g * self.get_mu(i+1) / (self.l[i+1]*np.sqrt(self.m[i]*self.m[i+1]))
                    k[i][i+1] = k[i+1][i]
        eigenvalues, eigenvectors = np.linalg.eig(k)
        # save the results in a property normal_modes, which is a list where each element is a two-element list,
        # first element is the natural frequency and the second element is a list of length N which is the associated eigenvector
        cur_normal_modes = []
        for i in range(len(eigenvalues)):
            cur_normal_modes.append([np.sqrt(eigenvalues[i]), eigenvectors[:,i]])
        cur_normal_modes.sort(key=lambda x : x[0])
        # convert normal modes from q space into theta space and save into normal_modes (eigenvectors) and normal_mode_frequencies (associated frequencies)
        self.normal_modes = self.modes_q_space_into_theta_space( [row[1] for row in cur_normal_modes] )
        self.normal_mode_frequencies = [row[0] for row in cur_normal_modes]
    
    def modal_frequency(self, mode):
        
        a = 0.0
        b = 0.0
        for i in range(self.N):
            a += self.l[i] * mode[i] * mode[i] * self.get_mu(i)
            c = 0.0
            for j in range(i+1):
                c += mode[j] * self.l[j]
            b += self.m[i] * c * c
        return(np.sqrt(self.g * a / b))
    
    def print_normal_modes(self):
        
        if self.N == 1:
            line = ""
            for i in range(self.N):
                line += f"v_{i+1} = ({self.normal_modes[i][i]:.3f}), omega_{i+1} = {self.normal_mode_frequencies[i]:.3f}; "
            print(line[:-2])
        else:
            main_i = int(np.ceil(self.N/2.0)-1)
            v_lens = []
            w_lens = []
            for i in range(self.N):
                w_lens.append(len(f"{self.normal_mode_frequencies[i]:.3f}"))
                v_lens.append(len(f"{self.normal_modes[i][0]:.3f}"))
                for j in range(1, self.N):
                    cur_len = len(f"{self.normal_modes[i][j]:.3f}")
                    if cur_len > v_lens[-1]:
                        v_lens[-1] = cur_len
            
            for i in range(self.N):
                if i == 0:
                    cur_left = '/'
                    cur_right = "\\"
                elif i == self.N - 1:
                    cur_left = "\\"
                    cur_right = "/"
                else:
                    cur_left = '|'
                    cur_right = '|'
                line = ""
                for j in range(self.N):
                    cur_v_len = len(f"{self.normal_modes[j][i]:.3f}")
                    #cur_w_len = len(f"{self.normal_mode_frequencies[j]:.3f}")
                    
                    if i == main_i:
                        line += f"v_{j+1} = {cur_left}{' ' * (v_lens[j] - cur_v_len)}{self.normal_modes[j][i]:.3f}{cur_right}, omega_{j+1} = {self.normal_mode_frequencies[j]:.3f}; "
                    else:
                        line += f"      {cur_left}{' ' * (v_lens[j] - cur_v_len)}{self.normal_modes[j][i]:.3f}{cur_right}       {' ' * w_lens[j]}       "
                print(line[:-2])
                
    
    # --------------- numerical integration methods ------------------------
    
    def step(self, dt, cur_external_force_list):
        
        # cur external force list is a list of exteral forces at time t, t+dt/2 and t+dt
        
        old_theta     = self.theta
        old_theta_dot = self.theta_dot
            
            
        k1_theta     = []
        k1_theta_dot = []
        acc_1 = self.get_acceleration(cur_external_force_list[0])
        for i in range(self.N):
            k1_theta.append(     dt * self.theta_dot[i] )
            k1_theta_dot.append( dt * acc_1[i]              )
        self.theta     = old_theta     + np.array(k1_theta)      / 2.0
        self.theta_dot = old_theta_dot + np.array(k1_theta_dot)  / 2.0
            
        k2_theta     = []
        k2_theta_dot = []
        acc_2 = self.get_acceleration(cur_external_force_list[1])
        for i in range(self.N):
            k2_theta.append(     dt * self.theta_dot[i] )
            k2_theta_dot.append( dt * acc_2[i]              )
        self.theta     = old_theta     + np.array(k2_theta)      / 2.0
        self.theta_dot = old_theta_dot + np.array(k2_theta_dot)  / 2.0
            
        k3_theta     = []
        k3_theta_dot = []
        acc_3 = self.get_acceleration(cur_external_force_list[1])
        for i in range(self.N):
            k3_theta.append(     dt * self.theta_dot[i] )
            k3_theta_dot.append( dt * acc_3[i]              )
        self.theta     = old_theta     + np.array(k3_theta)
        self.theta_dot = old_theta_dot + np.array(k3_theta_dot)
            
        k4_theta     = []
        k4_theta_dot = []
        acc_4 = self.get_acceleration(cur_external_force_list[2])
        for i in range(self.N):
            k4_theta.append(     dt * self.theta_dot[i] )
            k4_theta_dot.append( dt * acc_4[i]              )
            
        #print("huehue", (np.array(k1_theta)     + 2.0*np.array(k2_theta)     + 2.0*np.array(k3_theta)     + np.array(k4_theta)    )/6.0)
            
        self.theta     = old_theta     + (np.array(k1_theta)     + 2.0*np.array(k2_theta)     + 2.0*np.array(k3_theta)     + np.array(k4_theta)    )/6.0
        self.theta_dot = old_theta_dot + (np.array(k1_theta_dot) + 2.0*np.array(k2_theta_dot) + 2.0*np.array(k3_theta_dot) + np.array(k4_theta_dot))/6.0
    
    
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

