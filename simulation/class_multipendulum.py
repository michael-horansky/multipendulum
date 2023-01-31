
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend import Legend

import time

from linalg_funcs import *

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
        
        if inner_product(mode, mode) == 0.0:
            return(self.normal_mode_frequencies[0])
        a = 0.0
        b = 0.0
        for i in range(self.N):
            a += self.l[i] * mode[i] * mode[i] * self.get_mu(i)
            c = 0.0
            for j in range(i+1):
                c += mode[j] * self.l[j]
            b += self.m[i] * c * c
        return(np.sqrt(self.g * a / b))
    
    def get_constraint_forces(self, mode):
        
        if inner_product(mode, mode) == 0.0:
            return(mode.copy())
        
        # returns vec(Q) as a function of mode
        mode_freq = self.modal_frequency(mode)
        Q_res = []
        for i in range(self.N):
            cur_Q = self.g * mode[i] * self.get_mu(i)
            for x in range(self.N):
                cur_Q -= mode_freq * mode_freq * mode[x] * self.l[x] * self.get_mu(i, x)
            Q_res.append(self.l[i] * cur_Q)
        return(Q_res)
    
    def get_constraint_force_gradient_matrix(self, mode):
        # returns [[dQ_1/du_1, dQ_1/du_2... dQ_1/du_N], [dQ_2/du_1, dQ_2/du_2... dQ_2/du_N]...]
        P = []
        cur_modal_freq = self.modal_frequency(mode)
        for i in range(self.N):
            P.append([])
            for j in range(self.N):
                if i == j:
                    P[i].append(self.l[i] * (self.g - cur_modal_freq * cur_modal_freq * self.l[j]) * self.get_mu(i) )
                else:
                    P[i].append(- self.l[i] * cur_modal_freq * cur_modal_freq * self.l[j] * self.get_mu(i, j) )
        return(P)
    
    def get_P_element(self, mode, i, j):
        ass_freq = self.modal_frequency(mode)
        res = - ass_freq * ass_freq * self.l[j] * self.get_mu(i, j)
        if i == j:
            res += self.g * self.get_mu(i)
        return(self.l[i] * res)
    
    def get_corrected_resonant_frequencies(self, torque):
        cur_corrected_resonant_frequencies = []
        for n in range(len(self.normal_mode_frequencies)):
            cur_mode = self.normal_modes[n]
            """P = []
            for i in range(self.N):
                P.append([])
                for j in range(self.N):
                    P[-1].append()"""
            P = np.empty((self.N, self.N))
            for i in range(self.N):
                for j in range(self.N):
                    P[i][j] = self.get_P_element(cur_mode, i, j)
            print(P)
            
            print("GRAM SCHMIDT")
            zero_force_span = gram_schmidt(P[1:])
            print(zero_force_span)
            print("LOL", inner_product(zero_force_span[0], zero_force_span[1]))
            delta_u = P[0].copy()#perpendicularize([1]*self.N, P[1:-1])
            print(delta_u)
            print("delta u_0 dot grad 1 =", inner_product(P[0], delta_u))
            
            delta_u_coef = torque / inner_product(P[0], delta_u)
            
            print("Final delta_u =", scalar_product(delta_u, delta_u_coef))
            
            new_u = vector_sum(cur_mode, scalar_product(delta_u, delta_u_coef))
            
            """P_inv = np.linalg.inv(P)
            new_u = cur_mode.copy()
            for i in range(self.N):
                new_u[i] += torque * P_inv[i][0]"""
            cur_corrected_resonant_frequencies.append(self.modal_frequency(new_u))
        
        self.corrected_resonant_frequencies = cur_corrected_resonant_frequencies
    """def get_corrected_resonant_frequencies(self, torque): #TODO: this is incorrect
        cur_corrected_resonant_frequencies = []
        for n in range(len(self.normal_mode_frequencies)):
            omega_n = self.normal_mode_frequencies[n]
            cur_sum = 0
            for j in range(self.N):
                cur_sum += self.normal_modes[n][j] * self.l[j] * self.get_mu(j)
            cur_corrected_resonant_frequencies.append(omega_n - torque / (omega_n * self.l[0] * cur_sum))
        self.corrected_resonant_frequencies = cur_corrected_resonant_frequencies"""
    
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
    
    
    def dp_plot_mode_portrait(self, ranges=[[-1.0, 1.0], [-1.0, 1.0]], graining=256, force_graining = 25):
        # ranges[N][2] = [[u_1_min, u_1_max]...]. Default value -1.0:1.0 everywhere
        # graining = number of datapoints in meshgrid along each dimension
        
        self.get_normal_modes()
        self.print_normal_modes()
        # ASSUME N = 2
        def double_pendulum_modal_frequency(u1, u2):
            top = self.g * (self.l[0] * u1 * u1 * (self.m[0] + self.m[1]) + self.l[1] * u2 * u2 * self.m[1])
            bottom = self.m[0] * u1 * u1 * self.l[0] * self.l[0] + self.m[1] * (u1 * self.l[0] + u2 * self.l[1]) * (u1 * self.l[0] + u2 * self.l[1])
            return(np.sqrt(top / bottom))
        
        # the modal frequency scalar field
        
        x = np.linspace(ranges[0][0],ranges[0][1],graining)
        y = np.linspace(ranges[1][0],ranges[1][1],graining)
        x_q = np.linspace(ranges[0][0],ranges[0][1],force_graining)
        y_q = np.linspace(ranges[1][0],ranges[1][1],force_graining)
        
        xm, ym = np.meshgrid(x, y)
        zm = double_pendulum_modal_frequency(xm, ym)
        
        # the normal modes
        mode_points = []
        for i in range(len(self.normal_modes)):
            mode = self.normal_modes[i]
            mode_points.append([[ranges[i][0], ranges[i][1]], [ranges[i][0] * mode[1] / mode[0], ranges[i][1] * mode[1] / mode[0]]])
        
        
        # the constraint forces
        q1_results = []
        q2_results = []
        colors = []
        for y_val in y_q:
            q1_results.append([])
            q2_results.append([])
            colors.append([])
            for x_val in x_q:
                cur_constraint_forces = self.get_constraint_forces([x_val, y_val])
                u_mag = np.sqrt(x_val * x_val + y_val * y_val)
                q1_results[-1].append(cur_constraint_forces[0])
                q2_results[-1].append(cur_constraint_forces[0])
                #colors[-1].append(a_func(x_val, v_val)+v_val)
                colors[-1].append( 1.0 )
        
        # the constraint force gradients
        x_list = []
        y_list = []
        q_1_list = [[], []]
        q_2_list = [[], []]
        for i in range(len(self.normal_modes)):
            mode = self.normal_modes[i]
            x_list.append([])
            y_list.append([])
            q_1_list[0].append([])
            q_2_list[0].append([])
            q_1_list[1].append([])
            q_2_list[1].append([])
            for x_val_1 in x_q:
                x_val_2 = x_val_1 + 0.5 * (x_q[1] - x_q[0])
                y_val_1 = x_val_1 * mode[1] / mode[0]
                y_val_2 = x_val_2 * mode[1] / mode[0]
                x_list[-1].append(x_val_1)
                y_list[-1].append(y_val_1)
                cur_P_1 = self.get_constraint_force_gradient_matrix([x_val_1, y_val_1])
                cur_P_2 = self.get_constraint_force_gradient_matrix([x_val_2, y_val_2])
                q_1_list[0][-1].append(cur_P_1[0][0])
                q_1_list[1][-1].append(cur_P_1[0][1])
                q_2_list[0][-1].append(cur_P_2[1][0])
                q_2_list[1][-1].append(cur_P_2[1][1])
            
        fig, ax = plt.subplots()
        #plt.figure(figsize=(12, 8))
        plt.xlim(ranges[0][0], ranges[0][1])
        plt.ylim(ranges[1][0], ranges[1][1])
        plt.xlabel('$u_1$')
        plt.ylabel('$u_2$')
        pcm = plt.pcolormesh(xm, ym, zm, cmap='RdBu_r')
        
        isofrequenth_colors = ['#ff6699', '#66ff33']
        
        lines = []
        nm_labels = []
        
        for i in range(len(mode_points)):
            mode = mode_points[i]
            nm_labels.append(f'$v_{i+1}, \\omega_{i+1}={self.normal_mode_frequencies[i]:.2f}$ rad.s$^{{-1}}$')
            cur_plot, = plt.plot(mode[0], mode[1], linestyle='solid', color=isofrequenth_colors[i], label=nm_labels[-1])
            lines.append(cur_plot)
        plt.colorbar(pcm, label='$\omega(\\vec{u})$')
        #plt.clim(min(self.normal_mode_frequencies), max(self.normal_mode_frequencies))
        lines.append(plt.quiver(x_q, y_q, q1_results, q2_results, label='$\\vec{{Q}}$'))
        
        for i in range(len(self.normal_modes)):
            if i == 0:
                lines.append(plt.quiver(x_list[i], y_list[i], q_1_list[0][i], q_1_list[1][i], color='white', edgecolor='red', linewidth = 1, label='$\\nabla Q_1$'))
                lines.append(plt.quiver(x_list[i] + 0.5 * (x_q[1] - x_q[0]), y_list[i] + 0.5 * (x_q[1] - x_q[0]) * self.normal_modes[i][1] / self.normal_modes[i][0], q_2_list[0][i], q_2_list[1][i], color='white', edgecolor='black', linewidth = 1, label='$\\nabla Q_2$'))
            else:
                plt.quiver(x_list[i], y_list[i], q_1_list[0][i], q_1_list[1][i], color='white', edgecolor='red', linewidth = 1)
                plt.quiver(x_list[i] + 0.5 * (x_q[1] - x_q[0]), y_list[i] + 0.5 * (x_q[1] - x_q[0]) * self.normal_modes[i][1] / self.normal_modes[i][0], q_2_list[0][i], q_2_list[1][i], color='white', edgecolor='black', linewidth = 1)
        
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width, pos.height * 1.0])
        ax.legend(loc='center right', ncol=self.N + 3, frameon=False, columnspacing=0.8, bbox_to_anchor=(1.27, 1.1))
        #ax.legend(lines[:self.N], nm_labels, loc='center right', ncol=self.N, frameon=False, columnspacing=0.8, bbox_to_anchor=(1.27, 1.1))
        #leg = Legend(ax, lines[self.N:], ['$\\vec{{Q}}$', '$\\nabla Q_1$', '$\\nabla Q_2$'], loc='center right', ncol=3, frameon=False, columnspacing=0.8, bbox_to_anchor=(1.27, 1.05))
        #ax.add_artist(leg)
        plt.show()
    
    def tp_plot_mode_portrait(self, R=0.01, graining=256, force_graining = 25):
        # mode portrait for a triple pendulum - take a sphere of radius R and project onto a 2D surface
        # x = phi, y = theta
        
        # ASSUME N = 3
        def spherical_modal_frequency(phi, theta):
            u1 = R * np.sin(theta) * np.cos(phi)
            u2 = R * np.sin(theta) * np.sin(phi)
            u3 = R * np.cos(theta)
            
            top = self.g * (self.l[0] * u1 * u1 * (self.m[0] + self.m[1] + self.m[2]) + self.l[1] * u2 * u2 * (self.m[1] + self.m[2]) + self.l[2] * u3 * u3 * self.m[2])
            bottom = self.m[0] * u1 * u1 * self.l[0] * self.l[0] + self.m[1] * (u1 * self.l[0] + u2 * self.l[1]) * (u1 * self.l[0] + u2 * self.l[1]) + self.m[2] * (u1 * self.l[0] + u2 * self.l[1] + u3 * self.l[2]) * (u1 * self.l[0] + u2 * self.l[1] + u3 * self.l[2])
            return(np.sqrt(top / bottom))
        
        def c_to_s_transform(v):
            return([np.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]), np.arctan2(np.sqrt(v[0]*v[0]+v[1]*v[1]),v[2]), np.arctan2(v[1],v[0])])
        def s_to_c_transform(v):
            return([v[0] * np.sin(v[1]) * np.cos(v[2]), v[0] * np.sin(v[1]) * np.sin(v[2]), v[0] * np.cos(v[1])])
        
        def spherical_constraint_forces(phi, theta):
            cur_mode = s_to_c_transform(R, theta, phi)
            cur_constraint_forces = self.get_constraint_forces(cur_mode)
            c_forces = c_to_s_transform(cur_constraint_forces)
            return([c_forces[2], c_forces[1]]) #[phi, theta]
        
        def old_spherical_constraint_force_gradient(phi, theta):
            cur_mode = s_to_c_transform([R, theta, phi])
            print(cur_mode)
            cur_constraint_force_gradient_matrix = self.get_constraint_force_gradient_matrix(cur_mode)
            print(cur_constraint_force_gradient_matrix)
            c_force_gradient_matrix = []
            for i in range(len(cur_constraint_force_gradient_matrix)):
                cur_force_gradient = c_to_s_transform(cur_constraint_force_gradient_matrix[i])
                c_force_gradient_matrix.append([cur_force_gradient[2], cur_force_gradient[1]])
            return(c_force_gradient_matrix)
        
        def log_v(v):
            # rescales vector in log scale
            res = v.copy()
            v_size = np.sqrt(inner_product(v, v))
            v_log_size = np.log(v_size)
            for i in range(len(v)):
                res[i] = v[i] * v_log_size / v_size
            return(res)
            # res = v / |v| * log(|v|)
        
        def spherical_constraint_force_gradient(phi, theta):
            cur_mode = s_to_c_transform([R, theta, phi])
            cur_constraint_force_gradient_matrix = self.get_constraint_force_gradient_matrix(cur_mode)
            c_force_gradient_matrix = []
            for i in range(len(cur_constraint_force_gradient_matrix)):
                grad_q_i = cur_constraint_force_gradient_matrix[i]
                cur_phi = - grad_q_i[0] * np.sin(phi) + grad_q_i[1] * np.cos(theta)
                cur_theta = grad_q_i[0] * np.cos(phi) * np.cos(theta) + grad_q_i[1] * np.sin(phi) * np.cos(theta) - grad_q_i[2] * np.sin(theta)
                c_force_gradient_matrix.append(log_v([cur_phi, cur_theta]))
            return(c_force_gradient_matrix)
        
        def plus_pi(x):
            res = x + np.pi
            if res >= 2.0 * np.pi:
                res -= 2.0 * np.pi
            return(res)
        
        self.get_normal_modes()
        self.print_normal_modes()
            
        # the modal frequency scalar field
        
        phi_space   = np.linspace(0.0,2.0 * np.pi,graining)
        theta_space = np.linspace(0.0,np.pi      ,graining)
        #x_q = np.linspace(ranges[0][0],ranges[0][1],force_graining)
        #y_q = np.linspace(ranges[1][0],ranges[1][1],force_graining)
        
        phi_mesh, theta_mesh = np.meshgrid(phi_space, theta_space)
        omega_mesh = spherical_modal_frequency(phi_mesh, theta_mesh)
        
        # the normal modes and force gradient
        normal_mode_phi_list = []
        normal_mode_theta_list = []
        
        gq_phi_list = []
        gq_theta_list = []
        grad_q_phi_list = [[], [], []]
        grad_q_theta_list = [[], [], []]
        for i in range(len(self.normal_modes)):
            mode = self.normal_modes[i]
            s_mode = c_to_s_transform(mode)
            
            grad_q_m_1 = spherical_constraint_force_gradient(s_mode[2], s_mode[1])
            grad_q_m_2 = spherical_constraint_force_gradient(plus_pi(s_mode[2]), np.pi-s_mode[1])
            
            normal_mode_phi_list.append([s_mode[2], plus_pi(s_mode[2])])
            normal_mode_theta_list.append([s_mode[1], np.pi-s_mode[1]])
            
            gq_phi_list.append(s_mode[2])
            gq_theta_list.append(s_mode[1])
            gq_phi_list.append(plus_pi(s_mode[2]))
            gq_theta_list.append(np.pi-s_mode[1])
            
            for i in range(self.N):
                grad_q_phi_list[i].append(grad_q_m_1[i][0])
                grad_q_theta_list[i].append(grad_q_m_1[i][1])
                grad_q_phi_list[i].append(grad_q_m_2[i][0])
                grad_q_theta_list[i].append(grad_q_m_2[i][1])
        
        fig, ax = plt.subplots()
        plt.xlim(0.0, 2.0 * np.pi)
        plt.ylim(0.0, np.pi)
        plt.xlabel('$\\phi$ (equal to $\\arctan(u_y/u_x)$)')
        plt.ylabel('$\\theta$ (equal to $\\arccos(u_z/|\\vec{u}|)$)')
        
        pcm = plt.pcolormesh(phi_mesh, theta_mesh, omega_mesh, cmap='RdBu_r')
        plt.colorbar(pcm, label='$\omega(\\vec{u})$')
        
        normal_mode_colours = ['#ff6699', '#66ff33', '#009999']
        
        for i in range(len(self.normal_modes)):
            nm_label = f'$v_{i+1}, \\omega_{i+1}={self.normal_mode_frequencies[i]:.2f}$ rad.s$^{{-1}}$'
            plt.scatter(normal_mode_phi_list[i], normal_mode_theta_list[i], s=80, label=nm_label, color=normal_mode_colours[i], edgecolor='black', linewidth = 1)
        
        grad_q_colors = ['red', 'black', 'yellow']
        for i in range(self.N):
            """if i == 0:
                    plt.quiver(gq_phi_list, y_list[i], q_1_list[0][i], q_1_list[1][i], color='white', edgecolor='red', linewidth = 1, label='$\\nabla Q_1$')
                    plt.quiver(x_list[i] + 0.5 * (x_q[1] - x_q[0]), y_list[i] + 0.5 * (x_q[1] - x_q[0]) * self.normal_modes[i][1] / self.normal_modes[i][0], q_2_list[0][i], q_2_list[1][i], color='white', edgecolor='black', linewidth = 1, label='$\\nabla Q_2$')
                else:"""
            plt.quiver(gq_phi_list, gq_theta_list, grad_q_phi_list[i], grad_q_theta_list[i], color='white', edgecolor=grad_q_colors[i], linewidth = 1, scale = 50, label=f'$\\nabla Q_{i+1}$ (log scale)')
        
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width, pos.height * 1.0])
        ax.legend(loc='center right', ncol=self.N, frameon=False, columnspacing=0.8, bbox_to_anchor=(1.1, 1.08))
        plt.show()
                
    
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

