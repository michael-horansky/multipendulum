
"""

A framework to describe a general mechanical system using its
    1. Kinetic and potential term
    2. Equations of motion

What it does: simulates the trajectories; finds response curves to driven motion.

"""

from linalg_funcs import *

import numpy as np
import time


import matplotlib.pyplot as plt

class state:
    
    # The state of a physical system
    
    def __init__(self, q, q_dot):
        
        self.q = q
        self.q_dot = q_dot
        self.N = len(self.q)

def empty_state(N):
    new_empty_state = state(np.zeros(N), np.zeros(N))
    return(new_empty_state)



def force_mode_function(mode, frequency, phase=-1):
    N = len(mode)
    if phase == -1:
        phase = np.zeros(N)
    def force_functor(t):
        result_force_vector = np.zeros(N)
        for i in range(N):
            result_force_vector[i] = mode[i] * np.sin(t * frequency + phase[i])
        return(result_force_vector)
    return(force_functor)

uwu = force_mode_function([5, -2, 1], 0.5)
print(uwu(np.pi))

class mechanical_system:
    
    # ---------------- constructors, destructors, descriptors ----------------------
    
    def __init__(self, N, static_variables, energy_function, inertia_matrix, internal_force_vector):
        
        """
            N[int]: number of degrees of freedom
            static_variables[dict]: an array of variables that don't change during the simulation, but the user might wish to change between simulations
            energy_function[func]: a function that takes (state, static_variables) and returns a tuple (T, V)
            inertia_matrix[func]: a function that takes (state, static_variables) and returns the inertia matrix in the form M_ab = M[a][b]; a,b=0, 1... N-1
            internal_force_vector[func]: a function that takes (state, static_variables) and returns the internal force vector in the form S'_a = S[a]; a=0, 1... N-1
        """
        
        self.N = N
        
        self.static_variables = static_variables.copy()
        
        self.energy_function = energy_function
        self.inertia_matrix = inertia_matrix
        self.internal_force_vector = internal_force_vector
        
        self.set_to_state(empty_state(self.N))
    
    def set_to_state(self, new_state):
        if self.N == new_state.N:
            self.state = new_state
        else:
            print("ERROR: unapplicable state (number of degrees of freedom not matching)")
    
    def set_static_variables(self, new_static_variables):
        
        self.static_variables = new_static_variables.copy()
    
    # TODO add descriptors: print static variables, current state?
    
    
    # ---------------------------------------------------------
    # --------------- physical methods ------------------------
    # ---------------------------------------------------------
    
    def get_total_energy(self):
        T, V = self.energy_function(self.state, self.static_variables)
        return(T+V)
    
    def get_det_M(self):
        my_M = self.inertia_matrix(self.state, self.static_variables)
        return(np.linalg.det(my_M))
    
    def get_S(self, external_force_vector):
        # S_a = F_a(t) + S'_a
        # all in numpy babeyy
        
        cur_S = self.internal_force_vector(self.state, self.static_variables) + external_force_vector
        return(cur_S)
    
    
    def get_acceleration(self, external_force_vector):
        # M a = S
        cur_M = self.inertia_matrix(self.state, self.static_variables)
        cur_S = self.get_S(external_force_vector)
        
        # time complexity O(N^3), which matches that of gaussian elimination
        acceleration = np.linalg.solve(cur_M, cur_S)
        return(acceleration)
    
    
    
    # ---------------------------------------------------------
    # --------------- simulation methods ----------------------
    # ---------------------------------------------------------
    
    def step(self, dt, external_force_list):
        # updates the state according to acceleration
        # RK4 using 3 snapshots of external force
        
        k_0_q     = self.state.q
        k_0_q_dot = self.state.q_dot
        
        acc_1 = self.get_acceleration(external_force_list[0])
        k_1_q     = self.state.q_dot * dt
        k_1_q_dot = acc_1 * dt
        
        self.state.q     = k_0_q     + k_1_q     / 2.0
        self.state.q_dot = k_0_q_dot + k_1_q_dot / 2.0
        acc_2 = self.get_acceleration(external_force_list[1])
        k_2_q     = self.state.q_dot * dt
        k_2_q_dot = acc_2 * dt
        
        self.state.q     = k_0_q     + k_2_q     / 2.0
        self.state.q_dot = k_0_q_dot + k_2_q_dot / 2.0
        acc_3 = self.get_acceleration(external_force_list[1])
        k_3_q     = self.state.q_dot * dt
        k_3_q_dot = acc_3 * dt
        
        self.state.q     = k_0_q     + k_3_q
        self.state.q_dot = k_0_q_dot + k_3_q_dot
        acc_4 = self.get_acceleration(external_force_list[2])
        k_4_q     = self.state.q_dot * dt
        k_4_q_dot = acc_4 * dt
        
        
        self.state.q     = k_0_q     + (k_1_q     + 2.0 * k_2_q     + 2.0 * k_3_q     + k_4_q)    /6.0
        self.state.q_dot = k_0_q_dot + (k_1_q_dot + 2.0 * k_2_q_dot + 2.0 * k_3_q_dot + k_4_q_dot)/6.0
        
    
    def simulate(self, dt, t_thresh, t_max, force_function, return_trajectory = False):
        # returns time-averaged energy
        
        t = 0.0
        E_avg = 0.0
        
        self.set_to_state(empty_state(self.N))
        
        while(t < t_max):
            cur_external_force_list = [force_function(t), force_function(t+dt/2.0), force_function(t+dt)]
            self.step(dt, cur_external_force_list)
            t += dt
            
            if t > t_thresh:
                E_avg += self.get_total_energy() * dt
        
        E_avg /= (t_max - t_thresh)
        
        return(E_avg)
    
    def response_curve(self, omega_range, time_range, force_mode):
        # omega_range = ([omega_min=0], omega_max, N)
        # time_range = ([t_thresh = 0], t_max, dt)
        # force_mode = (mode[N], [phase[N]=0])
        
        # -------------------- overloading ------------------------
        if len(omega_range) == 2:
            omega_min = 0.0
            omega_max, N_datapoints = omega_range
        elif len(omega_range) == 3:
            omega_min, omega_max, N_datapoints = omega_range
        else:
            print("ERROR: omega_range shape invalid")
            return(-1)
        
        if len(time_range) == 2:
            t_thresh = 0.0
            t_max, dt = time_range
        elif len(time_range) == 3:
            t_thresh, t_max, dt = time_range
        else:
            print("ERROR: time_range shape invalid")
            return(-1)
        
        if type(force_mode[0]) != list and type(force_mode[0]) != np.ndarray:
            force_mode_amplitude = force_mode
            force_mode_phase     = [0]*len(force_mode_amplitude)
        elif len(force_mode) == 2:
            force_mode_amplitude = force_mode[0]
            force_mode_phase     = force_mode[1]
        else:
            print("ERROR: force_mode shape invalid")
            return(-1)
        
        # --------------------- response curve calculation --------------
        omega_space = np.linspace(omega_min, omega_max, N_datapoints)
        energy_space = np.zeros(omega_space.shape)
        
        start_time = time.time()
        progress = 0
            
        for i in range(len(omega_space)):
            
            if np.floor(i / N_datapoints * 100) > progress:
                progress = np.floor(i / N_datapoints * 100)
                print(  "  Analysis in progress: " + str(progress) + "%; est. time of finish: " + time.strftime("%H:%M:%S", time.localtime( (time.time()-start_time) * 100 / progress + start_time )), end='\r')
            
            E_avg = self.simulate(dt, t_thresh, t_max, force_mode_function(force_mode_amplitude, omega_space[i], force_mode_phase), False)
            energy_space[i] = E_avg
        
        return(omega_space, energy_space)
        


#------------------------------------------------------------------
#------------ Relevant mechanical system functors -----------------
#------------------------------------------------------------------

# multipendulum_static_variables = {g, l=[l_1, l_2... l_N], m=[m_1, m_2... m_N]}

def multipendulum_energy_function(cur_state, static_variables):
    
    U_g = 0.0
    T = 0.0
    for i in range(cur_state.N):
        U_g_ps = 0.0
        T_ps_1 = 0.0
        T_ps_2 = 0.0
        for j in range(0, i+1):
            U_g_ps += static_variables['l'][j] * np.cos(cur_state.q[j])
            T_ps_1 += cur_state.q_dot[j] * static_variables['l'][j] * np.cos(cur_state.q[j])
            T_ps_2 += cur_state.q_dot[j] * static_variables['l'][j] * np.sin(cur_state.q[j])
        U_g -= static_variables['m'][i] * U_g_ps
        T   += 0.5 * static_variables['m'][i] * (T_ps_1 * T_ps_1 + T_ps_2 * T_ps_2)
    U_g *= static_variables['g']
    return(T, U_g)

def multipendulum_inertia_matrix(cur_state, static_variables):
    M = np.zeros((cur_state.N, cur_state.N))
    for a in range(cur_state.N):
        for b in range(cur_state.N):
            M[a][b] = static_variables['l'][b] * np.cos( cur_state.q[a] - cur_state.q[b] ) * np.sum(static_variables['m'][max(a,b):cur_state.N])#self.get_mu(a, b)
    return(M)

def multipendulum_internal_force_vector(cur_state, static_variables):
    S = np.zeros(cur_state.N)
    #print("q", cur_state.q)
    #print("q_dot", cur_state.q_dot)
    for a in range(cur_state.N):
        cur_sum_i = 0.0
        for i in range(a, cur_state.N):
            cur_sum_j = 0.0
            for j in range(i+1):
                cur_sum_j += static_variables['l'][j] * cur_state.q_dot[j] * cur_state.q_dot[j] * np.sin(cur_state.q[a] - cur_state.q[j])
            cur_sum_i += static_variables['m'][i] * (static_variables['g'] * np.sin(cur_state.q[a]) + cur_sum_j)
        S[a] = -cur_sum_i
    return(S)

MS_p3_small_static = {'g' : 9.8, 'l' : np.array([5.0, 3.0, 3.0]), 'm' : np.array([5.0, 3.0, 3.0])}
MS_p3_small = mechanical_system(3, MS_p3_small_static, multipendulum_energy_function, multipendulum_inertia_matrix, multipendulum_internal_force_vector)

"""
from class_multipendulum import *


p3_small = multipendulum([5.0, 3.0, 3.0], [5.0, 3.0, 3.0], 9.8)

state_q = [0.1, -0.1, 0.2]
state_q_dot = [0.02, 0.05, -0.1]

p3_small.set_state(state_q, state_q_dot)
MS_p3_small.set_to_state(state(state_q, state_q_dot))"""

response_omega, response_energy = MS_p3_small.response_curve((0.1, 4.0, 50), (20.0, 100.0, 0.01), [2.5, 2.5, 2.5])


plt.plot(response_omega, response_energy)
plt.show()



    
    
    
