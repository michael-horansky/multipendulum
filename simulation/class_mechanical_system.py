
"""

A framework to describe a general mechanical system using its
    1. Kinetic and potential term
    2. Equations of motion

What it does: simulates the trajectories; finds response curves to driven motion.

"""

class state:
    
    # The state of a physical system
    
    def __init__(self, q, q_dot):
        
        self.q = q
        self.q_dot = q_dot
        self.N = len(self.q)

def empty_state(N):
    new_empty_state = state([0.0]*N, [0.0]*N)
    return(new_empty_state)

class mechanical_system:
    
    def __init__(self, N, static_variables, energy_function, inertia_matrix, internal_force_vector):
        
        """
            N: number of degrees of freedom
            static_variables: an array of variables that don't change during the simulation, but the user might wish to change between simulations
            energy_function: a function that takes (state, static_variables) and returns a tuple (K, V)
            inertia_matrix: a function that takes (state, static_variables) and returns the inertia matrix in the form M_ab = M[a][b]; a,b=0, 1... N-1
            internal_force_vector: a function that takes (state, static_variables) and returns the internal force vector in the form S'_a = S[a]; a=0, 1... N-1
        """
        
        self.N = N
        self.energy_function = energy_function
        self.inertia_matrix = inertia_matrix
        self.internal_force_vector = internal_force_vector
        
        self.set_to_state(empty_state(self.N))
    
    def set_to_state(self, new_state):
        if self.N == new_state.N:
            self.state = new_state
        else:
            print("ERROR: unapplicable state (number of degrees of freedom not matching)")
    
    def 
    
    
    
