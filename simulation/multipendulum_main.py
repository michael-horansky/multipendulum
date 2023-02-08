


from class_analyzer import *

# plotting presets
ang_f_analysis_preset = ['avg_E_graph', 'ang_f_theta_graph', 'ang_f_theta_space', 'hp_theta_std_graph', 'ang_f_phi_graph', 'ang_f_phi_space', 'hp_phi_std_graph']

mech_resonance_analysis_preset = ['avg_E_graph', 'avg_L_graph', 'max_theta_1_graph'] #how about avg_T? could be interesting. also max_L


my_analyzer = analyzer(0.01, 1.0, 0.0)

jarek_N = 5


p2_small = multipendulum([5.0, 3.0], [5.0, 3.0], 9.8)
p2_333   = multipendulum([5.0, 5.0], [5.0, 5.0], 9.8)
p3_small = multipendulum([5.0, 3.0, 3.0], [5.0, 3.0, 3.0], 9.8)
p3_big   = multipendulum([5.0, 5.0, 5.0], [5.0, 5.0, 5.0], 9.8)


# in state memory, first line is the driving frequency raneg etc, then every line is one state value (state of each pendulum separated by |)
"""
energy_vs_amplitude = analyzer(0.01, 1.0, 2.5, "middle_freq_decay_10dp")

#energy_vs_amplitude.add_pendulum(p2_small, "double-pendulum")
energy_vs_amplitude.add_pendulum(p3_small, "triple-pendulum")
#energy_vs_amplitude.add_pendulum(p3_big  , "p3 big"  )

#energy_vs_amplitude.load_state(25.0, 20.0)
#energy_vs_amplitude.driving_frequency_analysis(driving_frequency_range=(0.1, 4.0), datapoints=10, t_max = 20.0, t_threshold=15.0, overwrite_stored_states=False)

energy_vs_amplitude.driving_frequency_analysis(driving_frequency_range=(2.0, 2.6), datapoints = 5, t_max = 200.0, t_threshold=150.0)
#energy_vs_amplitude.save_resonance_analysis_data()
energy_vs_amplitude.plot_resonance_analysis_data(save_graph=False)

"""

"""avg_L_analyzer = analyzer(0.02, 1.0, 2.5, "avg_L_analyzer")
avg_L_analyzer.add_pendulum(p2_small, "small pendulum")

avg_L_analyzer.driving_frequency_analysis(driving_frequency_range=(1.9, 2.9), datapoints=100, t_max = 100.0, t_threshold = 50.0)
avg_L_analyzer.save_resonance_analysis_data()
avg_L_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset)"""



# we select mode 1 and create a small displacement in u space perpendicular to v_1:
"""perturbation_index = 0
displacement_scale = 0.3
mode_of_interest = p2_small.normal_modes[perturbation_index]
mode_displacement = [- mode_of_interest[1] * displacement_scale, mode_of_interest[0] * displacement_scale]
target_mode = vector_sum(mode_of_interest, mode_displacement)
target_frequency = p2_small.modal_frequency(target_mode)
print(f"New mode = ({mode_of_interest[0]}, {mode_of_interest[1]}) + ({mode_displacement[0]}, {mode_displacement[1]}) = ({mode_of_interest[0] + mode_displacement[0]}, {mode_of_interest[1] + mode_displacement[1]})")
print(f"Old frequency = {p2_small.normal_mode_frequencies[perturbation_index]}; new frequency = {target_frequency}")
target_force_mode = p2_small.get_constraint_forces(target_mode)
print(f"Associated constraint force mode = [{target_force_mode[0]:.02f}, {target_force_mode[1]:.02f}]")"""



#NOTE: looks like taking the perturbed force mode from mode 1 and calculating the perturbed frequency from mode 2 predicts well the resonance 2 with force 1. WHY

#force_magnitude = 2.5

def double_pendulum_mode_perturbation(double_pendulum, force_magnitude, analyze = True):
    
    torque_magnitude = force_magnitude * double_pendulum.l[0]
    print("Torque magnitude associated with inputted force magnitude =", torque_magnitude)
    
    dp_analyzer = analyzer(0.05, 1.0, 2.5, "dp_force_scalar")
    dp_analyzer.add_pendulum(double_pendulum, "2pend")
    if analyze:
        dp_analyzer.driving_frequency_analysis(driving_frequency_range=(0.9, 2.9), force_mode = torque_magnitude / double_pendulum.l[0], datapoints=100, t_max = 1000.0, t_threshold = 200.0)
        dp_analyzer.save_resonance_analysis_data()
    else:
        dp_analyzer.load_resonance_analysis_data()
        dp_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset)
    
    for n_m_i in range(2):
        target_frequency, target_torque_mode, used_displacement = double_pendulum.get_perturbed_force_mode(n_m_i, 0.1, print_values = False)
        
        mode_displacement_coef = torque_magnitude / magnitude(target_torque_mode)
        new_displacement = scalar_product(used_displacement, mode_displacement_coef)
        
        target_frequency, target_torque_mode, used_displacement = double_pendulum.get_perturbed_force_mode(n_m_i, new_displacement)

        dp_analyzer = analyzer(0.05, 1.0, 2.5, "dp_force_vector_%i" % (n_m_i + 1))
        dp_analyzer.add_pendulum(double_pendulum, "2pend")
        if analyze:
            dp_analyzer.driving_frequency_analysis(driving_frequency_range=(0.9, 2.9), force_mode = double_pendulum.get_force_from_torque(target_torque_mode), datapoints=100, t_max = 1000.0, t_threshold = 200.0)
            dp_analyzer.save_resonance_analysis_data()
        else:
            dp_analyzer.load_resonance_analysis_data()
            dp_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset, reference_frequencies=[target_frequency], reference_frequencies_label='displaced modal f.')
        print("-----------------------------------------------")

def double_pendulum_force_matching_mode(double_pendulum, force_magnitude, analyze = True):
    
    torque_magnitude = force_magnitude * double_pendulum.l[0]
    print("Torque magnitude associated with inputted force magnitude =", torque_magnitude)
    
    double_pendulum.modal_analysis()
    dp_analyzer = analyzer(0.05, 1.0, 2.5, "dp_v_n_force_scalar")
    dp_analyzer.add_pendulum(double_pendulum, "2pend")
    if analyze:
        dp_analyzer.driving_frequency_analysis(driving_frequency_range=(0.9, 2.9), force_mode = torque_magnitude / double_pendulum.l[0], datapoints=100, t_max = 1000.0, t_threshold = 200.0)
        dp_analyzer.save_resonance_analysis_data()
    else:
        dp_analyzer.load_resonance_analysis_data()
        dp_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset)
    
    for n_m_i in range(2):
        target_torque_mode = double_pendulum.normal_modes[n_m_i].copy()
        
        mode_displacement_coef = torque_magnitude / magnitude(target_torque_mode)
        target_torque_mode = scalar_product(target_torque_mode, mode_displacement_coef)
        print(f"Magnitude of torque mode = {magnitude(target_torque_mode):.3f}, magnitude of force mode = {magnitude(double_pendulum.get_force_from_torque(target_torque_mode)):.3f}, ")
        
        dp_analyzer = analyzer(0.05, 1.0, 2.5, "dp_v_n_force_vector_%i" % (n_m_i + 1))
        dp_analyzer.add_pendulum(double_pendulum, "2pend")
        if analyze:
            dp_analyzer.driving_frequency_analysis(driving_frequency_range=(0.9, 2.9), force_mode = double_pendulum.get_force_from_torque(target_torque_mode), datapoints=100, t_max = 1000.0, t_threshold = 200.0)
            dp_analyzer.save_resonance_analysis_data()
        else:
            dp_analyzer.load_resonance_analysis_data()
            dp_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset)
        print("-----------------------------------------------")



def triple_pendulum_mode_perturbation(triple_pendulum, force_magnitude, omega_range = (0.6, 4.2), analyze = True):

    triple_pendulum.modal_analysis()
    dp_analyzer = analyzer(0.05, 1.0, 2.5, "tp_force_scalar")
    dp_analyzer.add_pendulum(triple_pendulum, "3pend")
    if analyze:
        dp_analyzer.driving_frequency_analysis(driving_frequency_range=omega_range, force_mode = force_magnitude, datapoints=100, t_max = 500.0, t_threshold = 100.0)
        dp_analyzer.save_resonance_analysis_data()
    else:
        dp_analyzer.load_resonance_analysis_data()
        dp_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset)

    print("-----------------------------------------------")
    for n_m_i in range(3):
        initial_displacement = scalar_product(vector_product([1.0, 0.0, 0.0], triple_pendulum.normal_modes[n_m_i]), 1/5.0**(n_m_i))
        target_frequency, target_force_mode, used_displacement = triple_pendulum.get_perturbed_force_mode(n_m_i, initial_displacement, print_values = True)
        
        mode_displacement_coef = force_magnitude / magnitude(target_force_mode)
        new_displacement = scalar_product(used_displacement, mode_displacement_coef)
        
        target_frequency, target_force_mode, used_displacement = triple_pendulum.get_perturbed_force_mode(n_m_i, new_displacement)

        dp_analyzer = analyzer(0.05, 1.0, 2.5, "tp_force_vector_%i" % (n_m_i + 1))
        dp_analyzer.add_pendulum(triple_pendulum, "3pend")
        if analyze:
            dp_analyzer.driving_frequency_analysis(driving_frequency_range=omega_range, force_mode = target_force_mode, datapoints=100, t_max = 500.0, t_threshold = 100.0)
            dp_analyzer.save_resonance_analysis_data()
        else:
            dp_analyzer.load_resonance_analysis_data()
            dp_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset, reference_frequencies=[target_frequency], reference_frequencies_label='displaced modal f.')
        print("-----------------------------------------------")


"""dp_analyzer = analyzer(0.05, 1.0, 2.5, "3p_freq-mode_displacement")
dp_analyzer.add_pendulum(p3_small, "p3-small")
#dp_analyzer.add_pendulum(p3_big, "p3-big")
my_force_mode = [2.5, 0.0, 0.0]
dp_analyzer.driving_frequency_analysis(driving_frequency_range=(0.6, 4.2), force_mode = my_force_mode, datapoints=200, t_max = 500.0, t_threshold = 100.0)
dp_analyzer.save_resonance_analysis_data()

p3_small.get_corrected_resonant_frequencies(my_force_mode, scale = 10.0)
#p3_big.get_corrected_resonant_frequencies(my_force_mode, scale = 10.0)
dp_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset, reference_frequencies=p3_small.corrected_resonant_frequencies, reference_frequencies_label='corrected f.')
"""

"""
# TESTING whether corrected freq are scale invariant
c_space = np.linspace(-5.0, 5.1, 503)
omega_lists = [[], [], []]
p3_small.modal_analysis()
for c in c_space:
    cur_omega = p3_small.get_corrected_resonant_frequencies([2.5, 0.0, 0.5], c)
    for i in range(p3_small.N):
        omega_lists[i].append(cur_omega[i])

plt.xlabel('c')
plt.ylabel('$\\omega$ [rad.s$^{{-1}}$]')
plt.ylim((-0.1, 4.7))
plots = []
for i in range(p3_small.N):
    plots.append(plt.plot(c_space, omega_lists[i], label = f'$\\omega_{i+1} + \\delta s_{{\\omega}}(c\\vec{{w}})$'))
    plt.axhline(y=p3_small.normal_mode_frequencies[i], linestyle='dotted', label = f'$\\omega_{i+1}$', color = plots[i][0].get_color())
plt.legend()
plt.tight_layout()
plt.show()"""

#double_pendulum_mode_perturbation(p2_small, 2.5, analyze = False)
#double_pendulum_force_matching_mode(p2_small, 2.5, analyze = False)
#triple_pendulum_mode_perturbation(p3_small, 2.5, analyze=False)
p3_small.tp_plot_mode_portrait()

"""

a1 = p3_small.modal_frequency(p3_small.normal_modes[0])
a2 = p3_small.modal_frequency(p3_small.normal_modes[1])
a3 = p3_small.modal_frequency(p3_small.normal_modes[2])
print(f"modal freq. analysis: omega_1 = {a1}, omega_2 = {a2}, omega_3 = {a3}")

q1 = p3_small.get_constraint_forces(np.array(p3_small.normal_modes[0]))
q2 = p3_small.get_constraint_forces(p3_small.normal_modes[1])
q3 = p3_small.get_constraint_forces(p3_small.normal_modes[2])
print(f"constraint forces analysis: Q_1 = {q1}, Q_2 = {q2}, Q_3 = {q3}")

p3_small.get_corrected_resonant_frequencies(2.5)
print("CORRECTED", p3_small.corrected_resonant_frequencies)
"""

"""
print("Modal frequency on normal modes compared to normal mode frequencies sanity check:")
p3_big.get_normal_modes()

p3_big.print_normal_modes()





triple_p_analyzer = analyzer(0.05, 1.0, 2.5, "triple_pend_full_spectrum")
triple_p_analyzer.add_pendulum(p3_small, "3p-small")
triple_p_analyzer.add_pendulum(p3_big, "3p-big")
#triple_p_analyzer.driving_frequency_analysis(driving_frequency_range=(0.5, 4.1), datapoints=300, t_max = 500.0, t_threshold = 50.0)
triple_p_analyzer.load_resonance_analysis_data()
triple_p_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset)
#triple_p_analyzer.save_resonance_analysis_data()

"""


"""
p3_small.random_state(np.pi)
energy_vs_amplitude.animate(50.0, 40.0)"""

#resonant_frequencies = [1.065326633165829, 2.170854271356784, 3.5376884422110555]
#resonant frequencies = [1.0603015075376885, 2.1404522613065327, 2.2384422110552764, 3.549246231155779]
"""
for resonant_frequency in energy_vs_amplitude.resonant_frequencies[0]:
#for resonant_frequency in resonant_frequencies:
    energy_vs_amplitude.omega_F = resonant_frequency
    energy_vs_amplitude.animate(160.0, 150.0)
"""

"""
double_pendulum_resonance = analyzer(0.01, 1.0, 2.5)

double_pendulum_resonance.add_pendulum(p2_333, "p_333")

#double_pendulum_resonance.driving_frequency_analysis(driving_frequency_range=(0.0, 4.0), t_max = 100.0, t_threshold=20.0)
#double_pendulum_resonance.driving_frequency_analysis(driving_frequency_range=(0.1, 4.0), datapoints = 10, t_max = 50.0, t_threshold=20.0)
double_pendulum_resonance.load_resonance_analysis_data()
double_pendulum_resonance.plot_resonance_analysis_data()
"""
"""
for resonant_frequency in double_pendulum_resonance.resonant_frequencies[0]:
    double_pendulum_resonance.omega_F = resonant_frequency
    double_pendulum_resonance.animate(50.0, 40.0)
"""

"""

atak_bobra = "girth"

jacob = multipendulum(np.ones(4), [1, 10, 100, 1000], 9.8)


jacob.random_state()
print(jacob.get_acceleration(0.0))

p1 = multipendulum(np.ones(1), np.ones(1), 9.8)
p1b = multipendulum(np.ones(1)*2, np.ones(1)*2, 9.8)
p2 = multipendulum(np.ones(2), np.ones(2), 9.8)

my_analyzer.add_pendulum(p1)
my_analyzer.add_pendulum(p1b)
my_analyzer.add_pendulum(p2)

p1.set_state([np.pi / 2.0], [0.0])
p1b.set_state([1.0], [0.0])
p2.set_state([1.0, 0.0], [0.0, 0.0])
my_analyzer.animate(10)

balls = analyzer(0.01, 1.0, 0.0)
for i in range(5):
    balls.add_pendulum(multipendulum(np.ones(3), np.ones(3), 9.8))
    balls.pendulums[i].set_state([1.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    if i != 0:
        balls.pendulums[i].random_state(np.pi)
    
#balls.animate(20)

smak = analyzer(0.01, 1.0, 0.0)
smak.add_pendulum(multipendulum(np.ones(4), np.ones(4), 9.8))
smak.pendulums[0].random_state(np.pi)
smak.animate(20)
"""

"""
jacob.set_M_functor()

jacob.print_C_M()

jacob.random_state()
print(jacob.get_det_M())
print(jacob.get_det_M_opt())


# benchmark performance test

jarek = multipendulum(np.ones(jarek_N) + np.random.rand(jarek_N) * 1.0, np.ones(jarek_N) + np.random.rand(jarek_N) * 1.0, 9.8)
jarek.set_M_functor()
print("Number of terms for N=%i: %i" % (jarek_N, jarek.get_C_M_size()))

print("Comparative performance test: N=%i, 1000 random cycles" % (jarek_N))
start_time = time.time()
for i in range(1000):
    #if i % 100 == 0:
    #    print(i/10, "% done")
    jarek.random_state()
    jarek.get_det_M()
print("Direct determinant computation: %s seconds" % (time.time() - start_time))

start_time = time.time()
for i in range(1000):
    #if i % 100 == 0:
    #    print(i/10, "% done")
    jarek.random_state()
    jarek.get_det_M_opt()
print("Prepared determinant computation: %s seconds" % (time.time() - start_time))"""




