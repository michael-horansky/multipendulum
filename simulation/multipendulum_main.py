


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



p2_small.plot_mode_portrait()



"""
print("Modal frequency on normal modes compared to normal mode frequencies sanity check:")
p3_small.get_normal_modes()

p3_small.print_normal_modes()

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




triple_p_analyzer = analyzer(0.05, 1.0, 2.5, "middle_freq_no_decay")
triple_p_analyzer.add_pendulum(p3_small, "3p-small")
#triple_p_analyzer.add_pendulum(p3_big, "3p-big")
#triple_p_analyzer.driving_frequency_analysis(driving_frequency_range=(1.5, 2.3), datapoints=100, t_max = 500.0, t_threshold = 50.0)
triple_p_analyzer.load_resonance_analysis_data()
triple_p_analyzer.plot_resonance_analysis_data(mech_resonance_analysis_preset)"""


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




