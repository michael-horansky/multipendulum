


from class_analyzer import *


my_analyzer = analyzer(0.01, 1.0, 0.0)

jarek_N = 5


p2_small = multipendulum([5.0, 3.0], [5.0, 3.0], 9.8)
p3_small = multipendulum([5.0, 3.0, 3.0], [5.0, 3.0, 3.0], 9.8)


energy_vs_amplitude = analyzer(0.01, 1.0, 2.5)

#energy_vs_amplitude.add_pendulum(p2_small, "double-pendulum")
energy_vs_amplitude.add_pendulum(p3_small, "triple-pendulum")

#energy_vs_amplitude.driving_frequency_analysis(driving_frequency_range=(0.0, 4.0), t_max = 50.0, t_threshold=20.0)
#energy_vs_amplitude.plot_resonance_analysis_data(['max_theta_1_graph', 'avg_E_graph', 'avg_E_ratio_graph'])

resonant_frequencies = [1.065326633165829, 2.170854271356784, 3.5376884422110555]

#for resonant_frequency in energy_vs_amplitude.resonant_frequencies[0]:
for resonant_frequency in resonant_frequencies:
    energy_vs_amplitude.omega_F = resonant_frequency
    energy_vs_amplitude.animate(160.0, 150.0)


"""
double_pendulum_resonance = analyzer(0.01, 1.0, 2.5)

double_pendulum_resonance.add_pendulum(p2_small, "p_small")

#double_pendulum_resonance.driving_frequency_analysis((0.0, 2.5))
#double_pendulum_resonance.driving_frequency_analysis(driving_frequency_range=(0.0, 4.0), t_max = 100.0, t_threshold=20.0)
double_pendulum_resonance.driving_frequency_analysis(driving_frequency_range=(0.0, 4.0), t_max = 15.0, t_threshold=5.0)
double_pendulum_resonance.plot_resonance_analysis_data(['max_theta_1_graph', 'avg_E_graph', 'avg_E_ratio_graph'])

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




