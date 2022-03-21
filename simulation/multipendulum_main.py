


from class_analyzer import *


my_analyzer = analyzer(0.01, 1.0, 0.0)

jarek_N = 5


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
#my_analyzer.animate(10)

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




