
# Introduction and Structure

This software provides a bulk of functionalities to perform driving frequency analysis on multipendulums. This means we drive the upper segment of a multipendulum with a sine torque of a given frequency and look at aggregate properties of the system after a given amount of time, such as the time-average mechanical energy stored in the system or the maximal amplitude of the first-segment angle. This data can be then analyzed further to extract resonant frequencies and perform further analysis on them, like checking the mode coefficients.

The software provided consists of the following files: multipendulum_main.py, class_analyzer.py, class_multipendulum.py. These have the following purposes:
- **multipendulum_main.py:** This is the 'command file'. It contains instances of available classes and calls all their methods that define the experiment you're conducting. This is the file you need to edit to adjust your experiment, from parameters of the multipendulum to the range of driving frequency or other parameter values you're doing a comparative analysis on. **To use the program, you need to run this file with python3**, no flags necessary. This file shouldn't include any function or class definitions.
- **class_analyzer.py:** This is a class file for the _analyzer_ class. The _analyzer_ is a class described by a single bulk of pendulums and the global parameters independent of any pendulum, such as global time. Its purpose is to conduct comparative analysis on its pendulums, with the primary method to do so being _driving_frequency_analysis_. It can also analyze the output from the analysis to find resonant frequencies, check their mode coefficients, or generate a pendulum population that spans a given parameter space.
- **class_multipendulum.py:** This is a class file for the _multipendulum_ class. An instance of this class is a single multipendulum, described by its static properties (the lengths of its segments, the masses of its nodes, the gravitational acceleration) and its dynamic properties (angles and angular velocities of each segment). This class provides methods to obtain the equations of motion of a multipendulum and to perform numerical integration of these equations using RK4.

These files link to each other as listed above in a descending order. You need to place these files in a single folder to run the program. Then run multipendulum_main.py in the environment of your choice.

# Conducting an Experiment

An 'experiment' is a set of driving frequency analyses conducted on a specified population of multipendulums on a specified range of frequencies. First you need to set up the experiment's configuration, then run it, and finally, use the data.

## Setting up the configuration

First, you create any amount of analyzers as instances of the _analyzer_ class. The constructor takes 3 parameters:
1. my_dt (float): the step used in numerical integration, of dimension [s]
2. my_omega_F (float): the *default value* of the driving force angular frequency.
3. my_external_force_amplitude (float): the *default value* of the driving force amplitude.
The default values will be used in other methods called for this instance if a different value isn't specified.

### Manual population

Second, you create the pendulums in the experiment as instances of the _multipendulum_ class. The constructor takes 3 parameters:
1. my_l (array): list of lengths of the segments. The constructor will infer the amount of segments N from the length of this list.
2. my_m (array): list of masses of the nodes.
3. my_g (float): the gravitational acceleration acting on the multipendulum.

Now you add the pendulums to the analyzers they belong to, by calling the analyzers' method add_pendulum. This method takes two parameters:
1. pendulum (multipendulum): the instance of _multipendulum_ you're adding.
2. pendulum_name (string) [optional]: the name of the pendulum, used to mark it on plots and in filenames of exported data. If unspecified, the name will be "pendulum_#", where # is the current amount of pendulums added to the analyzer, which increments by 1 with each time this method's called.

### Automatic parameter-space population

Alternatively, you can skip the pendulum creation and addition and call the analyzer's method populate_parameter_space. This method creates a population of pendulums whose parameters l, m, g span the ranges provided. This method takes 4 parameters:
1. l_space (array of tuples or floats): list of tuples of length N, where the i-th element is a tuple in the form (start, end, amount[, pivot=0]), where _start_ is the lower boundary on the value of l_i, _end_ is the upper boundary, and _amount_ is the number of values we consider on this interval. If a single number is presented instead of the tuple, then l_i only takes this value and doesn't vary in the population. If pivot_expansion is flagged, then _pivot_ marks the index of the pivot value in the produced linspace.
2. m_space (array of tuples or floats): the same as l_space, but for masses of the nodes, so the i-th element describes the space of m_i.
3. g_space (tuple or float): if a tuple in the form (start, end, amount[, pivot=0]), g will take values from _start_ to _end_ with _amount_ datapoints. If a float, g always takes the value of the float. If pivot_expansion is flagged, _pivot_ marks the index of the pivot value in the produced linspace.
4. pivot_expansion (bool) [optional, False by default]: If false, the population will include all possible configurations given by the parameter spaces. If true, only configurations in which all values except one are equal to their respective pivot values are considered (so only one parameter is being varied at a time, which saves computational time and space greatly). If flagged true, _pivot_ takes the default value of 0 for each element unless specified different.

## Generating data

To perform the driving frequency analysis itself, call the analyzer's method driving_frequency_analysis. This takes 5 parameters:
1. driving_frequency_range (tuple): Specifies the range of frequencies of the external force. Is in the form (lower_bound, upper_bound)
2. cur_external_force_amplitude (float) [optional]: Amplitude of the driving force in N (the driving torque is then l_1\*F). If unspecified, it takes the default value given in the constructor of the _analyzer_.

## Saving and loading data and using dataset names

## Plotting data
Presets!
