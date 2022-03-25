
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

First, you create any amount of analyzers as instances of the _analyzer_ class. The constructor takes 3 arguments:

1. my_dt (float): the step used in numerical integration, of dimension [s]
2. my_omega_F (float): the *default value* of the driving force angular frequency.
3. my_external_force_amplitude (float): the *default value* of the driving force amplitude.

The default values will be used in other methods called for this instance if a different value isn't specified.

### Manual population

Second, you create the pendulums in the experiment as instances of the _multipendulum_ class. The constructor takes 3 arguments:

1. my_l (array): list of lengths of the segments. The constructor will infer the amount of segments N from the length of this list.
2. my_m (array): list of masses of the nodes.
3. my_g (float): the gravitational acceleration acting on the multipendulum.

Now you add the pendulums to the analyzers they belong to, by calling the analyzers' method add_pendulum. This method takes 2 arguments:

1. pendulum (multipendulum): the instance of _multipendulum_ you're adding.
2. pendulum_name (string) [optional]: the name of the pendulum, used to mark it on plots and in filenames of exported data. If unspecified, the name will be "pendulum_#", where # is the current amount of pendulums added to the analyzer, which increments by 1 with each time this method's called.

### Automatic parameter-space population

Alternatively, you can skip the pendulum creation and addition and call the analyzer's method populate_parameter_space. This method creates a population of pendulums whose parameters l, m, g span the ranges provided. This method takes 5 arguments:

1. l_space (array of tuples or floats): list of tuples of length N, where the i-th element is a tuple in the form (start, end, amount[, pivot=0]), where _start_ is the lower boundary on the value of l_i, _end_ is the upper boundary, and _amount_ is the number of values we consider on this interval. If a single number is presented instead of the tuple, then l_i only takes this value and doesn't vary in the population. If pivot_expansion is flagged, then _pivot_ marks the index of the pivot value in the produced linspace.
2. m_space (array of tuples or floats): the same as l_space, but for masses of the nodes, so the i-th element describes the space of m_i.
3. g_space (tuple or float): if a tuple in the form (start, end, amount[, pivot=0]), g will take values from _start_ to _end_ with _amount_ datapoints. If a float, g always takes the value of the float. If pivot_expansion is flagged, _pivot_ marks the index of the pivot value in the produced linspace.
4. pivot_expansion (bool) [optional, False by default]: If false, the population will include all possible configurations given by the parameter spaces. If true, only configurations in which all values except one are equal to their respective pivot values are considered (so only one parameter is being varied at a time, which saves computational time and space greatly). If flagged true, _pivot_ takes the default value of 0 for each element unless specified different.
5. parameter_set_name (string) [optional]: The prefix to the indexing name convention, separated by an underscore. 'p' by default.

This method handles the naming convention in a different matter to default name generation: rather than sequential indexing, it names each pendulum as '[parameter_set_name]\_[parameter indices]', where _parameter_indices_ is a list-turned-string of integers separated by commas, where the i-th integer refers to the index in the i-th dimension of the parameter space, omitting dimensions that only take a single value, indexing from 1. For example, is _l_space_=[(1.5, 2.7, 5), 2.0], _m_space_=[(5.5, 6.5, 3), (1.2, 6.9, 3)], _g_=9.8, then the pendulums will have names in the form 'p_[a],[b],[c]', where _a_ ranges from 1 to 5 and _b_,_c_ range from 1 to 3, describing the values of l_1, m_1, m_2 respectively.

## Generating data

To perform the driving frequency analysis itself, call the analyzer's method driving_frequency_analysis. This takes 5 arguments:

1. driving_frequency_range (tuple): Specifies the range of frequencies of the external force. Is in the form (lower_bound, upper_bound)
2. cur_external_force_amplitude (float) [optional]: Amplitude of the driving force in [N] (the driving torque is then l_1\*F). If unspecified, it takes the default value given in the constructor of the _analyzer_.
3. datapoints (int): Number of datapoints for the analysis. The program will iterate through values of omega_F from a linspace given by _driving_frequency_range_ and this value.
4. t_max (float): The value of simulation time at which it terminates for each datapoint, in [s].
5. t_threshold (float): The value of simulation time at which the aggregate data start being recorded. This is to exclude transient behaviour from the collected data.

This method stores the aggregate data for each pendulum in the list _pendulum_resonance_analysis_data_, and also perform a peak search that finds the local maxima on the average mechanical energy data, storing the corresponding driving frequencies in _resonant_frequencies_.

## Saving and loading data and using dataset names

The data can be saved by calling the _analyzer_'s method save_resonance_analysis_data, which takes 1 argument:

1. dataset_name (string) [optional]: The prefix added to every filename before the pendulum's name, separated by an underscore.
This method creates a .csv file _for every pendulum_ in a 'data' subfolder, with the following naming convention: '/data/[dataset_name]\_[pendulum_name].csv', and stores the aggregate data from _pendulum_resonance_analysis_data_ inside. If _dataset_name_ is unspecified, it will be omitted from the naming convention, so the filenames would be '/data/[pendulum_name].csv'.

This data can be loaded in any of the future runs of the program by calling the method load_resonance_analysis_data, which takes 2 arguments:

1. dataset_name (string) [optional]: The prefix added to every filename before the pendulum's name, separated by an underscore. In general, this should match the _dataset_name_ used to save the data you want to load.
2. analyze_self (bool) [optional]: Whether the program should find resonant frequencies of the loaded data. True by default.

This method should be called **instead of _driving_frequency_analysis_**, and is a direct equivalent of it.

## Plotting data

From the moment _pendulum_resonance_analysis_data_ is initialized, either by calling driving_frequency_analysis or load_resonance_analysis_data, the analyzer can plot the data on a graph with a set of subplots by calling the _analyzer_'s method plot_resonance_analysis_data, which takes 4 arguments:
1. graph_list (list of strings) [optional]: The list of names of graphs which should be included as subplots in the resulting graph. They will appear in the same order as in the list. If unspecified, _graph_list_ takes the value of _default_preset_.
2. names (list of strings) [optional]: Human-readable names of pendulums, ordered in the same way as the pendulums were added to the analyzer. If unspecified, _pendulum_names_ are used instead.
3. simple_comparison (tuple) [optional, (True, False) by default]: If enabled, the method will include harmonic oscillator spectra given by SHOs obtained by simplifying the multipendulums' parameters in relevant graphs ('max_theta_1_graph', 'avg_E_graph'). THe first value in the tuple will enable comparison with a SHO of length l_1 for every multipendulum; the second value in the tuple will enable comparison with a SHO of length l_1+l_2+...+l_N for every multipendulum. If multiple multipendulums share a SHO with the same property, it only gets plotted once, lumping their names together in its label.
4. save_graph (bool) [optional, False by default]: If enabled, the graph will be saved in 'data/[dataset_name].png'.

### Presets

Since _graph_list_ can be quite long and cloggy, you may create and use _presets_ instead. These are just prescribed combinations of graphs saved in lists at the beginning of multipendulum_main.py, which you can use as an argument in plot_resonance_analysis_data. Each preset specializes on a certain subset of relevant physical properties, which determines which graphs it includes. Feel free to add any presets you find yourself using regularly!

### List of graphs

There's a few
