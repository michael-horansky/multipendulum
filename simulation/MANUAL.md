
# Introduction and Structure

This software provides a bulk of functionalities to perform driving frequency analysis on multipendulums. This means we drive the upper segment of a multipendulum with a sine torque of a given frequency and look at aggregate properties of the system after a given amount of time, such as the time-average mechanical energy stored in the system or the maximal amplitude of the first-segment angle. This data can be then analyzed further to extract resonant frequencies and perform further analysis on them, like checking the mode coefficients.

The software provided consists of the following files: multipendulum_main.py, class_analyzer.py, class_multipendulum.py. These have the following purposes:
- **multipendulum_main.py:** This is the 'command file'. It contains instances of available classes and calls all their methods that define the experiment you're conducting. This is the file you need to edit to adjust your experiment, from parameters of the multipendulum to the range of driving frequency or other parameter values you're doing a comparative analysis on. **To use the program, you need to run this file with python3**, no flags necessary. This file shouldn't include any function or class definitions.
- **class_analyzer.py:** This is a class file for the _analyzer_ class. The _analyzer_ is a class described by a single bulk of pendulums and the global parameters independent of any pendulum, such as global time. Its purpose is to conduct comparative analysis on its pendulums, with the primary method to do so being _driving_frequency_analysis_. It can also analyze the output from the analysis to find resonant frequencies, check their mode coefficients, or generate a pendulum population that spans a given parameter space.
- **class_multipendulum.py:** This is a class file for the _multipendulum_ class. An instance of this class is a single multipendulum, described by its static properties (the lengths of its segments, the masses of its nodes, the gravitational acceleration) and its dynamic properties (angles and angular velocities of each segment). This class provides methods to obtain the equations of motion of a multipendulum and to perform numerical integration of these equations using RK4.

These files link to each other as listed above in a descending order. You need to place these files in a single folder to run the program. Then run multipendulum_main.py in the environment of your choice.


# Conducting an Experiment

An 'experiment' is a driving frequency analysis conducted on a specified population of multipendulums on a specified range of frequencies.

## Setting up the configuration
Pendulum names :)
### Default values

## Generating data

## Saving and loading data and using dataset names

## Plotting data
Presets!
