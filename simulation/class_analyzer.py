

import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import statistics

import csv
import io

import cv2
from cv2 import VideoWriter, VideoWriter_fourcc


from class_multipendulum import *

def fig_to_img(fig, dpi=180):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    buf.seek(0)
    img_arr = np.frombuffer(buf.getvalue(), dtype=np.uint8)
    buf.close()
    img = cv2.imdecode(img_arr, 1)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    return img

def subplot_dimensions(number_of_plots):
    # 1x1, 2x1, 3x1, 2x2, 3x2, 4x2, 3x3, 4x3, 5x3
    if number_of_plots == 1:
        return(1, 1)
    if number_of_plots == 2:
        return(2, 1)
    if number_of_plots == 3:
        return(3, 1)
    if number_of_plots == 4:
        return(2, 2)
    if number_of_plots <= 6:
        return(3, 2)
    if number_of_plots <= 8:
        return(4, 2)
    if number_of_plots <= 9:
        return(3, 3)
    if number_of_plots <= 12:
        return(4, 3)
    if number_of_plots <= 15:
        return(5, 3)
    return(np.ceil(np.sqrt(number_of_plots)), np.ceil(np.sqrt(number_of_plots)))

class analyzer:
    
    # ---------- constructors, destructors, descriptors ----------
    
    def __init__(self, my_dt, my_omega_F, my_external_force_amplitude, dataset_name = ''):
        
        self.pendulums = []
        self.pendulum_names = []
        self.pendulum_resonance_analysis_data = []
        
        self.max_total_length = 0.0
        
        self.t = 0
        self.dt = my_dt
        self.omega_F = my_omega_F
        self.external_force_amplitude = my_external_force_amplitude
        
        if dataset_name != '':
            dataset_name += '_'
        self.dataset_name = dataset_name
        self.number_of_videos = 0
        #self.trace_frame = np.zeros((self.o_c_h, self.o_c_w, 3), dtype=np.uint8) + 255
    
    
    def add_pendulum(self, pendulum, pendulum_name = "__UNNAMED__"):
        self.pendulums.append(pendulum)
        if pendulum_name == "__UNNAMED__":
            self.pendulum_names.append("pendulum_" + str(len(self.pendulums)))
        else:
            self.pendulum_names.append(pendulum_name)
        if sum(pendulum.l) > self.max_total_length:
            self.max_total_length = sum(pendulum.l)
    
    # ---------- numerical integration method ------------------
    
    def get_force_phasor(self):
        return(np.sin(self.t * self.omega_F))
    
    def get_force(self, t, omega_F, amplitude, decay_rate = -1):
        if decay_rate == -1:
            return(amplitude * np.sin(t * omega_F))
        return(amplitude * (1.0 - np.exp( - decay_rate * t)) * np.sin(t * omega_F))
        
    
    def step(self, cur_external_force_list=0):
        
        # cur external force list is a list of exteral forces at time t, t+dt/2 and t+dt
        if type(cur_external_force_list)!=list:
            cur_external_force_list=[self.external_force_amplitude*np.sin(self.t * self.omega_F), self.external_force_amplitude*np.sin((self.t + self.dt/2.0) * self.omega_F), self.external_force_amplitude*np.sin((self.t+self.dt) * self.omega_F)]
        
        # RK4
        
        #old_t = self.t
        for pendulum in self.pendulums:
            pendulum.step(self.dt, cur_external_force_list)
        
        self.t += self.dt
    
    # ----------- config memory ------------------------------
    
    def save_state_memory(self, state_memory_values_by_pendulum, driving_frequency_range, external_force_amplitude, datapoints, t_terminate):
        state_memory_file = open('data/config/' + self.dataset_name + 'state_memory.txt', 'w')
        
        omega_F_min, omega_F_max = driving_frequency_range
        config_string = '%f,%f,%f,%i,%f\n' % (omega_F_min, omega_F_max, external_force_amplitude, datapoints, t_terminate)
        state_memory_file.write(config_string)
        
        #print(state_memory_values_by_pendulum)
        
        for i in range(datapoints):
            cur_line = ""
            for p_i in range(len(self.pendulums)):
                for s_i in range(self.pendulums[p_i].N):
                    cur_line += (str(state_memory_values_by_pendulum[p_i][i][0][s_i]) + ',')
                cur_line = cur_line[:-1] + ';'
                for s_i in range(self.pendulums[p_i].N):
                    cur_line += (str(state_memory_values_by_pendulum[p_i][i][1][s_i]) + ',')
                cur_line = cur_line[:-1] + '|'
            state_memory_file.write(cur_line[:-1])
            if i != datapoints - 1:
                state_memory_file.write('\n')
        state_memory_file.close()
                
        
    
    
    def read_state_memory(self):
        print("Loading previous states for the dataset " + self.dataset_name)
        try:
            state_memory_file = open('data/config/' + self.dataset_name + 'state_memory.txt', 'r')
            state_memory_lines = state_memory_file.readlines()
            state_memory_file.close()
            #self.state_memory
            config_line = state_memory_lines[0]
            config_vals = config_line.split(",")
            omega_F_min = float(config_vals[0])
            omega_F_max = float(config_vals[1])
            external_force_amplitude = float(config_vals[2])
            datapoints  = int(config_vals[3])
            t_start     = float(config_vals[4])
            state_memory_lines = state_memory_lines[1:]
            
            print("  Dataset loaded with the following configuration:")
            print("  driving frequency range = (%0.2f, %0.2f), force amplitude = %0.2f, number of datapoints = %i, starting time = %0.2f" % (omega_F_min, omega_F_max, external_force_amplitude, datapoints, t_start))
            
            state_memory_values_by_pendulum = empty_list_list(len(self.pendulums))
            for i in range(datapoints):
                cur_vals_by_pendulum = state_memory_lines[i].split('|')
                for p_i in range(len(self.pendulums)):
                    cur_state = cur_vals_by_pendulum[p_i].split(';')
                    #print(cur_state)
                    cur_theta_state     = elementwise_casting(cur_state[0].split(','), float)
                    cur_theta_dot_state = elementwise_casting(cur_state[1].split(','), float)
                    state_memory_values_by_pendulum[p_i].append([cur_theta_state, cur_theta_dot_state])
            
            return(state_memory_values_by_pendulum, (omega_F_min, omega_F_max), external_force_amplitude, datapoints, t_start)
            
                
        except FileNotFoundError:
            print("  This dataset wasn't analyzed before.")
            return(0)
    
    
    # ----------- simple output methods ----------------------
    
    def initialize_video(self):
        #OpenCV visual output luggage
        self.h, self.w = 720, 1280
        self.FPS = 60
        self.fourcc = VideoWriter_fourcc(*'MP42')
        self.number_of_videos += 1
        #self.number_of_videos = 1
        self.video = VideoWriter('./' + self.dataset_name + '_video_' + str(self.number_of_videos) + '.avi', self.fourcc, float(self.FPS), (self.w, self.h))
        
        self.output_frames = []
        self.layout_frame = np.zeros((self.h, self.w, 3), dtype=np.uint8) + 255
    
    def release_video(self):
        print("  Number of frames =", len(self.output_frames))
        for frame in self.output_frames:
            #print(len(frame), len(frame[0]), len(frame[0][0]))
            self.video.write(frame)
        self.video.release()
        
        #cv2.destroyAllWindows()
        
        self.output_frames = []
        self.trace_frame = np.zeros((self.h, self.w, 3), dtype=np.uint8) + 255
    
    def animate(self, max_t, threshold_t = 0.0):
        
        print("Simulating the pendulums...")
        self.t = 0.0
        t_list = []
        theta_list = []
        position_list = []
        
        progress = 0
        
        t_list.append(self.t)
        theta_list.append(self.pendulums[0].theta[0])
        position_list.append([])
        for pendulum in self.pendulums:
            position_list[-1].append([])
            for segment_i in range(pendulum.N):
                position_list[-1][-1].append(pendulum.theta[segment_i])
        
        while(self.t < max_t):
            self.step()
            if self.t > threshold_t:
                t_list.append(self.t)
                theta_list.append(self.pendulums[0].theta[0])
                position_list.append([])
                
            for pendulum in self.pendulums:
                position_list[-1].append([])
                for segment_i in range(pendulum.N):
                    position_list[-1][-1].append(pendulum.theta[segment_i])
            
            if np.floor(self.t / max_t * 100) > progress:
                progress = np.floor(self.t / max_t * 100)
                print("  %i percent done" % progress, end='\r')
        
        print()
        # Video generation
        print("Generating the video output...")
        self.initialize_video()
        #scale = 100.0 # this will be a separate method
        scale = self.h / (2.0 * self.max_total_length)
        font = cv2.FONT_HERSHEY_SIMPLEX
        
        offset_list = np.linspace(self.w / (len(self.pendulums) + 1), self.w - self.w / (len(self.pendulums) + 1), len(self.pendulums), dtype=int)
        
        progress = 0
        for t_i in range(len(t_list)):
            if np.floor(t_i / len(t_list) * 100) > progress:
                progress = np.floor(self.t / max_t * 100)
                print("  %i percent done" % progress, end='\r')
            # create frame
            #print("Frame", t_i)
            cur_frame = self.layout_frame.copy()
            
            
            for p_i in range(len(self.pendulums)):
                # create offsets for each pendulum
                #print("Pendulum", p_i)
                offset_x = offset_list[p_i]
                offset_y = self.h / 4
                cur_x, cur_y = offset_x, offset_y
                for s_i in range(self.pendulums[p_i].N):
                    new_x = cur_x + self.pendulums[p_i].l[s_i] * scale * np.sin(position_list[t_i][p_i][s_i])
                    new_y = cur_y + self.pendulums[p_i].l[s_i] * scale * np.cos(position_list[t_i][p_i][s_i])
                    
                    cv2.line(cur_frame, (int(cur_x), int(cur_y)), (int(new_x), int(new_y)), (0,0,0), 2)
                    cv2.circle(cur_frame, (int(new_x), int(new_y)), int(scale * self.pendulums[p_i].m_r[s_i]), (0,0,0), 2 )
                    cv2.circle(cur_frame, (int(new_x), int(new_y)), int(scale * self.pendulums[p_i].m_r[s_i]), (0,0,255), -1 )
                    
                    cur_x = new_x
                    cur_y = new_y
            
            self.output_frames.append(cur_frame)
        print()
        print("Releasing video...")
        self.release_video()
        print("Video exported.")
        
        
        #plt.plot(t_list, theta_list)
        #plt.show()
    
    # ----------- comparative analysis methods ---------------
    
    def find_resonant_frequencies(self, use_max_theta = False, min_prominence=5.0):
        
        # by default we use avg_E for this thing
        self.resonant_frequencies = []
        for p_i in range(len(self.pendulums)):
            #self.resonant_frequencies.append([])
            peaks, properties = find_peaks(self.pendulum_resonance_analysis_data[p_i]['avg_E'], prominence=(min_prominence, None))
            self.resonant_frequencies.append(self.pendulum_resonance_analysis_data[p_i]['omega_F'][peaks])
            
    
    
    def driving_frequency_analysis(self, driving_frequency_range, cur_external_force_amplitude = -1, datapoints = 200, t_max = 50.0, t_threshold=10.0, overwrite_stored_states = True, starting_states=[], t_start=0.0):
        
        self.last_driving_frequency_range = driving_frequency_range
        if cur_external_force_amplitude == -1:
            cur_external_force_amplitude = self.external_force_amplitude
        
        # describe all pendulums in t
        print("F = " + str(cur_external_force_amplitude) + "; omega_F_range = " + str(driving_frequency_range))
        for i in range(len(self.pendulums)):
            print("Pendulum: " + self.pendulum_names[i])
            # Physical configuration
            cur_descriptor = self.pendulums[i].parameter_description()
            for line in cur_descriptor:
                print("  " + line)
            # Normal modes
            self.pendulums[i].get_normal_modes()
            line = "  normal mode frequencies [rad/s]: "
            for freq in self.pendulums[i].normal_mode_frequencies:
                line += f"{freq:.3f}, "
            print(line[:-2])
        
        print("-------------------------------------")
        
        omega_F_min, omega_F_max = driving_frequency_range
        omega_space = np.linspace(max(omega_F_min, 0.0), omega_F_max, datapoints)
        
        # analyze each pendulum separately, of course
        # The results will be stored in a matrix rather than a csv, so we can analyze them directly here.
        # Of course, they can be exported by using export_data(pendulum_name=False)
        # First purge the old data
        self.pendulum_resonance_analysis_data = []
        # The data is a dictionary where each value is a list of measured values, for easier data interpretation
        
        # If saving states, create an array for them
        if overwrite_stored_states:
            state_memory_values_by_pendulum = []
        
        for p_i in range(len(self.pendulums)):
            print("Analyzing pendulum " + self.pendulum_names[p_i])
            # initialize relevant data containers
            self.pendulum_resonance_analysis_data.append({})
            cur_data = self.pendulum_resonance_analysis_data[-1]
            cur_data['omega_F'     ] = []
            cur_data['max_theta_1' ] = []
            cur_data['avg_E'       ] = []
            cur_data['avg_L'       ] = []
            cur_data['ang_f_theta' ] = empty_list_list(self.pendulums[p_i].N)
            cur_data['ang_f_phi'   ] = empty_list_list(self.pendulums[p_i].N)
            cur_data['hp_theta_std'] = empty_list_list(self.pendulums[p_i].N)
            cur_data['hp_phi_std'  ] = empty_list_list(self.pendulums[p_i].N)
            
            if overwrite_stored_states:
                state_memory_values_by_pendulum.append([])
            
            start_time = time.time()
            progress = 0
            for omega_i in range(len(omega_space)):
                
                omega_val = omega_space[omega_i]
                
                if np.floor(omega_i / datapoints * 100) > progress:
                    progress = np.floor(omega_i / datapoints * 100)
                    print(  "  Analysis in progress: " + str(progress) + "%; est. time of finish: " + time.strftime("%H:%M:%S", time.localtime( (time.time()-start_time) * 100 / progress + start_time )), end='\r')
                
                self.t = t_start
                if starting_states == []:
                    self.pendulums[p_i].reset_state()
                else:
                    #print(len(starting_states[omega_i][0]), len(starting_states[omega_i][1]))
                    self.pendulums[p_i].set_state(starting_states[p_i][omega_i][0],starting_states[p_i][omega_i][1])
                
                max_theta_1 = 0.0
                avg_E = 0.0
                avg_L = 0.0
                total_halfperiods_theta     = np.zeros(self.pendulums[p_i].N)
                halfperiod_length_std_theta = np.zeros(self.pendulums[p_i].N)
                total_halfperiods_phi       = np.zeros(self.pendulums[p_i].N)
                halfperiod_length_std_phi   = np.zeros(self.pendulums[p_i].N)
                
                # Helping properties

                last_halfperiod_theta_t      = np.zeros(self.pendulums[p_i].N)
                last_halfperiod_phi_t        = np.zeros(self.pendulums[p_i].N)
                halfperiod_theta_length_list = empty_list_list(self.pendulums[p_i].N)
                halfperiod_phi_length_list   = empty_list_list(self.pendulums[p_i].N)
                
                # When to start counting the periods
                period_signs_theta = np.zeros(self.pendulums[p_i].N)
                period_signs_phi   = np.zeros(self.pendulums[p_i].N)
                period_signs_init = False
                t_period_threshold_theta = np.zeros(self.pendulums[p_i].N) + t_threshold#t_threshold
                t_period_threshold_phi   = np.zeros(self.pendulums[p_i].N) + t_threshold#t_threshold
                
                number_of_active_cycles = 0
                
                decay_rate = -1#1.0/30.0
                
                while self.t < t_max:
                    #self.pendulums[p_i].step(self.dt, [cur_external_force_amplitude*np.sin(self.t * omega_val), cur_external_force_amplitude*np.sin((self.t + self.dt/2.0) * omega_val), cur_external_force_amplitude*np.sin((self.t+self.dt) * omega_val)])
                    self.pendulums[p_i].step(self.dt, [self.get_force(self.t, omega_val, cur_external_force_amplitude, decay_rate), self.get_force((self.t + self.dt/2.0), omega_val, cur_external_force_amplitude, decay_rate), self.get_force((self.t + self.dt), omega_val, cur_external_force_amplitude, decay_rate)])
                    self.t += self.dt
                    
                    if self.t > t_threshold:
                        # record data
                        number_of_active_cycles += 1
                        if self.pendulums[p_i].theta[0] > max_theta_1:
                            max_theta_1 = self.pendulums[p_i].theta[0]
                        
                        T, U, E = self.pendulums[p_i].get_total_energy()
                        avg_E += E
                        avg_L += self.pendulums[p_i].get_angular_momentum()
                        
                        self.pendulums[p_i].get_phi()
                        if not period_signs_init:
                            period_signs_theta = nonzero_sign(base_angle(self.pendulums[p_i].theta))
                            period_signs_phi   = nonzero_sign(base_angle(self.pendulums[p_i].phi  ))
                            period_signs_init = True
                        
                        for s_i in range(self.pendulums[p_i].N):
                            if (period_signs_theta[s_i] * nonzero_sign(base_angle(self.pendulums[p_i].theta[s_i])) == -1):
                                period_signs_theta[s_i] *= -1
                                total_halfperiods_theta[s_i] += 1
                                if last_halfperiod_theta_t[s_i] != 0.0:
                                    halfperiod_theta_length_list[s_i].append(self.t - last_halfperiod_theta_t[s_i])
                                else:
                                    t_period_threshold_theta[s_i] = self.t
                                last_halfperiod_theta_t[s_i] = self.t
                            if (period_signs_phi[s_i] * nonzero_sign(base_angle(self.pendulums[p_i].phi[s_i])) == -1):
                                period_signs_phi[s_i] *= -1
                                total_halfperiods_phi[s_i] += 1
                                if last_halfperiod_phi_t[s_i] != 0.0:
                                    halfperiod_phi_length_list[s_i].append(self.t - last_halfperiod_phi_t[s_i])
                                else:
                                    t_period_threshold_phi[s_i] = self.t
                                last_halfperiod_phi_t[s_i] = self.t
                
                #print(total_halfperiods_theta, halfperiod_theta_length_list)
                
                avg_E /= number_of_active_cycles
                avg_L /= number_of_active_cycles
                ang_f_theta = total_halfperiods_theta / (2.0 * (t_max - t_period_threshold_theta)) * 2.0 * np.pi
                ang_f_phi   = total_halfperiods_phi   / (2.0 * (t_max - t_period_threshold_phi  )) * 2.0 * np.pi
                for s_i in range(self.pendulums[p_i].N):
                    if total_halfperiods_theta[s_i] >= 3: #you need 3 halfperiods for two intervals, two intervals to define std
                        halfperiod_length_std_theta[s_i] = statistics.stdev(halfperiod_theta_length_list[s_i])
                    else:
                        halfperiod_length_std_theta[s_i] = 0.0
                    if total_halfperiods_phi[s_i] >= 3:
                        halfperiod_length_std_phi[s_i]   = statistics.stdev(halfperiod_phi_length_list[s_i]  )
                    else:
                        halfperiod_length_std_phi[s_i]   = 0.0
                
                cur_data['omega_F'     ].append(omega_val  )
                cur_data['max_theta_1' ].append(max_theta_1)
                cur_data['avg_E'       ].append(avg_E      )
                cur_data['avg_L'       ].append(avg_L      )
                for s_i in range(self.pendulums[p_i].N):
                    cur_data['ang_f_theta' ][s_i].append(ang_f_theta[s_i]                )
                    cur_data['ang_f_phi'   ][s_i].append(ang_f_phi[s_i]                  )
                    cur_data['hp_theta_std'][s_i].append(halfperiod_length_std_theta[s_i])
                    cur_data['hp_phi_std'  ][s_i].append(halfperiod_length_std_phi[s_i]  )
                
                
                if overwrite_stored_states:
                    state_memory_values_by_pendulum[p_i].append([self.pendulums[p_i].theta, self.pendulums[p_i].theta_dot])
                
            cur_data['omega_F'    ] = np.array(cur_data['omega_F'    ])
            cur_data['max_theta_1'] = np.array(cur_data['max_theta_1'])
            cur_data['avg_E'      ] = np.array(cur_data['avg_E'      ])
            cur_data['avg_L'      ] = np.array(cur_data['avg_L'      ])
            for s_i in range(self.pendulums[p_i].N):
                cur_data['ang_f_theta' ][s_i] = np.array(cur_data['ang_f_theta' ][s_i])
                cur_data['ang_f_phi'   ][s_i] = np.array(cur_data['ang_f_phi'   ][s_i])
                cur_data['hp_theta_std'][s_i] = np.array(cur_data['hp_theta_std'][s_i])
                cur_data['hp_phi_std'  ][s_i] = np.array(cur_data['hp_phi_std'  ][s_i])
            print("  Pendulum " + self.pendulum_names[p_i] + " analysis finished at " + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

        print("Data collection done. Finding resonant frequencies...")
        self.find_resonant_frequencies()
        print("Resonant frequencies found:")
        for p_i in range(len(self.pendulums)):
            if len(self.resonant_frequencies[p_i]) > 0:
                print("  " + self.pendulum_names[p_i] + ': ' + list_to_str(self.resonant_frequencies[p_i]))
            else:
                print("  " + self.pendulum_names[p_i] + ': No resonant frequencies found')
        
        if overwrite_stored_states:
            self.save_state_memory(state_memory_values_by_pendulum, driving_frequency_range, cur_external_force_amplitude, datapoints, t_max)
    
    
    def load_state(self, t_max, t_threshold, overwrite_stored_states = True):
        
        state_memory_values_by_pendulum, driving_frequency_range, external_force_amplitude, datapoints, t_start = self.read_state_memory()
        self.driving_frequency_analysis(driving_frequency_range, external_force_amplitude, datapoints, t_max, t_threshold, overwrite_stored_states, state_memory_values_by_pendulum, t_start)
        
    
    # ---------------------------- analysis output methods --------------------------
        
    def plot_scalar_list(self, myplt, data, linestyle='dotted', label=False, direction='x', color = False):
        if color == False:
            # use the color of the last line, or default to blue
            if len(myplt.gca().lines) == 0:
                color = '#1f77b4'
            else:
                color = myplt.gca().lines[-1].get_color()
        if label != False:
            described = False
        for datapoint in data:
            if label != False:
                if described == False:
                    described = True
                    if 'x' in direction:
                        myplt.axvline(x=datapoint, linestyle=linestyle, label = label, color=color)
                        if 'y' in direction:
                            myplt.axhline(y=datapoint, linestyle=linestyle, color=color)
                    elif 'y' in direction:
                        myplt.axhline(y=datapoint, linestyle=linestyle, label = label, color=color)
                else:
                    if 'x' in direction:
                        myplt.axvline(x=datapoint, linestyle=linestyle, color=color)
                    if 'y' in direction:
                        myplt.axhline(y=datapoint, linestyle=linestyle, color=color)
            else:
                if 'x' in direction:
                    myplt.axvline(x=datapoint, linestyle=linestyle, color=color)
                if 'y' in direction:
                    myplt.axhline(y=datapoint, linestyle=linestyle, color=color)
        
    
    
    def plot_resonance_analysis_data(self, graph_list=['max_theta_1_graph', 'avg_E_graph', 'avg_L_graph', 'ang_f_theta_graph', 'hp_theta_std_graph', 'ang_f_phi_graph', 'hp_phi_std_graph'], names=-1, max_theta_simple_comparison = True, save_graph=False):
        
        if names == -1:
            names = self.pendulum_names
        
        p_indexes = []
        for name in names:
            p_i = self.pendulum_names.index(name)
            p_indexes.append(p_i)
        
        # create a superplot
        superplot_shape_x, superplot_shape_y = subplot_dimensions(len(graph_list))
        
        plt.figure(figsize=(15, 8))
        
        # initialize lims; if changed, apply them afterwards
        x_left = -1
        x_right = -1
        y_left = -1
        y_right = -1
        
        for i in range(len(graph_list)):
            plt.subplot(superplot_shape_x, superplot_shape_y, i + 1)
            if graph_list[i] == 'max_theta_1_graph':
                plt.title("Maximal amplitude")
                #plt.ylim(0, max_theta_val * 1.1)
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("maximal $\\theta_1$ [rad]")
                
                x_left, x_right = self.last_driving_frequency_range
                for p_i in p_indexes:
                    plt.plot(self.pendulum_resonance_analysis_data[p_i]['omega_F'], self.pendulum_resonance_analysis_data[p_i]['max_theta_1'], '.-', label=self.pendulum_names[p_i] + " data")
                    self.plot_scalar_list(plt, self.resonant_frequencies[p_i], label=f"{self.pendulum_names[p_i]} resonant f")
                    self.plot_scalar_list(plt, self.pendulums[p_i].normal_mode_frequencies, linestyle='solid', label=f"{self.pendulum_names[p_i]} normal mode f")
            if graph_list[i] == 'avg_E_graph':
                plt.title("Time-average mechanical energy")
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("$\langle E\\rangle$ [J]")
                x_left, x_right = self.last_driving_frequency_range
                for p_i in p_indexes:
                    plt.plot(self.pendulum_resonance_analysis_data[p_i]['omega_F'], self.pendulum_resonance_analysis_data[p_i]['avg_E'], '.-', label=self.pendulum_names[p_i] + " data")
                    self.plot_scalar_list(plt, self.resonant_frequencies[p_i], label=f"{self.pendulum_names[p_i]} resonant f")
                    self.plot_scalar_list(plt, self.pendulums[p_i].normal_mode_frequencies, linestyle='solid', label=f"{self.pendulum_names[p_i]} normal mode f")
            if graph_list[i] == 'avg_E_ratio_graph':
                plt.title("Time-average mechanical energy (fraction of max val.)")
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("$\langle E\\rangle /|\langle E\\rangle_{max}|$")
                x_left, x_right = self.last_driving_frequency_range
                for p_i in p_indexes:
                    #max_avg_E_index, max_avg_E_amp = max_item(np.abs(np.array(self.pendulum_resonance_analysis_data[p_i]['avg_E'])))
                    max_avg_E_index, max_avg_E_amp = max_item(np.abs(self.pendulum_resonance_analysis_data[p_i]['avg_E']))
                    plt.plot(self.pendulum_resonance_analysis_data[p_i]['omega_F'], self.pendulum_resonance_analysis_data[p_i]['avg_E'] / max_avg_E_amp, '.-', label=self.pendulum_names[p_i] + " data")
                    self.plot_scalar_list(plt, self.resonant_frequencies[p_i], label=f"{self.pendulum_names[p_i]} resonant f")
                    self.plot_scalar_list(plt, self.pendulums[p_i].normal_mode_frequencies, linestyle='solid', label=f"{self.pendulum_names[p_i]} normal mode f")
            if graph_list[i] == 'avg_L_graph':
                plt.title("Time-average angular momentum")
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("$\langle L_z\\rangle$ [J$\cdot$s]")
                x_left, x_right = self.last_driving_frequency_range
                for p_i in p_indexes:
                    plt.plot(self.pendulum_resonance_analysis_data[p_i]['omega_F'], self.pendulum_resonance_analysis_data[p_i]['avg_L'], '.-', label=self.pendulum_names[p_i] + " data")
                    self.plot_scalar_list(plt, self.resonant_frequencies[p_i], label=f"{self.pendulum_names[p_i]} resonant f")
                    self.plot_scalar_list(plt, self.pendulums[p_i].normal_mode_frequencies, linestyle='solid', label=f"{self.pendulum_names[p_i]} normal mode f")
            if graph_list[i] == 'ang_f_theta_graph':
                plt.title("Angular frequencies")
                #plt.ylim(0, max_ang_f_val * 1.1)
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("ang. freq. of seg. [rad.s$^{-1}$]")
                x_left, x_right = self.last_driving_frequency_range
                y_left, y_right = self.last_driving_frequency_range
                for p_i in p_indexes:
                    item = self.pendulum_resonance_analysis_data[p_i]
                    for s_i in range(self.pendulums[p_i].N):
                        plt.plot(item['omega_F'], item['ang_f_theta'][s_i], '.', label=self.pendulum_names[p_i] + " (segment " + str(s_i) + ")")
                    # add resonant frequencies and normal mode frequencies
                    self.plot_scalar_list(plt, self.resonant_frequencies[p_i], label=f"{self.pendulum_names[p_i]} resonant f", direction = 'xy')
                    self.plot_scalar_list(plt, self.pendulums[p_i].normal_mode_frequencies, linestyle='solid', label=f"{self.pendulum_names[p_i]} normal mode f", direction='xy')
            if graph_list[i] == 'ang_f_phi_graph':
                plt.title("Angular frequencies of angle differences")
                #plt.ylim(0, max_ang_f_val * 1.1)
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("ang. freq. of seg. ($\phi$) [rad.s$^{-1}$]")
                x_left, x_right = self.last_driving_frequency_range
                y_left, y_right = self.last_driving_frequency_range
                for p_i in p_indexes:
                    item = self.pendulum_resonance_analysis_data[p_i]
                    for s_i in range(self.pendulums[p_i].N):
                        plt.plot(item['omega_F'], item['ang_f_phi'][s_i], '.', label=self.pendulum_names[p_i] + " (segment " + str(s_i) + ")")
                    # add resonant frequencies and normal mode frequencies
                    self.plot_scalar_list(plt, self.resonant_frequencies[p_i], label=f"{self.pendulum_names[p_i]} resonant f", direction = 'xy')
                    self.plot_scalar_list(plt, self.pendulums[p_i].normal_mode_frequencies, linestyle='solid', label=f"{self.pendulum_names[p_i]} normal mode f", direction='xy')
            """if graph_list[i] == 'ang_f_theta_space':
                plt.title("Angular frequency space")
                #border_coefficient = 1.05
                #plt.xlim(min_ang_f_1_val - max_ang_f_1_val * (border_coefficient - 1.0), max_ang_f_1_val * border_coefficient)
                #plt.ylim(min_ang_f_2_val - max_ang_f_1_val * (border_coefficient - 1.0), max_ang_f_2_val * border_coefficient)
                plt.xlabel("ang. freq. of 1st pend. [rad.s$^{-1}$]")
                plt.ylabel("ang. freq. of 2nd pend. [rad.s$^{-1}$]")
                ang_f_space = np.linspace(0, max_ang_f_val * border_coefficient, 50)
                plt.plot(ang_f_space, ang_f_space, linestyle='dashed', label = "Equal frequency line")
                for i in range(len(my_filenames)):
                    item = df[my_filenames[i]]
                    plt.plot(item['ang_f_1'], item['ang_f_2'], '-', label=my_filenames[i] + " data")
                    # add resonant frequencies
                    for res_f_i in resonant_frequency_index_list[my_filenames[i]]:
                        plt.plot(item['ang_f_1'][res_f_i],item['ang_f_2'][res_f_i],'x',markersize=15, markeredgewidth=2)
                """
            if graph_list[i] == 'hp_theta_std_graph':
                plt.title("Standard deviation of half-periods")
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("$\sigma$ [dimensionless]")
                x_left, x_right = self.last_driving_frequency_range
                for p_i in p_indexes:
                    item = self.pendulum_resonance_analysis_data[p_i]
                    for s_i in range(self.pendulums[p_i].N):
                        plt.plot(item['omega_F'], item['hp_theta_std'][s_i], '.-', label=self.pendulum_names[p_i] + " (segment " + str(s_i) + ")")
                    # observed resonant frequencies and normal modes
                    self.plot_scalar_list(plt, self.resonant_frequencies[p_i], label=f"{self.pendulum_names[p_i]} resonant f")
                    self.plot_scalar_list(plt, self.pendulums[p_i].normal_mode_frequencies, linestyle='solid', label=f"{self.pendulum_names[p_i]} normal mode f")
            if graph_list[i] == 'hp_phi_std_graph':
                plt.title("Standard deviation of half-periods ($\phi$)")
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("$\sigma$ [dimensionless]")
                x_left, x_right = self.last_driving_frequency_range
                for p_i in p_indexes:
                    item = self.pendulum_resonance_analysis_data[p_i]
                    for s_i in range(self.pendulums[p_i].N):
                        plt.plot(item['omega_F'], item['hp_phi_std'][s_i], '.-', label=self.pendulum_names[p_i] + " (segment " + str(s_i) + ")")
                    # observed resonant frequencies and normal modes
                    self.plot_scalar_list(plt, self.resonant_frequencies[p_i], label=f"{self.pendulum_names[p_i]} resonant f")
                    self.plot_scalar_list(plt, self.pendulums[p_i].normal_mode_frequencies, linestyle='solid', label=f"{self.pendulum_names[p_i]} normal mode f")
                    # predicted resonant frequencies
                    #for pred_res_f in predicted_resonant_frequency_list[my_filenames[i]]:
                    #    plt.axvline(x=pred_res_f, linestyle='dotted', color=max_amp_plotline_color_list[my_filenames[i]], label = my_filenames[i] + " predicted res f")
            if x_left != -1:
                plt.xlim(x_left, x_right)
            if y_left != -1:
                plt.ylim(y_left, y_right)
            plt.legend()

        
        plt.tight_layout()
        
        if save_graph:
            plt.savefig("data/outputs/" + self.dataset_name + "graph_output.png")
        plt.show()
    
    
    def save_resonance_analysis_data(self):
        
        if len(self.pendulum_resonance_analysis_data) == 0:
            print("No data to save.")
            return(0)

        for p_i in range(len(self.pendulums)):
            
            cur_data = self.pendulum_resonance_analysis_data[p_i]
            
            output_file = open("data/" + self.dataset_name + self.pendulum_names[p_i] + '.csv', mode='w')
            output_writer = csv.writer(output_file,  delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            
            data_N = len(cur_data['omega_F'])
            
            """datakeys_raw = list(cur_data.keys())
            output_writer.writerow(datakeys)
            
            data_N = len(cur_data[datakeys[0]])
            # print aggregate data
            for i in range(data_N):
                cur_datarow = []
                for key in datakeys:
                    cur_datarow.append(cur_data[key][i])
                output_writer.writerow( cur_datarow )
            """
            regular_keys = ['omega_F', 'max_theta_1', 'avg_E', 'avg_L']
            ang_f_theta_keys  = []
            ang_f_phi_keys    = []
            hp_theta_std_keys = []
            hp_phi_std_keys   = []
            for i in range(self.pendulums[p_i].N):
                ang_f_theta_keys.append( 'ang_f_theta_seg_'  + str(i))
                ang_f_phi_keys.append(   'ang_f_phi_seg_'    + str(i))
                hp_theta_std_keys.append('hp_theta_std_seg_' + str(i))
                hp_phi_std_keys.append(  'hp_phi_std_seg_'   + str(i))
            output_writer.writerow(regular_keys + ang_f_theta_keys + ang_f_phi_keys + hp_theta_std_keys + hp_phi_std_keys)
            for i in range(data_N):
                cur_datarow = []
                for key in regular_keys:
                    cur_datarow.append(cur_data[key][i])
                for s_i in range(self.pendulums[p_i].N):
                    cur_datarow.append(cur_data['ang_f_theta' ][s_i][i])
                    cur_datarow.append(cur_data['ang_f_phi'   ][s_i][i])
                    cur_datarow.append(cur_data['hp_theta_std'][s_i][i])
                    cur_datarow.append(cur_data['hp_phi_std'  ][s_i][i])
                output_writer.writerow( cur_datarow )
            
            # loop through resonant frequency trajectories and for each print data into a separate file
            
            output_file.close()
        # edit config file - t_max, t_thresh, resonant frequencies maybe?
        
    # load resonance data
    def load_resonance_analysis_data(self, analyze_self = True):
        # before loading, you need to add the pendulums into the analyzer
        # the analyzer will match the pendulums' names with the csv files
        if len(self.pendulums) == 0:
            print("No pendulums to load data for.")
            return(0)
        
        # purge old data
        self.pendulum_resonance_analysis_data = []
        
        for p_i in range(len(self.pendulums)):
            
            self.pendulum_resonance_analysis_data.append({})
            cur_data = self.pendulum_resonance_analysis_data[-1]
            cur_data['omega_F'     ] = []
            cur_data['max_theta_1' ] = []
            cur_data['avg_E'       ] = []
            cur_data['avg_L'       ] = []
            cur_data['ang_f_theta' ] = empty_list_list(self.pendulums[p_i].N)
            cur_data['ang_f_phi'   ] = empty_list_list(self.pendulums[p_i].N)
            cur_data['hp_theta_std'] = empty_list_list(self.pendulums[p_i].N)
            cur_data['hp_phi_std'  ] = empty_list_list(self.pendulums[p_i].N)
            
            input_file = open("data/" + self.dataset_name + self.pendulum_names[p_i] + '.csv', newline='')
            input_reader = csv.reader(input_file, delimiter=',', quotechar='"')
            
            input_rows = list(input_reader)
            
            header_row = input_rows[0]
            segment_N = self.pendulums[p_i].N#0
            """segment_key_string = 'ang_f_theta_seg_0'
            while segment_key_string in header_row:
                segment_N += 1
                segment_key_string = 'ang_f_theta_seg_' + str(segment_N)"""
            
            for i in range(1, len(input_rows)):
                cur_data['omega_F'     ].append(float(input_rows[i][0]))
                cur_data['max_theta_1' ].append(float(input_rows[i][1]))
                cur_data['avg_E'       ].append(float(input_rows[i][2]))
                cur_data['avg_L'       ].append(float(input_rows[i][3]))
                row_index = 4
                for s_i in range(segment_N):
                    cur_data['ang_f_theta' ][s_i].append(float(input_rows[i][row_index +               s_i]))
                    cur_data['ang_f_phi'   ][s_i].append(float(input_rows[i][row_index +   segment_N + s_i]))
                    cur_data['hp_theta_std'][s_i].append(float(input_rows[i][row_index + 2*segment_N + s_i]))
                    cur_data['hp_phi_std'  ][s_i].append(float(input_rows[i][row_index + 3*segment_N + s_i]))
            cur_data['omega_F'    ] = np.array(cur_data['omega_F'    ])
            cur_data['max_theta_1'] = np.array(cur_data['max_theta_1'])
            cur_data['avg_E'      ] = np.array(cur_data['avg_E'      ])
            cur_data['avg_L'      ] = np.array(cur_data['avg_L'      ])
            for s_i in range(segment_N):
                cur_data['ang_f_theta' ][s_i] = np.array(cur_data['ang_f_theta' ][s_i])
                cur_data['ang_f_phi'   ][s_i] = np.array(cur_data['ang_f_phi'   ][s_i])
                cur_data['hp_theta_std'][s_i] = np.array(cur_data['hp_theta_std'][s_i])
                cur_data['hp_phi_std'  ][s_i] = np.array(cur_data['hp_phi_std'  ][s_i])
            
            input_file.close()
        
        omega_F_min = self.pendulum_resonance_analysis_data[0]['omega_F'][0]
        omega_F_max = self.pendulum_resonance_analysis_data[0]['omega_F'][-1]
        for p_i in range(1, len(self.pendulums)):
            cur_omega_F_min = self.pendulum_resonance_analysis_data[p_i]['omega_F'][0]
            cur_omega_F_max = self.pendulum_resonance_analysis_data[p_i]['omega_F'][-1]
            if cur_omega_F_min < omega_F_min:
                omega_F_min = cur_omega_F_min
            if cur_omega_F_max > omega_F_max:
                omega_F_max = cur_omega_F_max
        self.last_driving_frequency_range = (omega_F_min, omega_F_max)
        
        if analyze_self:
            self.analyze_resonance_analysis_data()
            
            
        
    # analyze resonance data (find resonant freq., simulate their trajectories)
    def analyze_resonance_analysis_data(self):
        # Measured resonant frequencies
        print("Finding resonant frequencies...")
        self.find_resonant_frequencies()
        print("Resonant frequencies found:")
        for p_i in range(len(self.pendulums)):
            print("  " + self.pendulum_names[p_i] + ': ' + list_to_str(self.resonant_frequencies[p_i]))
        
        # theoretical normal modes
        print("Theoretical normal modes:")
        for p_i in range(len(self.pendulums)):
            self.pendulums[p_i].get_normal_modes()
            line = f"  {self.pendulum_names[p_i]} normal mode frequencies [rad/s]: "
            for freq in self.pendulums[p_i].normal_mode_frequencies:
                line += f"{freq:.3f}, "
            print(line[:-2])
