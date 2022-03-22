

import matplotlib.pyplot as plt
from scipy.signal import find_peaks

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
    
    number_of_videos = 0
    
    def __init__(self, my_dt, my_omega_F, my_external_force_amplitude):
        
        self.pendulums = []
        self.pendulum_names = []
        self.pendulum_resonance_analysis_data = []
        
        self.max_total_length = 0.0
        
        self.t = 0
        self.dt = my_dt
        self.omega_F = my_omega_F
        self.external_force_amplitude = my_external_force_amplitude
        
        #self.trace_frame = np.zeros((self.o_c_h, self.o_c_w, 3), dtype=np.uint8) + 255
    
    
    def add_pendulum(self, pendulum, pendulum_name = "__UNNAMED__"):
        self.pendulums.append(pendulum)
        if pendulum_name == "__UNNAMED__":
            self.pendulum_names.append("pendulum_" + str(len(self.pendulums)))
        self.pendulum_names.append(pendulum_name)
        if sum(pendulum.l) > self.max_total_length:
            self.max_total_length = sum(pendulum.l)
    
    # ---------- numerical integration method ------------------
    
    def get_force_phasor(self):
        return(np.sin(self.t * self.omega_F))
    
    def step(self, cur_external_force_list=0):
        
        # cur external force list is a list of exteral forces at time t, t+dt/2 and t+dt
        if type(cur_external_force_list)!=list:
            cur_external_force_list=[self.external_force_amplitude*np.sin(self.t * self.omega_F), self.external_force_amplitude*np.sin((self.t + self.dt/2.0) * self.omega_F), self.external_force_amplitude*np.sin((self.t+self.dt) * self.omega_F)]
        
        # RK4
        
        #old_t = self.t
        for pendulum in self.pendulums:
            pendulum.step(self.dt, cur_external_force_list)
        
        self.t += self.dt    
    # ----------- simple output methods ----------------------
    
    def initialize_video(self):
        #OpenCV visual output luggage
        self.h, self.w = 720, 1280
        self.FPS = 60
        self.fourcc = VideoWriter_fourcc(*'MP42')
        analyzer.number_of_videos += 1
        #self.number_of_videos = 1
        self.video = VideoWriter('./multipendulum_video_output' + str(analyzer.number_of_videos) + '.avi', self.fourcc, float(self.FPS), (self.w, self.h))
        
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
                progress += 1
                print("%i percent done" % progress)
        
        
        # Video generation
        print("Commence video generation...")
        self.initialize_video()
        #scale = 100.0 # this will be a separate method
        scale = self.h / (2.0 * self.max_total_length)
        font = cv2.FONT_HERSHEY_SIMPLEX
        
        offset_list = np.linspace(self.w / (len(self.pendulums) + 1), self.w - self.w / (len(self.pendulums) + 1), len(self.pendulums), dtype=int)
        
        progress = 0
        for t_i in range(len(t_list)):
            if np.floor(t_i / len(t_list) * 100) > progress:
                progress += 1
                print("  %i percent done" % progress)
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
                    
        print("Video release commencing...")
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
            
    
    
    def driving_frequency_analysis(self, driving_frequency_range, cur_external_force_amplitude = -1, datapoints = 200, t_max = 50.0, t_threshold=10.0):
        
        if cur_external_force_amplitude == -1:
            cur_external_force_amplitude = self.external_force_amplitude
        
        
        # describe all pendulums in t
        print("F = " + str(cur_external_force_amplitude) + "; omega_F_range = " + str(driving_frequency_range))
        for i in range(len(self.pendulums)):
            print("Pendulum: " + self.pendulum_names[i])
            cur_descriptor = self.pendulums[i].parameter_description()
            for line in cur_descriptor:
                print("  " + line)
        
        print("-------------------------------------")
        
        omega_F_min, omega_F_max = driving_frequency_range
        omega_space = np.linspace(max(omega_F_min, 0.0), omega_F_max, datapoints)
        
        # analyze each pendulum separately, of course
        # The results will be stored in a matrix rather than a csv, so we can analyze them directly here.
        # Of course, they can be exported by using export_data(pendulum_name=False)
        # First purge the old data
        self.pendulum_resonance_analysis_data = []
        # The data is a dictionary where each value is a list of measured values, for easier data interpretation
        for p_i in range(len(self.pendulums)):
            print("Analyzing pendulum " + self.pendulum_names[p_i])
            # initialize relevant data containers
            self.pendulum_resonance_analysis_data.append({})
            cur_data = self.pendulum_resonance_analysis_data[-1]
            cur_data['omega_F'    ] = []
            cur_data['max_theta_1'] = []
            cur_data['avg_E'      ] = []
            #cur_data['ang_f_1'    ] = []
            #cur_data['ang_f_2'    ] = []
            
            start_time = time.time()
            progress = 0
            for omega_i in range(len(omega_space)):
                
                omega_val = omega_space[omega_i]
                
                if np.floor(omega_i / datapoints * 100) > progress:
                    progress += 1
                    print(  "Analysis in progress: " + str(progress) + "%; expected time of finish: " + time.strftime("%H:%M:%S", time.localtime( (time.time()-start_time) * 100 / progress + start_time )))
                
                self.t = 0.0
                self.pendulums[p_i].reset_state()
                
                max_theta_1 = 0.0
                avg_E = 0.0
                
                number_of_active_cycles = 0
                
                
                while self.t < t_max:
                    self.pendulums[p_i].step(self.dt, [cur_external_force_amplitude*np.sin(self.t * omega_val), cur_external_force_amplitude*np.sin((self.t + self.dt/2.0) * omega_val), cur_external_force_amplitude*np.sin((self.t+self.dt) * omega_val)])
                    self.t += self.dt
                    
                    if self.t > t_threshold:
                        # record data
                        number_of_active_cycles += 1
                        if self.pendulums[p_i].theta[0] > max_theta_1:
                            max_theta_1 = self.pendulums[p_i].theta[0]
                        
                        T, U, E = self.pendulums[p_i].get_total_energy()
                        avg_E += E
                
                avg_E /= number_of_active_cycles
                cur_data['omega_F'    ].append(omega_val  )
                cur_data['max_theta_1'].append(max_theta_1)
                cur_data['avg_E'      ].append(avg_E      )
            cur_data['omega_F'    ] = np.array(cur_data['omega_F'    ])
            cur_data['max_theta_1'] = np.array(cur_data['max_theta_1'])
            cur_data['avg_E'      ] = np.array(cur_data['avg_E'      ])
            print("Pendulum " + self.pendulum_names[p_i] + " analysis finished at " + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

        print("Data collection done. Finding resonant frequencies...")
        self.find_resonant_frequencies()
        print("Resonant frequencies found:")
        for p_i in range(len(self.pendulums)):
            print("  " + self.pendulum_names[p_i] + ': ' + list_to_str(self.resonant_frequencies[p_i]))
    
    
    
    # ---------------------------- analysis output methods --------------------------
        
    
    def plot_resonance_analysis_data(self, graph_list=['max_theta_1_graph', 'avg_E_graph'], names=-1, max_theta_simple_comparison = True):
        
        if names == -1:
            names = self.pendulum_names
        
        p_indexes = []
        for name in names:
            p_i = self.pendulum_names.index(name)
            p_indexes.append(p_i)
        
        # create a superplot
        superplot_shape_x, superplot_shape_y = subplot_dimensions(len(graph_list))
        
        plt.figure(figsize=(15, 8))
        for i in range(len(graph_list)):
            plt.subplot(superplot_shape_x, superplot_shape_y, i + 1)
            if graph_list[i] == 'max_theta_1_graph':
                plt.title("Maximal amplitude")
                #plt.ylim(0, max_theta_val * 1.1)
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("maximal $\\theta_1$ [rad]")
                for p_i in p_indexes:
                    plt.plot(self.pendulum_resonance_analysis_data[p_i]['omega_F'], self.pendulum_resonance_analysis_data[p_i]['max_theta_1'], '.-', label=self.pendulum_names[p_i] + " data")
                plt.legend()
            if graph_list[i] == 'avg_E_graph':
                plt.title("Time-average mechanical energy")
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("$\langle E\\rangle$ [J]")
                for p_i in p_indexes:
                    plt.plot(self.pendulum_resonance_analysis_data[p_i]['omega_F'], self.pendulum_resonance_analysis_data[p_i]['avg_E'], '.-', label=self.pendulum_names[p_i] + " data")
                plt.legend()
            if graph_list[i] == 'avg_E_ratio_graph':
                plt.title("Time-average mechanical energy (fraction of max val.)")
                plt.xlabel("$\omega_F$ [rad.s$^{-1}$]")
                plt.ylabel("$\langle E\\rangle /|\langle E\\rangle_{max}|$")
                for p_i in p_indexes:
                    #max_avg_E_index, max_avg_E_amp = max_item(np.abs(np.array(self.pendulum_resonance_analysis_data[p_i]['avg_E'])))
                    max_avg_E_index, max_avg_E_amp = max_item(np.abs(self.pendulum_resonance_analysis_data[p_i]['avg_E']))
                    plt.plot(self.pendulum_resonance_analysis_data[p_i]['omega_F'], self.pendulum_resonance_analysis_data[p_i]['avg_E'] / max_avg_E_amp, '.-', label=self.pendulum_names[p_i] + " data")
                plt.legend()
        
        plt.tight_layout()
        plt.show()
        
    

