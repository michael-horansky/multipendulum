

import matplotlib.pyplot as plt

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

class analyzer:
    
    # ---------- constructors, destructors, descriptors ----------
    
    number_of_videos = 0
    
    def __init__(self, my_dt, my_omega_F, my_external_force_amplitude):
        
        self.pendulums = []
        
        self.t = 0
        self.dt = my_dt
        self.omega_F = my_omega_F
        self.external_force_amplitude = my_external_force_amplitude
        
        #self.trace_frame = np.zeros((self.o_c_h, self.o_c_w, 3), dtype=np.uint8) + 255
    
    
    def add_pendulum(self, pendulum):
        self.pendulums.append(pendulum)
    
    # ---------- numerical integration method ------------------
    
    def get_force_phasor(self):
        return(np.sin(self.t * self.omega_F))
    
    def step(self, cur_external_force_amplitude):
        
        # RK4
        
        old_t = self.t
        for pendulum in self.pendulums:
            old_theta     = pendulum.theta
            old_theta_dot = pendulum.theta_dot
            self.t        = old_t
            
            
            k1_theta     = []
            k1_theta_dot = []
            acc_1 = pendulum.get_acceleration(cur_external_force_amplitude * self.get_force_phasor())
            for i in range(pendulum.N):
                k1_theta.append(     self.dt * pendulum.theta_dot[i] )
                k1_theta_dot.append( self.dt * acc_1[i]              )
            pendulum.theta     = old_theta     + np.array(k1_theta)      / 2.0
            pendulum.theta_dot = old_theta_dot + np.array(k1_theta_dot)  / 2.0
            self.t             = old_t + self.dt / 2.0
            
            k2_theta     = []
            k2_theta_dot = []
            acc_2 = pendulum.get_acceleration(cur_external_force_amplitude * self.get_force_phasor())
            for i in range(pendulum.N):
                k2_theta.append(     self.dt * pendulum.theta_dot[i] )
                k2_theta_dot.append( self.dt * acc_2[i]              )
            pendulum.theta     = old_theta     + np.array(k2_theta)      / 2.0
            pendulum.theta_dot = old_theta_dot + np.array(k2_theta_dot)  / 2.0
            self.t             = old_t + self.dt / 2.0
            
            k3_theta     = []
            k3_theta_dot = []
            acc_3 = pendulum.get_acceleration(cur_external_force_amplitude * self.get_force_phasor())
            for i in range(pendulum.N):
                k3_theta.append(     self.dt * pendulum.theta_dot[i] )
                k3_theta_dot.append( self.dt * acc_3[i]              )
            pendulum.theta     = old_theta     + np.array(k3_theta)
            pendulum.theta_dot = old_theta_dot + np.array(k3_theta_dot)
            self.t             = old_t + self.dt
            
            k4_theta     = []
            k4_theta_dot = []
            acc_4 = pendulum.get_acceleration(cur_external_force_amplitude * self.get_force_phasor())
            for i in range(pendulum.N):
                k4_theta.append(     self.dt * pendulum.theta_dot[i] )
                k4_theta_dot.append( self.dt * acc_4[i]              )
            
            #print("huehue", (np.array(k1_theta)     + 2.0*np.array(k2_theta)     + 2.0*np.array(k3_theta)     + np.array(k4_theta)    )/6.0)
            
            pendulum.theta     = old_theta     + (np.array(k1_theta)     + 2.0*np.array(k2_theta)     + 2.0*np.array(k3_theta)     + np.array(k4_theta)    )/6.0
            pendulum.theta_dot = old_theta_dot + (np.array(k1_theta_dot) + 2.0*np.array(k2_theta_dot) + 2.0*np.array(k3_theta_dot) + np.array(k4_theta_dot))/6.0
            #self.t = old_t
        self.t = old_t + self.dt
            
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
    
    def animate(self, max_t):
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
            self.step(self.external_force_amplitude, )
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
        scale = 100.0 # this will be a separate method
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
                offset_y = self.h / 3
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
    
    def driving_frequency_analysis(self, driving_frequency_range, cur_external_force_amplitude = -1):
        
        if cur_external_force_amplitude == -1:
            cur_external_force_amplitude = self.external_force_amplitude
        
        print("I AM IMPACT")

