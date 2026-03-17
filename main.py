
import sys
import numpy as np
import sympy as sp 
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt 
from scipy.integrate import odeint
from matplotlib import animation 
from PyQt5.QtWidgets import QApplication, QSizePolicy, QWidget, QMainWindow, QMenu, QVBoxLayout, QSlider, QPushButton, QHBoxLayout, QLabel
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure 

class Pendulum: 

    def __init__(self,m, l, t, theta0, v0, beta): 

        self.m = m
        self.l = l
        self.t = t
        self.beta = beta
        self.condt =  [theta0, v0]

    def solve(self):

        m,g,l,t,beta = sp.symbols(('m','g','l','t','beta')) # create m,g,l,t,beta as sp objects for symbolic operation
        theta = sp.Function('theta')(t) # theta as a function of time, so that sympy can derivate it wrt. t

        dtheta = theta.diff(t)
        ddtheta = dtheta.diff(t)

        x,y = l*sp.sin(theta), -l*sp.cos(theta) # defining x and y in dependency with theta, origin is the Pivot of the pendulum

        T = sp.Rational(1,2) * m * (x.diff(t)**2 + y.diff(t)**2) # E.Kin
        V = m*g*y # E.Pot
        L = T-V # the Lagrangian

        lhs = L.diff(theta)
        rhs = sp.diff(L.diff(dtheta), t)

        # Lagrange equation with Rayleigh dissipation
        eq = rhs- lhs + 2*m*beta*dtheta # beta = effective damping 


        eq = sp.solve(eq, ddtheta)[0] # solve the equation for Θ''(acceleration)

        dthetadt_num = sp.lambdify(dtheta, dtheta, modules='numpy') # change the symbolic dtheta function into a callable python function that can be later numerically evaluated
        domegadt_num = sp.lambdify((g,l,beta,theta, dtheta), eq, modules='numpy') # takes g, l, beta and dtheta as inputs and evaualtes eq numerically, which evaluates the angular acceleration 
        x_num,y_num = sp.lambdify((l,theta), x), sp.lambdify((l,theta), y, modules='numpy') # turns x and y into callable Python functions of (l,theta)

        del m,g,l,t,beta # delete the sp objects created earlier to prevent confusion

        g = 9.81
        l = self.l # pendulum length
        t = self.t # time array
        beta = self.beta
        condt = self.condt # initial conditions[theta0, v0]

        def dXdt(X,t,g,l,beta): 
            theta_num, omega_num = X # unpack current state
            return [dthetadt_num(omega_num), # dtheta/dt = angular velocity(omega)
                    domegadt_num(g,l,beta,theta_num, omega_num)] # domega/dt = angular acceleration
        
        sol = odeint(dXdt,y0 = condt, t=t, args= (g,l,beta)) 
        angle = sol.T[0]
        velocity = sol.T[1]
        acceleration = domegadt_num(g, l, beta, angle, velocity)

        KE = 0.5 * self.m * (self.l**2) * (velocity**2)
        PE = self.m * g * self.l * (1 - np.cos(angle))
        TE = KE + PE
        
        result = {
        'x': x_num(l, angle),
        'y': y_num(l, angle),
        'angle': angle,
        'velocity': velocity,
        'acceleration': acceleration,
        'KE': KE,
        'PE': PE,
        'TE': TE
        }
        return result
    

class MyMplCanvas(FigureCanvas): # inherits from FigureCanvas : Qt Widget to display matplotlib figures
    def __init__(self, parent=None, width=10, height=7, dpi=200):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax1 = fig.add_subplot(221) # pendulum
        self.ax2 = fig.add_subplot(222) # (x,t) graph
        self.ax3 = fig.add_subplot(223) # energy graph
        self.ax4 = fig.add_subplot(224) # phase plane
        fig.subplots_adjust(wspace=0.3, hspace=0.4)
        FigureCanvas.__init__(self, fig) # turns the figure into a Qt Widget
        self.setParent(parent) # embeds this canvas into parent Qt Widget(main_widget)
        FigureCanvas.setSizePolicy(self, # resizes the canvas to fill available horizontal and vertical space
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self) 


class MyStaticMplCanvas(MyMplCanvas):
        def update_figure(self,m,l,theta0,v0,beta, t):
            # clear all axes     
            self.ax1.cla()
            self.ax2.cla()
            self.ax3.cla()
            self.ax4.cla()
            # create time array, solve physics and unpack the result
            self.t = np.linspace(0,t,t*60) # 30s animation duration
            pend = Pendulum(m,l,self.t,theta0,v0,beta)
            result = pend.solve() # gets dictionary
            # unpacking the result dict
            for key, value in result.items():
                setattr(self, key, value)


            # Calculate system info
            self.omega0 = np.sqrt(9.81 / l) # natural frequency 
            self.period = 2 * np.pi / self.omega0 # time period of the pendulum
            self.tau = 1 / beta if beta > 0.01 else float('inf') # Time Constant
            self.zeta = beta / self.omega0 # Damping ratio
            
            # Classifying the system
            if self.zeta < 0.1:
                self.damping_type = "Lightly Damped"
            elif self.zeta < 1.0:
                self.damping_type = "Underdamped"
            elif 0.9 <= self.zeta <= 1.1:
                self.damping_type = "Critically Damped"
            else:
                self.damping_type = "Overdamped"

            # Defining the Stop threshold
            initial_amplitude = np.max(np.abs(self.x[:100]))
            self.stop_threshold = max(0.005, 0.01 * initial_amplitude) # stops at 1% or 5 mm, whichever is larger

            # Setting Up Plots
            # Pendulum Objects
            self.string, = self.ax1.plot([0,self.x[0]], [0,self.y[0]], lw=1, c='green')
            bob_radius = 0.1
            self.circle = self.ax1.add_patch(plt.Circle((self.x[0],self.y[0]), bob_radius, fc='r', zorder = 3))

            # Pendulum Plot
            self.ax1.axis('off')
            max_length = 2.0
            self.ax1.set_xlim([-max_length - 0.2, max_length + 0.2])
            self.ax1.set_ylim([-max_length - 0.2, max_length + 0.2])
            self.ax1.set_aspect('equal')
            self.ax1.grid(True, alpha=0.3) # alpha = grid transparency
            self.ax1.set_xlabel('x (m)')
            self.ax1.set_ylabel('y (m)')
            self.ax1.set_title('Pendulum Animation')

            # x,t graph
            self.line_x, = self.ax2.plot([],[], lw=1, c = 'blue', label = 'x(t)')
            self.ax2.set_xlim([self.t[0], self.t[-1]])
            self.ax2.set_ylim([-max_length - 0.2, max_length + 0.2])
            self.ax2.set_xlabel('Time (s)')
            self.ax2.set_ylabel('Displacement x (m)')
            self.ax2.set_title('Displacement vs Time')
            self.ax2.grid(True, alpha=0.3)
            self.ax2.legend()

            # energy graph
            self.line_KE, = self.ax3.plot([], [], 'b-', label='Kinetic', lw=1, alpha=0.8)
            self.line_PE, = self.ax3.plot([], [], 'r-', label='Potential', lw=1, alpha=0.8)
            self.line_TE, = self.ax3.plot([], [], 'k--', label='Total', lw=1)
            self.ax3.set_xlim([self.t[0], self.t[-1]])
            self.ax3.set_ylim([0, np.max(self.TE) * 1.1])
            self.ax3.set_xlabel('Time (s)')
            self.ax3.set_ylabel('Energy (J)')
            self.ax3.set_title('Energy vs Time')
            self.ax3.grid(True, alpha=0.3)
            self.ax3.legend(loc='upper right')
            
            # Phase Plot
            self.angle_wrapped, self.velocity_wrapped = self.wrap_phase_data(self.angle, self.velocity) # wrap angles to [-pi, pi] and add NaN in jumps (in case of full rotations)
            self.line_phase, = self.ax4.plot([],[], lw=1, c = 'blue', label = 'ω(θ)')

            # current state marker - red dot
            self.phase_marker, = self.ax4.plot([],[], 'ro', markersize = 1, zorder = 10, label = 'Current State') # zorder is the drawing order
            
            # Vector arrows in Phase space
            # Horizontal Arrow 
            '''
            Length = Magnitude of Angular Velocity
            Points right when angle is increasing 
            Points left when angle is decreasing
            '''
            self.phase_arrow_h = self.ax4.quiver([0], [0], [0], [0],
                                               color ='blue',
                                               scale = 2, # smaller = longer arrows
                                               scale_units = 'xy',
                                               angles='xy',
                                               width = 0.005,
                                               headwidth = 4,
                                               headlength=5,
                                               zorder=9,
                                               label='Velocity'                                   

            )
            
            # Vertical Arrow 
            '''
            Length = Magnitude of Angular Acceleration
            Points up when velocity is increasing 
            Points down when velocity is decreasing
            '''
            self.phase_arrow_v = self.ax4.quiver([0], [0], [0], [0],
                                               color ='green',
                                               scale = 2, 
                                               scale_units = 'xy',
                                               angles='xy',
                                               width = 0.005,
                                               headwidth = 4,
                                               headlength=5,
                                               zorder=9,
                                               label='Acceleration'                                   

            )

            # Equilibrium point
            self.ax4.plot(0,0, 'go', markersize = 1, label = 'Equilibrium', zorder=5)
            # Reference lines
            self.ax4.axhline(0,color='k',linewidth=0.5, alpha=0.3)
            self.ax4.axvline(0,color='k', linewidth=0.5, alpha=0.3)

            theta_max = np.pi * 1.1 # padding
            omega_max = 10 * 1.1 # padding
            self.ax4.set_xlim([-theta_max, theta_max])
            self.ax4.set_ylim([-omega_max, omega_max])
            self.ax4.set_xlabel('Theta(θ)')
            self.ax4.set_ylabel('Angular Velocity(ω)')
            self.ax4.set_title('Phase Space (θ vs ω)')
            self.ax4.grid(True, alpha=0.3)
            self.ax4.legend(loc='upper right')

            self.animated_artists =(
                self.string,
                self.circle,
                self.line_x,
                self.line_KE, 
                self.line_PE,
                self.line_TE,
                self.line_phase,
                self.phase_marker,
                self.phase_arrow_h,
                self.phase_arrow_v
            )

            self.draw()   

        def wrap_phase_data(self, angle, velocity):
            """Wrap angle to [-pi,pi] and insert NaN at discontinuities"""
            # Check if there are any full rotations
            if np.max(angle) - np.min(angle) <= 2*np.pi:
                # No full rotations, returns original values
                return angle, velocity
            # Wrap Angles to [-pi,pi]
            angle_wrapped = np.mod(angle + np.pi, 2*np.pi) -np.pi

            # Find Jumps
            diff = np.diff(angle_wrapped)
            jump_indices = np.where(np.abs(diff) > np.pi)[0] # np.where returns a one element tuple here so [0]

            # Insert NaN at jump points to break the line
            angle_plot = list(angle_wrapped)
            velocity_plot = list(velocity)

            # Insert from end to beginning to preserve indices
            for i in reversed(jump_indices):
                angle_plot.insert(i + 1, np.nan)
                velocity_plot.insert(i + 1, np.nan)
            
            return np.array(angle_plot), np.array(velocity_plot) 


# GUI
class ApplicationWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Pendulum Simulator")
        self.setStyleSheet('QLabel {font-size: 12px; }')
        # menu bar
        self.file_menu = QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.close, Qt.CTRL + Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        # main container
        self.main_widget = QWidget()
        layout = QVBoxLayout(self.main_widget)

        # Info Label - shows calculated system parameters 
        self.info_label = QLabel()
        self.info_label.setStyleSheet('''
            QLabel {
                font-size: 14px;
                font-family: monospace;
                padding: 10px;
                background-color: #f0f0f0;
                border: 1px solid #ccc;
                border-radius: 5px;
            }
        ''')
        self.info_label.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        
        self.description_label = QLabel() # parameter description
        description_text = """
        <b>PARAMETER EXPLANATIONS</b><br>
        <b>Effective Damping β:</b> Rate at which energy is lost due to friction/air resistance. 
        Higher values mean faster decay of oscillations.<br>
        0.01-0.1 s⁻¹ = realistic air damping<br>
        0.5-2.0 s⁻¹ = water-like damping<br>
        <b>Damping Coefficient:</b> Physical damping constant (b = 2mβ) relating 
        damping force to velocity.<br>
        <b>Natural Frequency ω₀:</b> How fast the pendulum would oscillate without damping. 
        Determined by √(g/l).<br>
        <b>Period T₀:</b> Time for one complete oscillation without damping (2π/ω₀).<br>
        <b>Time Constant τ:</b> Time for amplitude to decay to ~37% of initial value. 
        Larger τ = slower decay.<br>
        <b>Damping Ratio ζ:</b> Compares actual damping to critical damping.<br>
        - ζ &lt; 1: Underdamped (oscillates)<br>
        - ζ = 1: Critically damped (fastest return without oscillation)<br>
        - ζ &gt; 1: Overdamped (slow return, no oscillation)<br>
        """
        self.description_label.setText(description_text)
        self.description_label.setStyleSheet('''
            QLabel {
                font-size: 14px;
                font-family: monospace;
                padding: 10px;
                background-color: #f0f0f0;
                border: 1px solid #ccc;
                border-radius: 5px;
            }
        ''')
        self.description_label.setWordWrap(True)
        self.description_label.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.description_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)



        # Sliders
        self.m_slider = QSlider(Qt.Horizontal, minimum=1, maximum=50, value=10)
        self.m_slider.sliderReleased.connect(self.update_plot)

        self.l_slider = QSlider(Qt.Horizontal, minimum=10, maximum=20, singleStep=1, value=10)
        self.l_slider.sliderReleased.connect(self.update_plot)

        self.theta0_slider = QSlider(Qt.Horizontal, minimum=-170, maximum=170, singleStep=1, value=30)
        self.theta0_slider.sliderReleased.connect(self.update_plot)

        self.v0_slider = QSlider(Qt.Horizontal, minimum=-10, maximum=10, singleStep=1, value=0)
        self.v0_slider.sliderReleased.connect(self.update_plot)

        self.beta_slider = QSlider(Qt.Horizontal, minimum=0, maximum=50, singleStep=1, value=5)
        self.beta_slider.sliderReleased.connect(self.update_plot)

        self.t_slider = QSlider(Qt.Horizontal, minimum=10, maximum=120, value=30)
        self.t_slider.sliderReleased.connect(self.update_plot)

        # Buttons
        self.play_pause_button = QPushButton('Start')
        self.play_pause_button.setFixedWidth(150)
        self.play_pause_button.clicked.connect(self.toggle_animation)

        self.reset_button = QPushButton('Reset')
        self.reset_button.setFixedWidth(150)
        self.reset_button.clicked.connect(self.reset_animation)

        self.sc = MyStaticMplCanvas(self.main_widget) # creates the canvas with 4 subplots


        controls_layout = QHBoxLayout()

        sliders_layout = QVBoxLayout()
        sliders_layout.addLayout(self.make_slider_row("Mass (Kg):", self.m_slider, divide_by_10=True))
        sliders_layout.addLayout(self.make_slider_row("Length (m):", self.l_slider, divide_by_10=True))
        sliders_layout.addLayout(self.make_slider_row("Angle (deg):", self.theta0_slider))
        sliders_layout.addLayout(self.make_slider_row("Angular Velocity (rad/s):", self.v0_slider))
        sliders_layout.addLayout(self.make_slider_row("Damping Parameter β (1/s) :", self.beta_slider, divide_by_10=True))
        sliders_layout.addLayout(self.make_slider_row("Simulation Duration (s):", self.t_slider))


        # Layout assembly - 1,0,2 are the strech factors
        controls_layout.addLayout(sliders_layout, 1) 
        controls_layout.addWidget(self.info_label, 0)
        controls_layout.addWidget(self.description_label, 2)

        # Buttons
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.play_pause_button)
        button_layout.addWidget(self.reset_button)
        button_layout.addStretch()

        # Add to main layout
        layout.addLayout(controls_layout)
        layout.addLayout(button_layout)
        layout.addWidget(self.sc)

        self.setCentralWidget(self.main_widget) # set as main content
        self.update_plot() # initializes the plot with default values

    def stop_animation(self):
        """Safely stop any running animation"""
        if hasattr(self, 'ani') and self.ani is not None:
            try:
                self.ani.event_source.stop()
            except (AttributeError, RuntimeError):
                pass
            self.ani = None

    def update_plot(self):
        # Stop animation safely
        self.stop_animation()
        
        # Reset button state
        self.play_pause_button.setText('Start')
        self.is_paused = False
        
        # Update plot

        # Get UI Values
        m = self.m_slider.value() / 10
        l = self.l_slider.value() / 10
        theta0 = np.radians(self.theta0_slider.value())
        v0 = self.v0_slider.value()
        beta = self.beta_slider.value() / 10
        t = self.t_slider.value()

        # Update Canvas
        self.sc.update_figure(m, l, theta0, v0, beta, t)

        # Calculate Info label
        b_linear = 2 * beta * m 
    
        tau_str = f"{self.sc.tau:.1f} s" if self.sc.tau != float('inf') else "∞"
        info_text = f"""<b>SYSTEM PARAMETERS</b><br>
    Mass m: {m:.2f} kg<br>
    Length l: {l:.2f} m<br>
    Effective Damping β: {beta:.2f} s⁻¹<br>
    Damping Coefficient: {b_linear:.3f} kg/s<br>
    <br>
    <b>CALCULATED VALUES</b><br>
    Natural freq ω₀: {self.sc.omega0:.2f} rad/s<br>
    Period T₀: {self.sc.period:.2f} s<br>
    Time constant τ: {tau_str}<br>
    Damping ratio ζ: {self.sc.zeta:.3f}<br>
    Type: <b>{self.sc.damping_type}</b>"""
        

        self.info_label.setText(info_text)
        self.info_label.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Preferred)
        self.info_label.adjustSize()

    def start_animation(self):
        # stop any existing animation
        self.stop_animation()

        check_window = 100 # Considers the last 100 frames to see if the oscillation has decayed enough

        def animate(i):
            # update pendulum
            self.sc.string.set_data([0,self.sc.x[i]], [0,self.sc.y[i]])
            self.sc.circle.set_center((self.sc.x[i], self.sc.y[i]))
            
            # update (x,t) line
            self.sc.line_x.set_data(self.sc.t[:i], self.sc.x[:i])
            
            # update Energy plot
            self.sc.line_KE.set_data(self.sc.t[:i], self.sc.KE[:i])
            self.sc.line_PE.set_data(self.sc.t[:i], self.sc.PE[:i])
            self.sc.line_TE.set_data(self.sc.t[:i], self.sc.TE[:i])

            # update phase plot
            # local state variables at time i
            theta = self.sc.angle_wrapped[i]
            omega = self.sc.velocity_wrapped[i]
            alpha = self.sc.acceleration[i]

            self.sc.line_phase.set_data(self.sc.angle_wrapped[:i], self.sc.velocity_wrapped[:i]) # draws phase line till ith time, and uses wrapped values with NaN to
                                                                                                 # draw discontinous open wavy lines in case of full rotations
            self.sc.phase_marker.set_data([theta], [omega]) # moves the phase marker to the new position
            self.sc.phase_arrow_h.set_offsets([[theta, omega]]) # moves arrow base to the new position
            self.sc.phase_arrow_h.set_UVC([omega], [0]) # velocity vector
            self.sc.phase_arrow_v.set_offsets([[theta,omega]]) 
            self.sc.phase_arrow_v.set_UVC([0],[alpha]) # acceleration   

            if i > check_window and i % 10 == 0: # check every 10th frame after frame 100
                recent_max = np.max(np.abs(self.sc.x[max(0, i-check_window):i])) # max displacement in last 100 frames

                if recent_max < self.sc.stop_threshold: # if oscillation amplitude is tiny 
                    self.ani.event_source.stop()
                    self.ani = None
                    self.play_pause_button.setText('Start')
                    self.is_paused = False

            return self.sc.animated_artists

        nsteps = len(self.sc.x) # total data points
        nframes = nsteps # one frame per data point
        dt = self.sc.t[1] - self.sc.t[0] # time step between frames
        interval = dt * 1000 # convert to miliseconds for plot
        self.ani = animation.FuncAnimation(self.sc.figure, animate, frames=nframes, repeat = False, interval=interval, blit = True)
        self.sc.draw()

    def toggle_animation(self):
        # create an animation if there isnt one already
        if not hasattr(self, 'ani') or self.ani is None:
            self.start_animation()
            self.play_pause_button.setText('Pause') # changes the button to pause when the animation is running
            self.is_paused = False
        
        else :
            try:
                # Toggle Play/Pause
                if self.is_paused: 
                    self.ani.event_source.start()
                    self.is_paused = False
                    self.play_pause_button.setText('Pause')
                
                else:
                    # Pause
                    self.ani.event_source.stop()
                    self.is_paused = True
                    self.play_pause_button.setText('Resume')
            
            except (AttributeError, RuntimeError):
                # Animation finished or in bad state, restart
                self.start_animation()
                self.play_pause_button.setText('Pause')
                self.is_paused = False

    def reset_animation(self):
        if hasattr(self, 'ani') and self.ani is not None:
            try : 
                self.ani.event_source.stop()
            except (AttributeError, RuntimeError):
                pass     
            self.ani = None

        # Reset button state
        self.play_pause_button.setText('Start')
        self.is_paused = False

        # reset sliders to default
        self.m_slider.setValue(10)
        self.l_slider.setValue(10)
        self.theta0_slider.setValue(-30)
        self.v0_slider.setValue(0)
        self.beta_slider.setValue(5)
        self.t_slider.setValue(30)
        self.update_plot()  # redraws with default values

    def make_slider_row(self, label_text, slider, divide_by_10=False):
        row = QHBoxLayout()
        label = QLabel(label_text)
        label.setFixedWidth(150)

        # slider.setMaximumWidth(300)

        if divide_by_10:
            value_label = QLabel(str(slider.value()/10))
            slider.valueChanged.connect(lambda v: value_label.setText(str(v/10)))

        else:
            value_label = QLabel(str(slider.value()))
            slider.valueChanged.connect(lambda v: value_label.setText(str(v)))
            
        value_label.setFixedWidth(50)
        row.addWidget(label)
        row.addWidget(value_label)
        row.addWidget(slider)
        return row

# setting global matplotlib parameters
params = {'axes.labelsize': 6,
        'axes.titlesize': 8, 
        'xtick.labelsize': 4,
        'ytick.labelsize': 4,
        'legend.fontsize': 4}
plt.rcParams.update(params)

app = QApplication(sys.argv) # creates Qt Application Instance
win = ApplicationWindow() # Sets up the UI, sliders, and canvas
win.show()
sys.exit(app.exec_()) # event loop with app.exec_(), which gets terminated when the program finishes and window is closed