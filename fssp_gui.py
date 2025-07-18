import tkinter as tk
from tkinter import ttk, messagebox
from tkinter.filedialog import asksaveasfilename
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv

# -=/O/=- JimTheoch -=/O/=-

# ==============================================================================
# SECTION 1: DEFINITION OF ALGORITHM STRATEGIES
# ==============================================================================
class FSAlgorithm:
    # Abstract Base Class for any FSSP algorithm.
    # Defines what every algorithm we create must provide.
    def get_time_complexity(self, n):
        return (2 * n) + 2 # Default value, can be overridden
    def get_states(self): raise NotImplementedError
    def get_rules(self): raise NotImplementedError
    def get_colors(self): raise NotImplementedError
    def get_display_names(self): return {v: k for k, v in self.get_states().items()}

class MazoyerOptimalAlgorithm(FSAlgorithm):
    # Class that implements the optimal solution by Mazoyer with time T=2n-2.
    def get_states(self): 
        # Defines all possible states of the automaton.
        return {'Q': 0, 'G': 1, 'A': 2, 'B': 3, 'B1': 4, 'B2': 5, 'F': 6, 'BDR': 7}
    
    def get_colors(self): 
        # Defines the colors for each state for the plot.
        return {0: 'white', 1: 'blue', 2: 'cyan', 3: 'lime', 4: 'green', 5: 'yellowgreen', 6: 'red', 7: 'black'}
    
    def get_display_names(self):
        # Names to be displayed in the plot's colorbar.
        return {
            0: 'Q (Quiescent)', 1: 'G (General)', 2: 'A (Fast Signal)', 3: 'B (Slow Signal)',
            4: 'B1', 5: 'B2', 6: 'F (Fire!)', 7: 'BDR (Boundary)'
        }
        
    def get_rules(self):
        # The heart of the algorithm: All transition rules.
        S = self.get_states(); rules = {}
        Q, G, A, B, B1, B2, F, BDR = S['Q'], S['G'], S['A'], S['B'], S['B1'], S['B2'], S['F'], S['BDR']
        
        # Propagation rules for the fast waves
        rules[(BDR, G, Q)] = A
        rules[(G, Q, Q)] = A
        rules[(A, Q, Q)] = A
        
        # Reflection rule on the right boundary and creation of the slow wave 'B'.
        rules[(A, Q, BDR)] = B
        # The slow wave 'B' uses 3 states (B, B1, B2) to move at a 1/3 speed.
        for l_state in [Q, BDR, G, A]:
            rules[(l_state, B, Q)] = B1
            rules[(l_state, B1, Q)] = B2
            rules[(l_state, B2, Q)] = B

        # ---- The optimization for 2n-2 time ----
        # Direct collision rules: When A sees B next to it, it becomes F.
        # This saves the "lost" step.
        rules[(A, A, B)] = F
        rules[(A, A, B1)] = F
        rules[(A, A, B2)] = F
        
        # Rules for creating a new General 'G' when there is an intermediate 'Q'.
        rules[(A, Q, B)] = G
        rules[(A, Q, B1)] = G
        rules[(A, Q, B2)] = G
        
        # Firing rules for when a 'G' gets trapped.
        rules[(A, G, B)] = F
        rules[(A, G, B1)] = F
        rules[(A, G, B2)] = F

        # Special rules for small N (e.g., N=2)
        rules[(G, Q, BDR)] = F
        rules[(A, A, BDR)] = F
        
        # Rule that makes the 'F' (Fire) state propagate to everyone.
        for s1 in S.values():
            for s2 in S.values():
                if s1 != F: rules[(F, s1, s2)] = F
                if s2 != F: rules[(s1, s2, F)] = F
        return rules

# ==============================================================================
# SECTION 2: The Central Automaton Engine
# ==============================================================================
class FiringSquadAutomaton:
    # This class is the "engine" that runs the cellular automaton.
    # It accepts a "strategy" (an algorithm) as an argument.
    def __init__(self, n, algorithm_strategy: FSAlgorithm): 
        self.n = n # Number of cells
        self.algorithm = algorithm_strategy # The algorithm to be used
        self.states = self.algorithm.get_states() # Gets states from the algorithm
        self.state_names = {v: k for k, v in self.states.items()} # Reverse dictionary
        self.rules = self.algorithm.get_rules() # Gets rules from the algorithm
        self.reset() # Resets to the initial state

    def reset(self):
        # Reset to the initial state.
        self.time = 0 # Time reset
        self.is_firing = False # No firing has occurred yet
        self.cells = np.full(self.n, self.states['Q'], dtype=int) # Fills the array with 'Q'
        self.cells[0] = self.states['G'] # Places the General at the start
        self.history = [self.cells.copy()] # Stores the first step in the history

    def step(self):
        # Executes one time step.
        if self.is_firing: return # If we have already fired, stop
        
        next_cells = self.cells.copy() # Temporary array for the next step
        bdr_state = self.states['BDR'] # Shortcut for the boundary state
        fire_state = self.states['F'] # Shortcut for the firing state
        
        # Loop through each cell
        for i in range(self.n):
            if self.cells[i] == fire_state: continue # If a cell is F, it remains F
            
            # Find the neighbors. If we are at the edge, the neighbor is BDR.
            left = self.cells[i-1] if i > 0 else bdr_state
            curr = self.cells[i]
            right = self.cells[i+1] if i < self.n - 1 else bdr_state
            
            # Look for a rule for this triplet.
            # If none exists, the cell remains in the same state.
            next_cells[i] = self.rules.get((left, curr, right), curr)
            
        self.cells = next_cells # Update the array
        self.history.append(self.cells.copy()) # Store in history
        self.time += 1 # Increment time
        
        # Check if everyone fired simultaneously
        if np.all(self.cells == fire_state): 
            self.is_firing = True

# ==============================================================================
# SECTION 3: Graphical User Interface (GUI)
# ==============================================================================
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Firing Squad Synchronization Problem")
        self.geometry("1400x800")
        
        # GUI state variables
        self.automaton = None
        self.simulation_running = False
        self.simulation_paused = False
        self.simulation_finished = False
        
        # Parameters for display and zoom
        self.max_display_rows = 25
        self.min_display_rows = 5
        self.zoom_step = 5
        self.current_scroll_pos = 0
        self.steps_per_update = 1
        
        # Create the figure for the Matplotlib plot
        self.fig = plt.figure(figsize=(10, 8), dpi=100)
        self.ax = None
        
        # Dictionary that maps algorithm names to their classes
        self.algorithm_map = {"Mazoyer (Optimal)" : MazoyerOptimalAlgorithm}
        
        # Create all graphical elements
        self._configure_styles()
        self._create_widgets()
        self._reset_application_state()
        self.draw_spacetime_diagram()
        self.protocol("WM_DELETE_WINDOW", self._on_closing)

    def _configure_styles(self):
        # Configuration of the style for buttons and other widgets
        s = ttk.Style()
        s.theme_use('clam')
        s.configure('Green.TButton', background='#4CAF50', foreground='black', font=('Arial', 10, 'bold'))
        s.map('Green.TButton', background=[('active', '#8BC34A'), ('disabled', '#A5D6A7')])
        s.configure('Red.TButton', background='#F44336', foreground='black', font=('Arial', 10, 'bold'))
        s.map('Red.TButton', background=[('active', '#EF9A9A'), ('disabled', '#E57373')])

    def _create_widgets(self):
        # Creation of all buttons, entry fields, sliders etc.
        top_control_frame = ttk.Frame(self, padding="10"); top_control_frame.pack(side=tk.TOP, fill=tk.X)
        left_controls = ttk.Frame(top_control_frame); left_controls.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        ttk.Label(left_controls, text="Algorithm:").pack(side=tk.LEFT, padx=5)
        self.algo_combo = ttk.Combobox(left_controls, values=list(self.algorithm_map.keys()), width=25)
        self.algo_combo.pack(side=tk.LEFT, padx=5)
        self.algo_combo.current(0)
        self.algo_combo.config(state='readonly') # Read-only

        ttk.Label(left_controls, text="Soldiers (N):").pack(side=tk.LEFT, padx=(10,0))
        self.n_entry = ttk.Entry(left_controls, width=5)
        self.n_entry.pack(side=tk.LEFT, padx=(0,5))

        ttk.Label(left_controls, text="Steps/Update:").pack(side=tk.LEFT, padx=(5,0))
        self.steps_entry = ttk.Entry(left_controls, width=5)
        self.steps_entry.pack(side=tk.LEFT, padx=(0,10))

        # Initial default values
        self.n_entry.insert(0, "0")
        self.steps_entry.insert(0, "0")

        self.run_button = ttk.Button(left_controls, text="RUN", command=self._handle_run)
        self.run_button.pack(side=tk.LEFT, padx=(5, 5))
        self.stop_button = ttk.Button(left_controls, text="STOP", command=self._handle_stop)
        self.stop_button.pack(side=tk.LEFT, padx=5)
        self.reset_button = ttk.Button(left_controls, text="RESET", command=self._handle_reset)
        self.reset_button.pack(side=tk.LEFT, padx=5)

        right_controls = ttk.Frame(top_control_frame); right_controls.pack(side=tk.RIGHT, fill=tk.Y, padx=(10, 0))
        
        self.zoom_in_button = ttk.Button(right_controls, text="Zoom In", style='Green.TButton', command=self._zoom_in)
        self.zoom_in_button.pack(side=tk.LEFT, padx=5)
        self.zoom_out_button = ttk.Button(right_controls, text="Zoom Out", style='Green.TButton', command=self._zoom_out)
        self.zoom_out_button.pack(side=tk.LEFT, padx=5)
        
        ttk.Label(right_controls, text="Export as:").pack(side=tk.LEFT, padx=(15,5))
        self.export_combo = ttk.Combobox(right_controls, values=["Image (PNG)", "Data (CSV)"], width=15)
        self.export_combo.pack(side=tk.LEFT, padx=5)
        self.export_combo.current(0)
        self.export_button = ttk.Button(right_controls, text="Export", command=self._handle_export)
        self.export_button.pack(side=tk.LEFT, padx=5)
        
        self.status_label = ttk.Label(right_controls, text="Status: Idle", font=('Arial', 10, 'bold'), width=25, anchor=tk.W)
        self.status_label.pack(side=tk.LEFT, padx=(10, 0))
        
        middle_controls = ttk.Frame(top_control_frame); middle_controls.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        ttk.Label(middle_controls, text="Speed (ms):").pack(side=tk.LEFT, padx=(10, 5))
        self.current_speed_label = ttk.Label(middle_controls, text="", width=6)
        self.speed_scale = ttk.Scale(middle_controls, from_=1, to=500, orient=tk.HORIZONTAL, command=self._update_speed_label)
        self.speed_scale.set(50)
        self.speed_scale.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        self.current_speed_label.pack(side=tk.LEFT, padx=5)
        self._update_speed_label(50)
        
        plot_container_frame = ttk.Frame(self); plot_container_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_container_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.scrollbar = ttk.Scrollbar(plot_container_frame, orient=tk.VERTICAL, command=self._scroll_plot)
        # Bind mouse wheel to the canvas for scrolling
        self.canvas_widget.bind("<MouseWheel>", self._on_mouse_wheel)  # Windows/macOS
        self.canvas_widget.bind("<Button-4>", self._on_mouse_wheel)    # Linux (up)
        self.canvas_widget.bind("<Button-5>", self._on_mouse_wheel)    # Linux (down)

    def _update_button_states(self):
        # Updates the state of buttons (enabled/disabled) based on what the app is doing.
        is_simulating = self.simulation_running
        is_user_control_active = not is_simulating
        
        if self.simulation_finished or is_simulating:
            self.run_button.config(state=tk.DISABLED, style='Red.TButton')
        elif self.simulation_paused:
            self.run_button.config(text="CONTINUE", state=tk.NORMAL, style='Green.TButton')
        else:
            self.run_button.config(text="RUN", state=tk.NORMAL, style='Green.TButton')
        
        self.stop_button.config(state=tk.NORMAL if is_simulating else tk.DISABLED, style='Green.TButton' if is_simulating else 'Red.TButton')
        self.reset_button.config(state=tk.NORMAL if not is_simulating else tk.DISABLED, style='Green.TButton' if not is_simulating else 'Red.TButton')
        
        entry_state = tk.NORMAL if is_user_control_active else tk.DISABLED
        self.n_entry.config(state=entry_state)
        self.steps_entry.config(state=entry_state)
        
        # The export button is enabled ONLY after a simulation is finished
        export_state = tk.NORMAL if self.simulation_finished else tk.DISABLED
        self.export_button.config(state=export_state, style='Green.TButton' if self.simulation_finished else 'Red.TButton')
        self.export_combo.config(state='readonly' if self.simulation_finished else 'disabled')

        enable_zoom = bool(self.automaton and self.automaton.history) and is_user_control_active
        self.zoom_in_button.config(state=tk.NORMAL if enable_zoom and self.max_display_rows > self.min_display_rows else tk.DISABLED)
        total_rows = len(self.automaton.history) if self.automaton else 0
        self.zoom_out_button.config(state=tk.NORMAL if enable_zoom and self.max_display_rows < total_rows else tk.DISABLED)

    def _handle_run(self):
        # What happens when we press the RUN or CONTINUE button.
        if self.simulation_paused:
            # If paused, we just continue.
            self.simulation_running = True
            self.simulation_paused = False
            self.status_label.config(text="Status: Running")
            self._update_button_states()
            self._simulation_step()
            return

        # If we start a new simulation, we validate the inputs.
        error_messages = []
        n, steps = None, None

        # Unified check for N
        try:
            n_val = int(self.n_entry.get())
            if n_val < 2:
                error_messages.append("• The number of soldiers (N) must be 2 or greater.")
            else:
                n = n_val
        except ValueError:
            error_messages.append("• The value for 'N' must be a valid integer.")
        
        # Unified check for Steps
        try:
            steps_val = int(self.steps_entry.get())
            if steps_val < 1:
                error_messages.append("• The steps per update must be 1 or greater.")
            else:
                steps = steps_val
        except ValueError:
            error_messages.append("• The value for 'Steps' must be a valid integer.")
        
        # If there are errors, show ONE window with all the errors.
        if error_messages:
            final_message = "Please correct the following errors:\n\n" + "\n".join(error_messages)
            messagebox.showerror("Invalid Input", final_message)
            return

        # If everything is correct, we begin.
        selected_algo_name = self.algo_combo.get()
        algo_class = self.algorithm_map[selected_algo_name]
        self.automaton = FiringSquadAutomaton(n, algo_class())
        self.steps_per_update = steps
        
        self.simulation_running = True
        self.simulation_paused = False
        self.simulation_finished = False
        self.current_scroll_pos = 0
        self.status_label.config(text="Status: Running")
        
        self._update_button_states()
        self.draw_spacetime_diagram()
        self._simulation_step()

    def _handle_stop(self):
        # What happens when we press STOP.
        if self.simulation_running:
            self.simulation_running = False
            self.simulation_paused = True
            self.status_label.config(text="Status: Paused")
            self._update_button_states()
            self._update_scrollbar_visibility()

    def _handle_reset(self):
        # What happens when we press RESET.
        self.simulation_running = False
        self._reset_application_state()
        self.draw_spacetime_diagram()
        self.status_label.config(text="Status: Idle")
    
    def _reset_application_state(self):
        # Resets all variables to their initial values.
        self.automaton = None
        self.simulation_running = False
        self.simulation_paused = False
        self.simulation_finished = False
        self.current_scroll_pos = 0
        self.max_display_rows = 25  
        self.scrollbar.pack_forget()
        self._update_button_states()
        
    def _simulation_step(self):
        # The main "loop" of the simulation.
        if not self.simulation_running: return # If STOP was pressed, stop
        
        # Execute a batch of steps for speed
        for _ in range(self.steps_per_update):
            if self.automaton.is_firing: break
            self.automaton.step()
        
        # Check for Timeout, as a safety net
        max_time = self.automaton.algorithm.get_time_complexity(self.automaton.n)
        if self.automaton.time >= max_time and not self.automaton.is_firing:
            self.simulation_finished = True
            self.simulation_running = False
            messagebox.showwarning("Timeout", f"Simulation stopped due to exceeding max time (T ≥ {max_time}).")
            self.status_label.config(text=f"Status: Timeout at T={self.automaton.time}")
            self._update_button_states()
            self._update_scrollbar_visibility()
            self.draw_spacetime_diagram()
            return
        
        self._update_gui()
        
        # If everyone has fired...
        if self.automaton.is_firing:
            self.simulation_finished = True
            self.simulation_running = False
            self.status_label.config(text=f"FIRING COMPLETE! (T={self.automaton.time})")
            self._update_button_states()
            self._update_scrollbar_visibility()
            self.draw_spacetime_diagram()
        else:
            # If not, call ourselves again after a delay
            self.after(int(self.speed_scale.get()), self._simulation_step)

    def _update_gui(self):
        # Updates the graphics.
        self._update_scrollbar_visibility()
        total_rows = len(self.automaton.history)
        # Auto-scroll to the end
        if not self.simulation_paused:
            self.current_scroll_pos = max(0, total_rows - self.max_display_rows)
        self._update_scrollbar_view()
        self.draw_spacetime_diagram()
        if self.simulation_running:
            self.status_label.config(text=f"Status: Running (T={self.automaton.time})")
    
    def _update_scrollbar_visibility(self):
        # Logic for when the scrollbar should be visible
        is_needed = self.automaton and len(self.automaton.history) > self.max_display_rows
        if not self.simulation_running and is_needed:
            self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        else:
            self.scrollbar.pack_forget()

    def draw_spacetime_diagram(self):
        # The method that does all the work for drawing the diagram
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        if self.automaton is None or not self.automaton.history:
            self.ax.set_facecolor('white')
            self.ax.set_title("Firing Squad Spacetime Diagram")
            self.ax.set_xlabel("Cell (Soldier)")
            self.ax.set_ylabel("Time (Step)")
            self.canvas.draw()
            return
            
        colors = self.automaton.algorithm.get_colors()
        display_names = self.automaton.algorithm.get_display_names()
        history_matrix = np.array(self.automaton.history)
        total_rows, n_cols = history_matrix.shape
        start_row = self.current_scroll_pos
        end_row = min(start_row + self.max_display_rows, total_rows)
        display_matrix = history_matrix[start_row:end_row, :]
        num_displayed_rows = display_matrix.shape[0]
        
        if num_displayed_rows <= 0: self.canvas.draw(); return
        
        colors_list = [colors[i] for i in sorted(colors.keys())]
        cmap = ListedColormap(colors_list)
        bounds = np.arange(len(colors) + 1) - 0.5
        norm = BoundaryNorm(bounds, cmap.N)
        
        im = self.ax.imshow(display_matrix, cmap=cmap, norm=norm, interpolation='nearest', aspect='equal', extent=[-0.5, n_cols - 0.5, num_displayed_rows - 0.5, -0.5])
        
        self.ax.set_title(f"FSSP Spacetime Diagram (N={n_cols})")
        self.ax.set_xlabel("Cell (Soldier)")
        self.ax.set_ylabel("Time (Step)")
        self.ax.set_xlim(-0.5, n_cols - 0.5)
        self.ax.set_ylim(num_displayed_rows - 0.5, -0.5)
        
        y_ticks = np.linspace(0, num_displayed_rows-1, min(num_displayed_rows, 10), dtype=int)
        self.ax.set_yticks(y_ticks)
        self.ax.set_yticklabels([yt + start_row for yt in y_ticks])
        self.ax.set_xticks(np.arange(0, n_cols, max(1, n_cols // 10)))
        self.ax.set_xticks(np.arange(-.5, n_cols, 1), minor=True)
        self.ax.set_yticks(np.arange(-.5, num_displayed_rows, 1), minor=True)
        self.ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5)
        self.ax.tick_params(which='minor', size=0)
        
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = self.fig.colorbar(im, cax=cax, ticks=np.arange(len(colors)))
        cbar.ax.set_yticklabels([display_names[i] for i in sorted(display_names.keys())])
        cbar.set_label("State")
        
        self.fig.tight_layout(rect=[0, 0.03, 1, 0.97])
        self.canvas.draw()
        self._update_button_states()

    def _update_speed_label(self, value):
        self.current_speed_label.config(text=f"{int(float(value))} ms")
        
    def _update_scrollbar_view(self):
        if not self.automaton or not self.automaton.history: return
        total_rows = len(self.automaton.history)
        if total_rows <= self.max_display_rows: self.scrollbar.set(0, 1); return
        view_start = self.current_scroll_pos / total_rows
        thumb_size = self.max_display_rows / total_rows
        self.scrollbar.set(view_start, view_start + thumb_size)
        
    def _on_mouse_wheel(self, event):
        # Logic for scrolling with the mouse wheel
        if self.simulation_running or not self.automaton or not self.automaton.history: return
        
        if event.num == 4 or event.delta > 0: scroll_amount = -3
        elif event.num == 5 or event.delta < 0: scroll_amount = 3
        else: return
        
        self.current_scroll_pos += scroll_amount
        total_rows = len(self.automaton.history)
        max_scroll = max(0, total_rows - self.max_display_rows)
        self.current_scroll_pos = max(0, min(self.current_scroll_pos, max_scroll))
        
        self._update_scrollbar_view()
        self.draw_spacetime_diagram()
        
    def _scroll_plot(self, *args):
        # Logic for scrolling with the scrollbar
        if not self.automaton or not self.automaton.history or self.simulation_running: return
        total_rows = len(self.automaton.history)
        max_scroll = max(0, total_rows - self.max_display_rows)
        if args[0] == 'moveto': self.current_scroll_pos = int(float(args[1]) * max_scroll)
        elif args[0] == 'scroll': self.current_scroll_pos += int(args[1])
        self.current_scroll_pos = max(0, min(self.current_scroll_pos, max_scroll))
        self._update_scrollbar_view()
        self.draw_spacetime_diagram()

    def _zoom_in(self):
        if self.max_display_rows > self.min_display_rows:
            self.max_display_rows = max(self.min_display_rows, self.max_display_rows - self.zoom_step)
            self._handle_zoom_update()

    def _zoom_out(self):
        if not self.automaton or not self.automaton.history: return
        total_rows = len(self.automaton.history)
        if self.max_display_rows < total_rows:
            self.max_display_rows += self.zoom_step
            self._handle_zoom_update()
            
    def _handle_zoom_update(self):
        # Updates the view after a zoom in/out
        if self.automaton and self.automaton.history:
            total_rows = len(self.automaton.history)
            self.current_scroll_pos = min(self.current_scroll_pos, max(0, total_rows - self.max_display_rows))
        self._update_scrollbar_visibility()
        self._update_scrollbar_view()
        self.draw_spacetime_diagram()

    def _on_closing(self):
        # Interrupts the simulation before closing the window
        self.simulation_running = False
        self.destroy()

    def _handle_export(self):
        # Logic for the "Export" button
        if not self.automaton or not self.simulation_finished:
            messagebox.showwarning("Export", "There is no completed simulation to export.")
            return
        
        export_format = self.export_combo.get()
        
        if export_format == "Image (PNG)":
            file_path = asksaveasfilename(title="Save Image", defaultextension=".png", filetypes=[("PNG Image", "*.png"), ("All Files", "*.*")])
            if file_path:
                try:
                    self.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                    messagebox.showinfo("Success", f"Image successfully saved to:\n{file_path}")
                except Exception as e:
                    messagebox.showerror("Save Error", f"Could not save the image:\n{e}")
                    
        elif export_format == "Data (CSV)":
            file_path = asksaveasfilename(title="Save Data", defaultextension=".csv", filetypes=[("CSV File", "*.csv"), ("All Files", "*.*")])
            if file_path:
                try:
                    state_map = {v: k for k, v in self.automaton.states.items()}
                    header = [f"Cell {i}" for i in range(self.automaton.n)]
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:
                        writer = csv.writer(f)
                        writer.writerow(["Time"] + header)
                        for time_step, row in enumerate(self.automaton.history):
                            string_row = [state_map.get(state, 'UNK') for state in row]
                            writer.writerow([time_step] + string_row)
                    messagebox.showinfo("Success", f"Data successfully saved to:\n{file_path}")
                except Exception as e:
                    messagebox.showerror("Save Error", f"Could not save the data:\n{e}")

if __name__ == "__main__":
    app = App()
    app.mainloop()
