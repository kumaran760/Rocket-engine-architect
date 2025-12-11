import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Rectangle
from scipy.optimize import fsolve
import csv

# --- Configuration & Constants ---
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")
G0 = 9.80665  # m/s^2

# --- 1. Physics Engine ---
class RocketEngine:
    """
    Handles thermodynamic, geometric, and structural calculations.
    """
    def __init__(self):
        self.results = {}
        self.geometry = {}

    def prandtl_meyer(self, M, gamma):
        if M <= 1.0: return 0.0
        term1 = np.sqrt((gamma + 1) / (gamma - 1))
        term2 = np.arctan(np.sqrt((gamma - 1) / (gamma + 1) * (M**2 - 1)))
        term3 = np.arctan(np.sqrt(M**2 - 1))
        return term1 * term2 - term3

    def solve_mach_from_pressure_ratio(self, pr_ratio, gamma):
        def equation(M):
            return (1 + (gamma - 1)/2 * M**2)**(gamma / (gamma - 1)) - pr_ratio
        M_sol = fsolve(equation, 2.0)
        return M_sol[0]

    def calculate_design(self, inputs):
        # --- A. Unpack Inputs ---
        F = inputs['thrust'] * 1000      # kN -> N
        Pc = inputs['pc'] * 1e5          # Bar -> Pa
        Pa = inputs['pa'] * 1e5          # Bar -> Pa
        gamma = inputs['gamma']
        R = inputs['R']
        Tc = inputs['Tc']
        SF = inputs['sf']
        sigma_y = inputs['sigma_y'] * 1e6 # MPa -> Pa
        
        # Chamber specific inputs
        L_star = inputs['L_star']         # m
        contraction_ratio = inputs['contraction_ratio'] # Ac / At

        # --- B. Thermodynamics ---
        # 1. Optimal Expansion Condition (Pe = Pa)
        Pe = Pa 
        pr_ratio = Pc / Pe
        Me = self.solve_mach_from_pressure_ratio(pr_ratio, gamma)

        # 2. Performance Parameters
        # Specific Impulse (Theoretical)
        # Isp = c* * Cf / g0
        # C_star (Characteristic Velocity)
        gamma_term = np.sqrt(gamma) * (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
        c_star = np.sqrt(R * Tc) / gamma_term
        
        # Cf (Thrust Coefficient)
        term_a = (2 * gamma**2) / (gamma - 1)
        term_b = (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))
        term_c = 1 - (Pe / Pc) ** ((gamma - 1) / gamma)
        Cf = np.sqrt(term_a * term_b * term_c)
        
        Isp = (c_star * Cf) / G0

        # 3. Sizing (Throat & Exit)
        At = F / (Pc * Cf)
        Dt = 2 * np.sqrt(At / np.pi)
        
        # Exit Area calculation using Mach-Area relation
        def area_mach(M, g):
            return (1/M) * ((2/(g+1)) * (1 + (g-1)/2 * M**2))**((g+1)/(2*(g-1)))
        epsilon = area_mach(Me, gamma) # Ae / At
        Ae = At * epsilon
        De = 2 * np.sqrt(Ae / np.pi)

        # Mass Flow
        mdot = F / (Isp * G0)

        # --- C. Combustion Chamber Sizing ---
        # Vc = L* * At
        Vc = L_star * At
        
        # Chamber Area (Ac) determined by Contraction Ratio
        Ac = contraction_ratio * At
        Dc = 2 * np.sqrt(Ac / np.pi)
        
        # Length of the cylindrical section (L_cyl)
        # Simplified: Vc = Ac * Lc_cyl + V_converging. 
        # Ignoring converging volume for safety margin (slightly larger chamber is safer)
        # Lc_cyl = Vc / Ac
        Lc_cyl = Vc / Ac 

        # --- D. Structural (Hoop Stress) ---
        # t = (P * r) / (sigma_y / SF)
        t_chamber = (Pc * (Dc/2)) / (sigma_y / SF)
        t_nozzle = (Pe * (De/2)) / (sigma_y / SF)
        t_wall_min = max(t_chamber, t_nozzle)

        # --- E. Geometry Generation ---
        nozzle_type = inputs['type']
        
        # 1. Chamber Geometry (Common for Bell types)
        # Cylinder from x = -(Lc_cyl + Converging_Len) to ...
        # Converging length approx: assume 45 deg converging angle
        L_conv = (Dc/2 - Dt/2) / np.tan(np.radians(45))
        
        x_cham = np.linspace(-(Lc_cyl + L_conv), -L_conv, 20)
        y_cham = [Dc/2] * len(x_cham)
        
        x_conv = np.linspace(-L_conv, 0, 20)
        # Linear converging section
        y_conv = np.linspace(Dc/2, Dt/2, 20)

        x_geom = []
        y_geom = []
        
        if nozzle_type == "Aerospike":
            # Aerospike Logic (Plug Nozzle)
            # Length approx based on expansion
            L_spike = 0.8 * (Dt/2) * epsilon 
            
            mach_steps = np.linspace(1.0, Me, 100)
            nu_Me = self.prandtl_meyer(Me, gamma)
            
            x_spike = []
            y_spike = []
            
            for M_loc in mach_steps:
                nu = self.prandtl_meyer(M_loc, gamma)
                # Angle relative to vertical axis for plug
                theta = nu_Me - nu 
                
                # Geometric approximation for Spike contour
                # x grows, y shrinks from Rt to 0 (ideally)
                val = (M_loc - 1) / (Me - 1)
                x_curr = val * L_spike
                
                # Simple contour shaping based on expansion area
                # A_local / At ~ (y^2) for plug? 
                # Visual approximation:
                y_curr = (Dt/2) * (1 - val)
                
                x_spike.append(x_curr)
                y_spike.append(y_curr)
                
            self.geometry = {
                "type": "Aerospike",
                "x_spike": x_spike,
                "y_spike": y_spike,
                "Rt": Dt/2,
                "Re": De/2, # Represents Cowl Radius here
                "Dc": Dc # For relative scale
            }
            
        else:
            # Bell / Dual Bell
            Rt = Dt / 2
            Re = De / 2
            
            if nozzle_type == "Dual-Bell":
                # 15 degree cone length reference
                L_cone = (Re - Rt) / np.tan(np.radians(15))
                Ln = 0.85 * L_cone # slightly longer for dual bell
                
                # Inflection point at 60%
                L_1 = 0.6 * Ln
                y_1 = Rt + (Re - Rt)*0.5 # Arbitrary inflection height
                
                x_n1 = np.linspace(0, L_1, 40)
                y_n1 = Rt + (y_1 - Rt) * (x_n1/L_1)**0.8
                
                x_n2 = np.linspace(L_1, Ln, 40)
                y_n2 = y_1 + (Re - y_1) * ((x_n2 - L_1)/(Ln - L_1))**0.7
                
                x_nozz = np.concatenate([x_n1, x_n2])
                y_nozz = np.concatenate([y_n1, y_n2])
                
            else: # Standard Bell
                # 15 deg cone length
                L_cone = (Re - Rt) / np.tan(np.radians(15))
                Ln = 0.8 * L_cone
                
                x_nozz = np.linspace(0, Ln, 50)
                # Rao parabola: y = Rt + (Re-Rt)*sqrt(x/Ln) is a rough visual approx
                # More accurate: Quadratic fit matching angles
                y_nozz = Rt + (Re - Rt) * np.sqrt(x_nozz / Ln)

            self.geometry = {
                "type": "Bell",
                "x": np.concatenate([x_cham, x_conv, x_nozz]),
                "y": np.concatenate([y_cham, y_conv, y_nozz]),
                "t_wall": t_wall_min,
                "Dc": Dc,
                "Lc": Lc_cyl
            }

        # --- F. Results Packaging ---
        self.results = {
            "Thrust (kN)": F / 1000,
            "Isp (s)": Isp,
            "Mass Flow (kg/s)": mdot,
            "Throat Dia (mm)": Dt * 1000,
            "Exit Dia (mm)": De * 1000,
            "Expansion Ratio": epsilon,
            "Chamber Dia (mm)": Dc * 1000,
            "Chamber Len (mm)": Lc_cyl * 1000,
            "Wall Thick (mm)": t_wall_min * 1000
        }

        return self.results

# --- 2. GUI Application ---
class RocketDesignApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        
        self.title("AeroSpace Engine Architect v2.0")
        self.geometry("1280x850")
        
        self.engine = RocketEngine()
        
        # Grid Layout
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.create_sidebar()
        self.create_main_view()

    def create_sidebar(self):
        # Make Sidebar Scrollable to fit all new inputs
        self.sidebar = ctk.CTkScrollableFrame(self, width=320, corner_radius=0, label_text="Design Parameters")
        self.sidebar.grid(row=0, column=0, sticky="nsew")
        
        # --- A. Mission ---
        self.add_section_header("1. Mission Requirements")
        self.entry_thrust = self.add_labeled_entry("Target Thrust (kN):", "100")
        self.entry_pc = self.add_labeled_entry("Chamber Pressure (Bar):", "50")
        self.entry_pa = self.add_labeled_entry("Ambient Pressure (Bar):", "1.013")
        self.entry_sf = self.add_labeled_entry("Safety Factor:", "1.5")

        # --- B. Propellant ---
        self.add_section_header("2. Propellant Chemistry")
        self.combo_prop = ctk.CTkComboBox(self.sidebar, values=["LOX/RP-1", "LOX/LH2", "N2O4/MMH", "Custom"],
                                          command=self.update_prop_props)
        self.combo_prop.pack(padx=20, pady=5, fill="x")
        
        self.frame_prop = ctk.CTkFrame(self.sidebar, fg_color="transparent")
        self.frame_prop.pack(padx=10, pady=5, fill="x")
        self.entry_gamma = self.add_grid_entry(self.frame_prop, "Gamma:", "1.22", 0, 0)
        self.entry_tc = self.add_grid_entry(self.frame_prop, "Tc (K):", "3400", 0, 1)
        self.entry_r = self.add_grid_entry(self.frame_prop, "R (J/kgK):", "300", 1, 0)

        # --- C. Chamber Geometry ---
        self.add_section_header("3. Chamber Sizing")
        self.entry_lstar = self.add_labeled_entry("Char. Length L* (m):", "1.1")
        self.entry_cr = self.add_labeled_entry("Contraction Ratio (Ac/At):", "3.0")

        # --- D. Material ---
        self.add_section_header("4. Material Properties")
        self.combo_mat = ctk.CTkComboBox(self.sidebar, values=["Inconel 718", "SS 304", "Ti-6Al-4V", "Copper C18200"],
                                         command=self.update_mat_props)
        self.combo_mat.pack(padx=20, pady=5, fill="x")
        
        self.frame_mat = ctk.CTkFrame(self.sidebar, fg_color="transparent")
        self.frame_mat.pack(padx=10, pady=5, fill="x")
        self.entry_yield = self.add_grid_entry(self.frame_mat, "Yield (MPa):", "1035", 0, 0)
        self.entry_rho = self.add_grid_entry(self.frame_mat, "Density (kg/m3):", "8190", 0, 1)

        # --- E. Nozzle Type ---
        self.add_section_header("5. Nozzle Topology")
        self.radio_var = tk.StringVar(value="Bell")
        ctk.CTkRadioButton(self.sidebar, text="Bell (Rao)", variable=self.radio_var, value="Bell").pack(padx=20, pady=2, anchor="w")
        ctk.CTkRadioButton(self.sidebar, text="Dual-Bell", variable=self.radio_var, value="Dual-Bell").pack(padx=20, pady=2, anchor="w")
        ctk.CTkRadioButton(self.sidebar, text="Aerospike", variable=self.radio_var, value="Aerospike").pack(padx=20, pady=2, anchor="w")

        # --- Buttons ---
        self.btn_calc = ctk.CTkButton(self.sidebar, text="CALCULATE & PLOT", height=40, fg_color="#1f538d", 
                                      font=("Arial", 14, "bold"), command=self.run_calculation)
        self.btn_calc.pack(padx=20, pady=(30, 10), fill="x")

        self.btn_exp = ctk.CTkButton(self.sidebar, text="Export Geometry .CSV", fg_color="transparent", 
                                     border_width=1, command=self.export_csv)
        self.btn_exp.pack(padx=20, pady=10, fill="x")

    def create_main_view(self):
        self.main_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.main_frame.grid(row=0, column=1, sticky="nsew", padx=20, pady=20)
        
        # Results Grid
        self.results_frame = ctk.CTkFrame(self.main_frame)
        self.results_frame.pack(fill="x", pady=(0, 20))
        
        # We will dynamically update these
        self.res_labels = {}
        keys = ["Throat Dia (mm)", "Exit Dia (mm)", "Expansion Ratio", 
                "Chamber Dia (mm)", "Chamber Len (mm)", "Wall Thick (mm)",
                "Mass Flow (kg/s)", "Isp (s)"]
        
        for i, key in enumerate(keys):
            f = ctk.CTkFrame(self.results_frame, fg_color="transparent")
            row = 0 if i < 4 else 1
            col = i % 4
            f.grid(row=row, column=col, sticky="nsew", padx=10, pady=10)
            self.results_frame.grid_columnconfigure(col, weight=1)
            
            ctk.CTkLabel(f, text=key, font=("Arial", 11, "bold"), text_color="gray").pack()
            l = ctk.CTkLabel(f, text="--", font=("Arial", 16, "bold"), text_color="#4facfe")
            l.pack()
            self.res_labels[key] = l

        # Plot
        self.fig, self.ax = plt.subplots(figsize=(6, 5), dpi=100)
        self.fig.patch.set_facecolor('#2b2b2b')
        self.ax.set_facecolor('#2b2b2b')
        self.ax.tick_params(colors='white', which='both')
        self.ax.xaxis.label.set_color('white')
        self.ax.yaxis.label.set_color('white')
        for spine in self.ax.spines.values(): spine.set_color('white')

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.main_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    # --- GUI Helpers ---
    def add_section_header(self, text):
        ctk.CTkLabel(self.sidebar, text=text, font=("Arial", 14, "bold"), anchor="w").pack(padx=10, pady=(20, 5), fill="x")

    def add_labeled_entry(self, text, default):
        ctk.CTkLabel(self.sidebar, text=text, anchor="w", font=("Arial", 12)).pack(padx=20, pady=(2, 0), fill="x")
        e = ctk.CTkEntry(self.sidebar)
        e.insert(0, default)
        e.pack(padx=20, pady=(0, 5), fill="x")
        return e

    def add_grid_entry(self, parent, text, default, r, c):
        f = ctk.CTkFrame(parent, fg_color="transparent")
        f.grid(row=r, column=c, padx=5, pady=5, sticky="ew")
        ctk.CTkLabel(f, text=text, font=("Arial", 11)).pack(anchor="w")
        e = ctk.CTkEntry(f, width=80)
        e.insert(0, default)
        e.pack(fill="x")
        return e

    # --- Logic ---
    def update_prop_props(self, choice):
        props = {
            "LOX/RP-1": ("1.22", "3400", "300"),
            "LOX/LH2": ("1.26", "3600", "520"),
            "N2O4/MMH": ("1.20", "3100", "290"),
        }
        if choice in props:
            vals = props[choice]
            self.entry_gamma.delete(0, "end"); self.entry_gamma.insert(0, vals[0])
            self.entry_tc.delete(0, "end"); self.entry_tc.insert(0, vals[1])
            self.entry_r.delete(0, "end"); self.entry_r.insert(0, vals[2])

    def update_mat_props(self, choice):
        mat_data = {
            "Inconel 718": ("1035", "8190"),
            "SS 304": ("205", "8000"),
            "Ti-6Al-4V": ("880", "4430"),
            "Copper C18200": ("350", "8940")
        }
        if choice in mat_data:
            vals = mat_data[choice]
            self.entry_yield.delete(0, "end"); self.entry_yield.insert(0, vals[0])
            self.entry_rho.delete(0, "end"); self.entry_rho.insert(0, vals[1])

    def run_calculation(self):
        try:
            inputs = {
                "thrust": float(self.entry_thrust.get()),
                "pc": float(self.entry_pc.get()),
                "pa": float(self.entry_pa.get()),
                "sf": float(self.entry_sf.get()),
                "gamma": float(self.entry_gamma.get()),
                "Tc": float(self.entry_tc.get()),
                "R": float(self.entry_r.get()),
                "L_star": float(self.entry_lstar.get()),
                "contraction_ratio": float(self.entry_cr.get()),
                "sigma_y": float(self.entry_yield.get()),
                "rho_mat": float(self.entry_rho.get()),
                "type": self.radio_var.get()
            }

            results = self.engine.calculate_design(inputs)
            geo = self.engine.geometry

            # Update Labels
            for key, val in results.items():
                if key in self.res_labels:
                    if "Dia" in key or "Len" in key or "Thick" in key:
                         self.res_labels[key].configure(text=f"{val:.1f}")
                    elif "Flow" in key:
                         self.res_labels[key].configure(text=f"{val:.2f}")
                    else:
                         self.res_labels[key].configure(text=f"{val:.1f}")

            # Plotting
            self.ax.clear()
            self.ax.set_title(f"Engine Profile: {inputs['type']}", color='white', pad=20)
            self.ax.set_xlabel("Axial Length (m)")
            self.ax.set_ylabel("Radial Height (m)")

            if geo['type'] == "Aerospike":
                # Plot Spike
                self.ax.plot(geo['x_spike'], geo['y_spike'], color='#ff0066', linewidth=2, label="Spike")
                self.ax.plot(geo['x_spike'], [-y for y in geo['y_spike']], color='#ff0066', linewidth=2)
                # Plot Cowl Ref
                self.ax.plot([0, 0], [geo['Rt'], geo['Re']], color='#00ffcc', linewidth=4, label="Cowl")
                self.ax.plot([0, 0], [-geo['Rt'], -geo['Re']], color='#00ffcc', linewidth=4)
                
            else:
                # Bell / Dual Bell
                self.ax.plot(geo['x'], geo['y'], color='#00ffcc', linewidth=2, label="Inner Wall")
                self.ax.plot(geo['x'], [-y for y in geo['y']], color='#00ffcc', linewidth=2)
                
                # Outer Wall (Thickness)
                t = geo['t_wall']
                y_outer = [y + t for y in geo['y']]
                self.ax.fill_between(geo['x'], geo['y'], y_outer, color='#444444', alpha=0.5, label=f"Wall ({t*1000:.1f}mm)")
                
                # Highlight Chamber
                if 'Lc' in geo:
                    rect = Rectangle((-geo['Lc'], -geo['Dc']/2), geo['Lc'], geo['Dc'], 
                                     linewidth=1, edgecolor='yellow', facecolor='none', linestyle='--', label="Chamber Vol")
                    self.ax.add_patch(rect)

            self.ax.axis('equal')
            self.ax.legend(loc='upper right', facecolor='#2b2b2b', labelcolor='white')
            self.ax.grid(True, linestyle=':', alpha=0.3)
            self.canvas.draw()

        except ValueError:
            messagebox.showerror("Input Error", "Please ensure all fields contain valid numbers.")

    def export_csv(self):
        geo = self.engine.geometry
        if not geo: return
        filename = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if filename:
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                if geo['type'] == "Aerospike":
                    writer.writerow(["X", "Y_Spike"])
                    for x, y in zip(geo['x_spike'], geo['y_spike']): writer.writerow([x, y])
                else:
                    writer.writerow(["X", "Y_Inner", "Thickness"])
                    for x, y in zip(geo['x'], geo['y']): writer.writerow([x, y, geo['t_wall']])

if __name__ == "__main__":
    app = RocketDesignApp()
    app.mainloop()