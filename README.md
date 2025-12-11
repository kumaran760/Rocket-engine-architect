# 🚀 Rocket Engine Architect

**Rocket Engine Architect** is a specialized desktop application designed for aerospace engineers and enthusiasts. It provides a streamlined interface to calculate, visualize, and export geometric specifications for liquid rocket engines based on high-level mission requirements.

![Status](https://img.shields.io/badge/Status-Active-success)
![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## 🌟 Key Features

### 1. Parametric Design Interface
- **Mission Inputs:** Define Target Thrust, Chamber Pressure ($P_c$), and Ambient Pressure ($P_a$).
- **Propellant Chemistry:** Built-in presets for LOX/RP-1, LOX/LH2, and N2O4/MMH with customizable Gamma ($\gamma$), Combustion Temp ($T_c$), and Gas Constant ($R$).
- **Material Selection:** Auto-loads material properties (Yield Strength, Density) for Inconel 718, Ti-6Al-4V, SS-304, and Copper.

### 2. Advanced Nozzle Topologies
The physics engine supports three distinct nozzle types:
- **Bell Nozzle:** Uses Rao's parabolic approximation for optimal contouring.
- **Dual-Bell Nozzle:** Simulates an altitude-compensating design with an inflection point for flow separation control.
- **Aerospike (Plug) Nozzle:** Utilizes **Prandtl-Meyer expansion** logic to generate the isentropic spike contour.

### 3. Physics & Structure
- **Thermodynamics:** Solves area-Mach relations ($A/A^*$) and Isentropic flow equations to determine Throat/Exit dimensions and Specific Impulse ($I_{sp}$).
- **Combustion Chamber:** Calculates chamber volume ($V_c$) and length ($L_c$) based on Characteristic Length ($L^*$) and Contraction Ratio.
- **Structural Integrity:** Computes required wall thickness using **Hoop Stress** formulas based on local pressure and material yield strength (with Safety Factor).

### 4. Visualization & Export
- **Real-time Plotting:** Embedded **Matplotlib** canvas renders 2D cross-sections of the engine profile.
- **CAD Ready:** Export the generated geometry point cloud to `.CSV` for import into SolidWorks, Fusion 360, or ANSYS.

## 🛠️ Tech Stack

* **Language:** Python 3.10+
* **GUI Framework:** [CustomTkinter](https://github.com/TomSchimansky/CustomTkinter) (Modern, Dark-mode native UI)
* **Computation:** NumPy, SciPy (fsolve for Mach calculations)
* **Visualization:** Matplotlib

## 📦 Installation & Usage

1.  **Clone the Repository**
    ```bash
    git clone [https://github.com/YourUsername/Rocket-Engine-Architect.git](https://github.com/YourUsername/Rocket-Engine-Architect.git)
    cd Rocket-Engine-Architect
    ```

2.  **Install Dependencies**
    ```bash
    pip install customtkinter numpy scipy matplotlib
    ```

3.  **Run the Application**
    ```bash
    python main.py
    ```

4.  **Build Executable (Optional)**
    To create a standalone `.exe`:
    ```bash
    pyinstaller --noconsole --onefile --collect-all customtkinter --name "RocketEngineApp" main.py
    ```

## 📐 Physics Validation
The tool assumes ideal gas behavior and isentropic flow for initial sizing.
* **Bell Contours:** Approximated using Rao's method (80% bell length).
* **Aerospike:** Generated via Method of Characteristics approximation (Prandtl-Meyer function).

## 📄 License
This project is open-source and available under the MIT License.
