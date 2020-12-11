# import numpy as np
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Simulation Definition:
# Type of Simulation:
EIS_sim = 0  # Flag for EIS simulation
Polarization_sim = 1  # Flag for EIS simulation
# EIS Simulation Parameters
freq = 1  # [rad/s]
# freq = np.array[0.001, 0.1, 1, 10]  # [rad/s]
CurrentAmplitude = 0.01  # [Amps]
# Polarization Simulation Parameters
C_rate = 0.1  # How many charges per hour????%%%%%%%%%%
# C_rate = np.array[0.1, 1, 2, 5, 8]
charge_frac = 0.9  # How deep do we want to charge/discharge?
# Overall Simulation Parameters
SOC_start = 100  # [%], Initial state of charge of the cell
Temp_start = 298  # K, Initial temperature of the cell
Temp_Dep = 0  # Flag for temperature dependence
VoltageMax = 3.5  # [V]
VoltageMin = 2.25  # [V]
AnodeFormation_X = 0.2  # [-]
CathodeFormation_X = 0.9  # [-]

# Cell Geometry
# Anode
del_x_an = 0.1  # [m], Direction without tab, length of electrode
del_y_an = 0.1  # [m], Direction with tab, height of electrode
del_z_an = 0.1  # [m], Direction normal to current flow, thickness of electrode
r_p_an = 1e-6  # [m], outer radius of electrode particle
eps_an = 0.7  # [-], volume fraction of active material
tau_an = 1  # [-?], tortuosity of active material
density_an = 2260  # [kg/m3], density of the active material
# Cathode
del_x_ca = 0.1  # [m], Direction without tab, length of electrode
del_y_ca = 0.1  # [m], Direction with tab, height of electrode
del_z_ca = 0.1  # [m], Direction normal to current flow, thickness of electrode
r_p_ca = 1e-6  # [m], outer radius of electrode particle
eps_ca = 0.7  # [-], volume fraction of active material in the electrode cell
tau_ca = 1  # [-?], tortuosity of active material in the electrode cell
density_ca = 2292  # [kg/m3], density of the active material

# Cell Performance
# Anode
C_dl_an = 1e4  # [F/m2]
R_sei_an = 0.1  # [units], solid electrolyte interface resistance
i_o_an = 4.0  # [A/m2], exchange current density
n_an = -1  # [equivalents], charge transferred to the electrode
beta_an = 0.5  # [-], Symmetry factor?????
capacity_graphite = 350  # [Ah/kg]
D_o_an = 1e-16  # [?], Solid-state diffusion coefficient
# Separator
# conductivity?
# sigma_sep = 0.1  # [?], ionic conductivity
# Cathode
C_dl_ca = 1e4  # [F/m2]
R_sei_ca = 0.1  # [units], solid electrolyte interface resistance
i_o_ca = 100  # [A/m2], exchange current density
n_ca = -1  # [equivalents], charge transferred to the electrode
beta_ca = 0.5  # [-], Symmetry factor?????
capacity_LFP = 175  # [Ah/kg]
D_o_ca = 1e-16  # [?], Solid-state diffusion coefficient

# Numerical Parameters
N_part = 1  # Number of nodes in solid particle
N_an = 1  # Number of nodes in the anode
N_sep = 1  # Number of nodes in the separator, Currently is set for only 1 sep node
N_ca = 1  # Number of nodes in the cathode
