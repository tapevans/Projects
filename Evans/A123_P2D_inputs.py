# import numpy as np

# Simulation Definition:
# Type of Simulation:
EIS_sim = 0  # Flag for EIS simulation
Polarization_sim = 1  # Flag for EIS simulation
#   - EIS Simulation Parameters
freq = 1  # [rad/s]
# freq = np.array[0.001, 0.1, 1, 10]  # [rad/s]
CurrentAmplitude = 0.01  # [Amps]
#   - Polarization Simulation Parameters
C_rate = 1  # How many charges per hour
# C_rate = np.array[0.1, 1, 2, 5, 8]
charge_frac = 0.9  # How deep do we want to charge/discharge?
ChargeOrDischarge = 1  # -1 if Charge, 1 if Discharge
# Overall Simulation Parameters
SOC_start = 100  # [%], Initial state of charge of the cell
Temp_start = 298  # K, Initial temperature of the cell
Temp_Dep = 0  # Flag for temperature dependence
VoltageMax = 3.60  # [V]
VoltageMin = 2.25  # [V]
AnodeFormation_X = 0.9  # [-]
CathodeFormation_X = 0.2  # [-]

# Cell Geometry
#   - Anode
del_x_an        = 690.98e-3  # [m], Direction without tab, length of electrode
del_y_an        = 57e-3      # [m], Direction with tab, height of electrode
H_an            = 76e-6      # [m], Direction normal to current flow, thickness of electrode
r_p_an          = 5e-7       # [m], outer radius of electrode particle
eps_an          = 0.414      # [-], volume fraction of active material
tau_prefac_an   = 2.25       # [-], Bruggeman pre-exponential multiplier
bruggeman_an    = 0.5        # [-], Bruggeman exponential factor
density_an      = 2250       # [kg/m^3], density of the active material
SA_gram_an      = 5          # [m^2/s], surface area of active material per gram of active material
#   - Cathode
del_x_ca        = 653.98e-3 # [m], Direction without tab, length of electrode
del_y_ca        = 55.5e-3   # [m], Direction with tab, height of electrode
H_ca            = 158e-6    # [m], Direction normal to current flow, thickness of electrode
r_p_ca          = 5e-7      # [m], outer radius of electrode particle
eps_ca          = 0.562     # [-], volume fraction of active material in the electrode cell
tau_prefac_ca   = 2.25      # [-], Bruggeman pre-exponential multiplier
bruggeman_ca    = 0.5       # Bruggeman exponential factor
density_ca      = 3600      # [kg/m^3], density of the active material
SA_gram_ca      = 10        # [m^2/s], surface area of active material per gram of active material
#   - Separator
eps_sep         = 0.4       # [-] porosity
tau_prefac_sep  = 1         # [-], Bruggeman pre-exponential multiplier
bruggeman_sep   = 0.5       # Bruggeman exponential factor
density_sep     = 950       # [kg/m^3] density of the solid material
H_sep           = 25e-6     # [m] thickness

# Cell Performance
#   - Anode
capacity_graphite = 370  # [Ah/kg]
C_max_graphite = 17000  # [mol/m^3] Guess value
el_cond_an = 100000  # [S/m]
C_dl_an = 0.01  # [F/m^2]
R_sei_an = 0  # [ohm m^2], solid electrolyte interface resistance, 0.027 <-- BDS
i_o_an = 1.5  # [A/m^2], exchange current density
n_an = -1  # [equivalents], charge transferred to the electrode
beta_an = 0.5  # [-], Symmetry factor
D_o_an = 3e-16  # [m^2/s], Solid-state diffusion coefficient
#   - Separator
# conductivity?
# sigma_sep = 0.1  # [?], ionic conductivity
#   - Cathode
capacity_LFP = 170  # [Ah/kg]
C_max_LFP = 22800  # [mol/m^3] ,found from paper Representative Volume... by Ali Ghorbani Kashkooli (osti.gov)
el_cond_ca = 100000  # [S/m]
C_dl_ca = 0.01  # [F/m^2]
R_sei_ca = 0  # [ohm m^2], solid electrolyte interface resistance, 0.012 <-- BDS
i_o_ca = 0.2  # [A/m^2], exchange current density
n_ca = -1  # [equivalents], charge transferred to the electrode
beta_ca = 0.5  # [-], Symmetry factor
D_o_ca = 1e-14  # [m^2/s], Solid-state diffusion coefficient
#   - Electrolyte
D_0_Li_ion = 1.75e-10  # [m^2/s], Li^+ liquid diffusion coefficient


# Numerical Parameters
N_part  = 3  # Number of nodes in solid particle
N_an    = 3  # Number of nodes in the anode
N_sep   = 1  # Number of nodes in the separator
N_ca    = 3  # Number of nodes in the cathode
