from A123_P2D_inputs import *
import numpy as np
import math
from A123_P2D_function import SOCtoX
from A123_P2D_function import VoltEquibAnode
from A123_P2D_function import VoltEquibCathode
#a = [[1, 2, 4], [3, 7, 8]]
#file = open("test.txt", "w")
#file.write(str(a))
#file.close()

# If multiple simulation implemented then create a loop that creates a matrix of important sim info
#   - Add those calcs here (Future or maybe have an ext file loop through parameters)

# Calculate initial parameters from the inputs
#   - SOC --> X of anode and cathode
X_an, X_ca = SOCtoX(SOC_start)
#   - Concentration of the electrolyte
X_elyte = 1  # [-] %%%%%%%%%%%%%%%%%%%%
C_Li_ion_ref = 1  # [M]
#   - Voltage of outer volume of the electrode
phi_an = VoltEquibAnode(X_an)  # V
phi_ca = VoltEquibCathode(X_ca)  # V
#   - phi_elyte
phi_elyte = 1.0  # V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   - temperature (Just use the variable imported)

# Initialize State Variables (SV) array
#   - Initialize component block arrays
Np_ones = np.ones(N_part)
X_an_vec = X_an*Np_ones
X_ca_vec = X_ca*Np_ones
q_0 = 0
#anode_block = np.concatenate([[X_elyte], [phi_an - phi_elyte], X_an_vec, [Temp_start]])
anode_block = np.concatenate([[X_elyte], [q_0], X_an_vec, [Temp_start]])
#sep_block = np.concatenate([[X_elyte], [0], 0*X_an_vec, [Temp_start]])
sep_block = np.concatenate([[X_elyte], [0], 0*X_an_vec, [Temp_start]])
#cathode_block = np.concatenate([[X_elyte], [phi_ca - phi_elyte], X_ca_vec, [Temp_start]])
cathode_block = np.concatenate([[X_elyte], [q_0], X_ca_vec, [Temp_start]])

#   - Initialize the SV_O using the first anode_block
N_tot = N_an + N_sep + N_ca
SV_0 = np.zeros((N_tot, len(anode_block)))
for i in range(N_tot):
    if i < (N_an + N_sep):
        if i < N_an:
            SV_0[i][:] = anode_block
        else:
            SV_0[i][:] = sep_block
    else:
        SV_0[i][:] = cathode_block

    #print(SV_0)
    #SV_0 = np.reshape(SV_0, (1, -1))
SV_0 = SV_0.flatten()
    #print(SV_0)

# Determine Simulation run time
A_geo_an = del_x_an * del_y_an
A_geo_ca = del_x_ca * del_y_ca
A_geo_min = min(A_geo_an, A_geo_ca)

capacity_an = capacity_graphite*(H_an*A_geo_an)*eps_an*density_an  # [Ahr]
capacity_ca = capacity_LFP*(H_ca*A_geo_ca)*eps_ca*density_ca  # [Ahr]
capacity = min(capacity_an, capacity_ca)

I_ext = capacity*C_rate*ChargeOrDischarge # This would change if doing EIS simulation
#print('I_ext', I_ext)
t_final = charge_frac*3600./C_rate



# Create vector of dVol for electrode solid particles
N_tot_ss_part_an = eps_an*(del_x_an*del_y_an*H_an)/((4/3)*math.pi*r_p_an**3)
N_tot_ss_part_ca = eps_ca*(del_x_ca*del_y_ca*H_ca)/((4/3)*math.pi*r_p_ca**3)
N_ss_an_block = N_tot_ss_part_an/N_an
N_ss_ca_block = N_tot_ss_part_ca/N_ca

delta_r_an = r_p_an/N_part
delta_r_ca = r_p_ca/N_part
r_outer_an = np.zeros(N_part+1)
r_outer_ca = np.zeros(N_part+1)
    #print('delta_r_an', delta_r_ca)
for i in range(N_part):
    r_outer_an[i] = r_p_an - delta_r_an * i
    r_outer_ca[i] = r_p_ca - delta_r_ca * i
    #print('r_outer_an', r_outer_ca)
#%%%%%%%%%%%%%%%%%%%%%%May try to compress this to not use a for loop and won't need to initialize
#SA_part_an = np.zeros(N_part)
#V_part_an = np.zeros(N_part)
#SA_part_ca = np.zeros(N_part)
#V_part_ca = np.zeros(N_part)
#for i in range(N_part):
#    SA_part_an[i] = 4 * math.pi * r_outer_an[i] ** 2
#    SA_part_ca[i] = 4 * math.pi * r_outer_ca[i] ** 2
#    V_part_an[i] = (4 / 3) * math.pi * (r_outer_an[i] ** 3 - r_outer_an[i+1] ** 3)
#    V_part_ca[i] = (4 / 3) * math.pi * (r_outer_ca[i] ** 3 - r_outer_ca[i + 1] ** 3)
    #print('SA_part_an', SA_part_ca)
    #print('V_part_an', V_part_ca)
SA_part_an = (4 * math.pi * r_outer_an[0:N_part] ** 2) * N_ss_an_block
V_part_an = ((4 / 3) * math.pi * (r_outer_an[0:N_part] ** 3 - r_outer_an[1:N_part+1] ** 3)) * N_ss_an_block
SA_part_ca = (4 * math.pi * r_outer_ca[0:N_part] ** 2) * N_ss_ca_block
V_part_ca = ((4 / 3) * math.pi * (r_outer_ca[0:N_part] ** 3 - r_outer_ca[1:N_part+1] ** 3)) * N_ss_ca_block
#print(V_part_an_test)
#print(V_part_an)
#a = np.array([1 , 2, 3])
#print('a', 3*a)
# Create block position, thickness, and liquid diffusion vector
#   - Thickness calculations
del_z_an = H_an/N_an
del_z_ca = H_ca/N_ca
del_z_sep = H_sep/N_sep
#   - Tortuosity calculations
tau_an = tau_prefac_an*eps_an**(1-bruggeman_an)
tau_ca = tau_prefac_ca*eps_ca**(1-bruggeman_ca)
tau_sep = tau_prefac_sep*eps_sep**(1-bruggeman_sep)
#   - Effective diffusion coefficient calculation
D_eff_li_ion_an = ((1-eps_an)/tau_an) * D_0_Li_ion
D_eff_li_ion_ca = ((1-eps_ca)/tau_ca) * D_0_Li_ion
D_eff_li_ion_sep = ((1-eps_sep)/tau_sep) * D_0_Li_ion
#   - Initialize vectors
pos_vec = np.zeros(N_tot)
del_z_vec = np.zeros(N_tot)
D_eff_Li_ion_vec = np.zeros(N_tot)

for i in range(N_tot):
    if i < (N_an + N_sep):
        if i < N_an:
            pos_vec[i] = i*del_z_an + 0.5*del_z_an
            del_z_vec[i] = del_z_an
            D_eff_Li_ion_vec[i] = D_eff_li_ion_an
        else:
            pos_vec[i] = (i-N_an) * del_z_sep + 0.5 * del_z_sep + H_an
            del_z_vec[i] = del_z_sep
            D_eff_Li_ion_vec[i] = D_eff_li_ion_sep
    else:
        pos_vec[i] = (i - N_an - N_sep) * del_z_ca + 0.5 * del_z_ca + H_an + H_sep
        del_z_vec[i] = del_z_ca
        D_eff_Li_ion_vec[i] = D_eff_li_ion_ca
#print(pos_vec)
#print(del_z_vec)
#print(D_eff_Li_ion_vec)
# Create block particle surface area vector
# %%%%%%%Furture implementation


# Create class to point to the correct variable location in the SV
class ptr:
    X_Li_ion_ptr = 0
    #delta_phi_ptr = X_Li_ion_ptr + 1
    charge_dl_ptr = X_Li_ion_ptr + 1
    X_Li_ptr = delta_phi_ptr + 1
    T_ptr = X_Li_ptr + N_part


# Load inputs and other parameters into 'pars' class:
class pars:
    T = Temp_start  # This could change if temperature dependence is implemented

    i_ext = I_ext/A_geo_min  # This could change if doing EIS simulation and not be calculated here
    N_part = N_part
    N_an = N_an
    N_sep = N_sep
    N_ca = N_ca
    N_tot = N_tot
    nu_Li = -1
    nu_Li_ion = 1
    A_geo_min = A_geo_min
    pos_vec = pos_vec
    del_z_vec = del_z_vec
    D_eff_Li_ion_vec = D_eff_Li_ion_vec

    # Anode
    C_dl_an_inv = 1 / C_dl_an
    R_sei_an = R_sei_an
    i_o_an = i_o_an
    n_an = n_an
    beta_an = beta_an
    tau_an = tau_an
    D_ss_an = D_o_an
    del_r_an_inv = 1 / delta_r_an
    SA_part_an = SA_part_an
    eps_an = eps_an
    V_part_an = V_part_an
    C_max_an = C_max_graphite
    A_geo_an = A_geo_an
    A_surf_an = SA_part_an[0]
    #A_fac_an = r_p_an / 3 / H_an / eps_an

    # Cathode
    C_dl_ca_inv = 1 / C_dl_ca
    R_sei_ca = R_sei_ca
    i_o_ca = i_o_ca
    n_ca = n_ca
    beta_ca = beta_ca
    tau_ca = tau_ca
    D_ss_ca = D_o_ca
    del_r_ca_inv = 1 / delta_r_ca
    SA_part_ca = SA_part_ca
    eps_ca = eps_ca
    V_part_ca = V_part_ca
    C_max_ca = C_max_LFP
    A_geo_ca = A_geo_ca
    A_surf_ca = SA_part_ca[0]
    #A_fac_ca = r_p_ca/3/H_ca/eps_ca

    # Separator
    phi_elyte_equil = phi_elyte
    eps_sep = eps_sep

    # Electrolyte
    C_Li_ion_ref = C_Li_ion_ref
