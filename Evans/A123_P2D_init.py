from A123_P2D_inputs import *
import numpy as np
import math

# If multiple simulation implemented then create a loop that creates a matrix of important sim info
#   - Add those calcs here (Future or maybe have an ext file loop through parameters)

# Calculate initial parameters from the inputs
#   - SOC --> X of anode and cathode
X_an = 0.2  # [-] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_ca = 0.9  # [-] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   - Concentration of the electrolyte
X_sep = 1  # [?]
#   - Voltage of outer volume of the electrode
phi_an = 2  # V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_ca = 4  # V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   - phi_elyte
phi_elyte = 2.5  # V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   - temperature (Just use the variable imported)

# Initialize State Variables (SV) array
#   - Initialize component block arrays
Np_ones = np.ones(N_part)

voltage_an = np.array([phi_an, phi_elyte])
X_an_vec = X_an*Np_ones
anode_block = np.concatenate([voltage_an, X_an_vec, [Temp_start]])
# print('anode_block = ', anode_block)
#voltage_sep = np.array([phi_elyte])
sep_block = np.array([phi_elyte, X_sep, Temp_start])
# print('sep_block = ', sep_block)
voltage_ca = np.array([phi_ca, phi_elyte])
X_ca_vec = X_ca*Np_ones
cathode_block = np.concatenate([voltage_ca, X_ca_vec, [Temp_start]])
# print('cathode_block = ', anode_block)

#   - Initialize the SV_O using the first anode_block
SV_0 = anode_block
#   - Add additional anode_blocks if N_an is greater than 1
if N_an > 1:
    # print('In an if statement')
    for i in range(N_an-1):
        # print('in anode for loop')
        SV_0 = np.concatenate((SV_0, anode_block))
#   - Add sep_blocks
for i in range(N_sep):
    # print('in sep for loop')
    SV_0 = np.concatenate((SV_0, sep_block))
#   - Add cathode_blocks
for i in range(N_ca):
    # print('in cathode for loop')
    SV_0 = np.concatenate((SV_0, cathode_block))

# print('SV_0', SV_0)

# Determine Simulation run time
H_an = del_z_an
H_ca = del_z_ca
capacity_an = capacity_graphite*H_an*eps_an*density_an
capacity_ca = capacity_LFP*H_ca*eps_ca*density_ca
capacity_area = min(capacity_an, capacity_ca)

t_final = charge_frac*3600./C_rate


# Create class to point to the correct variable location in the SV
class ptr:
    phi_an_ptr = []
    for i in range(N_an):
        phi_an_ptr.append(i*(3+N_part))   # index of the solid state anode voltage (outer volume of the particle)
    #print(phi_an_ptr)

    phi_elyte_an_ptr = phi_an_ptr + np.ones_like(phi_an_ptr)  # index for the voltage of the electrolyte in an anode_block
    #print(phi_elyte_an_ptr)

    X_an_ptr = phi_elyte_an_ptr + np.ones_like(phi_elyte_an_ptr)  # index for the concentration of the outer volume for
    #print(X_an_ptr)

    T_an_ptr = X_an_ptr + N_part*np.ones_like(X_an_ptr)  # index of the temperature for each block
    #print(T_an_ptr)
    #print(T_an_ptr[-1])

    phi_sep_ptr = [T_an_ptr[-1] + 1]
    #print('phi_sep_ptr',phi_sep_ptr)
    X_sep_ptr = [phi_sep_ptr[0] + 1]
    #print(X_sep_ptr)
    T_sep_ptr = [X_sep_ptr[0] + 1]
    #print(T_sep_ptr)

    # phi_ca_ptr = T_sep_ptr + [1]
    phi_ca_ptr = [T_sep_ptr[-1] + 1]
    # print(phi_ca_ptr)
    if N_ca > 1:
        for i in range(N_ca-1):
            phi_ca_ptr.append(phi_ca_ptr[-1] + (3+N_part))

    # print('phi_ca_ptr',phi_ca_ptr)
    phi_elyte_ca_ptr = phi_ca_ptr + np.ones_like(phi_ca_ptr)
    X_ca_ptr = phi_elyte_ca_ptr + np.ones_like(phi_elyte_ca_ptr)  # index for the concentration of the outer volume for
    T_ca_ptr = X_ca_ptr + N_part*np.ones_like(X_ca_ptr)  # index of the temperature for each block
    # print(T_ca_ptr)
    # print(np.shape(SV_0))


# Load inputs and other parameters into 'pars' class:
class pars:
    T = Temp_start  # This could change if temperature dependence is implemented

    i_ext = C_rate * capacity_area  # This could change if doing EIS simulation
    N_part = N_part
    # Anode
    C_dl_an_inv = 1 / C_dl_an
    R_sei_an = R_sei_an
    i_o_an = i_o_an
    n_an = n_an
    beta_an = beta_an
    tau_an = tau_an
    A_geo_an = del_x_an * del_y_an
    A_surf_an = 4*math.pi*r_p_an**2
    A_fac_an = r_p_an / 3 / H_an / eps_an

    # Cathode
    C_dl_ca_inv = 1 / C_dl_ca
    R_sei_ca = R_sei_ca
    i_o_ca = i_o_ca
    n_ca = n_ca
    beta_ca = beta_ca
    tau_ca = tau_ca
    A_geo_ca = del_x_ca * del_y_ca
    A_surf_ca = 4*math.pi*r_p_ca**2
    A_fac_ca = r_p_ca/3/H_ca/eps_ca
