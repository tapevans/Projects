import numpy as np
from math import exp

# Constants
F = 96485
R = 8.3145


def residual(t, SV, pars, ptr):
    # Initialize dSV_dt
    SV = np.reshape(SV, (pars.N_tot, -1))
    print('t=', t)
    print('SV=',SV)
    dSV_dt = np.zeros_like(SV)
    #print(dSV_dt)

    # Anode
    for i in range(pars.N_an):
        block = i

        # Charge Transfer
        X_an_avg = np.dot(SV[block][ptr.X_Li_ptr:ptr.X_Li_ptr+pars.N_part], pars.V_part_an)/np.sum(pars.V_part_an)
        #print('X_an_avg', X_an_avg)
        phi_ed_equil = VoltEquibAnode(X_an_avg)
        #print('phi_ed_equil', phi_ed_equil)
        phi_ed_surf = VoltEquibAnode(SV[block][ptr.X_Li_ptr])
        #print('phi_ed_surf', phi_ed_surf)
        phi_elyte = phi_ed_surf - SV[block][ptr.delta_phi_ptr]
        #print('phi_elyte', phi_elyte)
        eta = (phi_ed_surf - phi_elyte) - (phi_ed_equil - pars.phi_elyte_equil)
        print('eta_an', eta)
        RTinv = 1 / R / SV[block][ptr.T_ptr]
        #print('T',SV[block][ptr.T_ptr])
        i_Far = pars.i_o_an * (exp(-pars.n_an * F * pars.beta_an * eta * RTinv) - exp(pars.n_an * F * (1 - pars.beta_an) * eta * RTinv))
        #print('i_Far_an', i_Far)
        s_dot_Li = -pars.nu_Li * i_Far / pars.n_an / F
        print('s_dot_Li_an', s_dot_Li)
        s_dot_Li_ion = -pars.nu_Li_ion * i_Far / pars.n_an / F
        print('s_dot_Li_ion_an', s_dot_Li_ion)
        i_dl = pars.i_ext * (pars.A_geo_min / pars.A_surf_an) - i_Far
        #print('i_dl', i_dl)
            #A_surf[block] use this if start making particles different sizes throughout the electrode

        # Mass Transport
        #   - Li^+ diffusion coming in from the left block
        if block == 0: # If current block is the first anode block (boundary condition)
            N_dot_left = 0
        else:
            D_left = (1/(pars.del_z_vec[block-1] + pars.del_z_vec[block])) * (pars.del_z_vec[block-1]*pars.D_eff_Li_ion_vec[block]+ pars.del_z_vec[block]*pars.D_eff_Li_ion_vec[block])
            gradC_left = (SV[block-1][ptr.X_Li_ion_ptr] - SV[block][ptr.X_Li_ion_ptr])/np.absolute(pars.pos_vec[block-1]-pars.pos_vec[block])
            C_intf_left = SV[block-1][ptr.X_Li_ion_ptr] + gradC_left * pars.del_z_vec[block-1] * 0.5
            phi_elyte_left = VoltEquibAnode(SV[block-1][ptr.X_Li_ptr]) - SV[block-1][ptr.delta_phi_ptr]
            #gradphi_left = (phi_elyte_left - phi_elyte)/np.absolute(pars.pos_vec[block-1]-pars.pos_vec[block])
            gradphi_left = 0
            N_dot_left = -(- D_left * gradC_left - D_left * C_intf_left * F * RTinv * gradphi_left)  # Changed signs

        #   - Li^+ diffusion coming in from the right block
        D_right = (1 / (pars.del_z_vec[block + 1] + pars.del_z_vec[block])) * (pars.del_z_vec[block + 1] * pars.D_eff_Li_ion_vec[block] + pars.del_z_vec[block] *pars.D_eff_Li_ion_vec[block])
        gradC_right = (SV[block + 1][ptr.X_Li_ion_ptr] - SV[block][ptr.X_Li_ion_ptr]) / np.absolute(pars.pos_vec[block + 1] - pars.pos_vec[block])
        #print('gradC_right', gradC_right)
        C_intf_right = SV[block][ptr.X_Li_ion_ptr] + gradC_right * pars.del_z_vec[block] * 0.5
        phi_elyte_right = VoltEquibAnode(SV[block + 1][ptr.X_Li_ptr]) - SV[block + 1][ptr.delta_phi_ptr]
        #print('phi_elyte_right', phi_elyte_right)
        #gradphi_right = (phi_elyte_right - phi_elyte) / np.absolute(pars.pos_vec[block + 1] - pars.pos_vec[block])
        gradphi_right = 0
        #print('gradphi_right', gradphi_right)
        N_dot_right = -(- D_right * gradC_right - D_right * C_intf_right * F * RTinv * gradphi_right)  # Changed signs

        #print('N_dot_left', N_dot_left)
        #print('N_dot_right', N_dot_right)
        N_dot_Li_ion = N_dot_left + N_dot_right
        #print('N_dot_Li', N_dot_Li_ion)
        # Governing Equations
        #   - Li^+ mole fraction
        dSV_dt[block][ptr.X_Li_ion_ptr] = (1 / (pars.C_Li_ion_ref * (1 - pars.eps_an) * pars.del_z_vec[block])) * ((pars.A_surf_an / pars.A_geo_min) + N_dot_Li_ion)
        #   - Double layer voltage
        dSV_dt[block][ptr.delta_phi_ptr] = i_dl*pars.C_dl_an_inv
        #   - Li mole fraction (Solid state diffusion)
        for j in range(pars.N_part):
            #print('j', j)
            if j == 0:
                gradC_inside = pars.C_max_an * (SV[block][ptr.X_Li_ptr+1] - SV[block][ptr.X_Li_ptr]) * pars.del_r_an_inv
                #print('gradC_inside', gradC_inside)
                gradphi_inside = (VoltEquibAnode(SV[block][ptr.X_Li_ptr + 1]) - VoltEquibAnode(SV[block][ptr.X_Li_ptr])) * pars.del_r_an_inv
                #print('gradphi_inside', gradphi_inside)
                C_intf_inside = pars.C_max_an * (SV[block][ptr.X_Li_ptr+1] + SV[block][ptr.X_Li_ptr]) * 0.5
                N_dot_inside = -((- pars.D_ss_an * gradC_inside ) * pars.SA_part_an[j+1])
                #print('N_dot_inside', N_dot_inside)
                N_dot_outside = s_dot_Li * pars.A_surf_an
                #print('N_dot_outside', N_dot_outside)
            elif j == (pars.N_part-1):
                #print('elif')
                N_dot_inside = 0

                gradC_outside = pars.C_max_an * (SV[block][(ptr.X_Li_ptr+j)-1] - SV[block][(ptr.X_Li_ptr+j)]) * pars.del_r_an_inv
                #print('gradC_outside', gradC_outside)
                gradphi_outside = (VoltEquibAnode(SV[block][(ptr.X_Li_ptr+j) - 1]) - VoltEquibAnode(SV[block][(ptr.X_Li_ptr+j)])) * pars.del_r_an_inv
                #print('gradphi_outside', gradphi_outside)
                C_intf_outside = pars.C_max_an * (SV[block][(ptr.X_Li_ptr+j)-1] + SV[block][(ptr.X_Li_ptr+j)]) * 0.5
                N_dot_outside = -((- pars.D_ss_an * gradC_outside ) * pars.SA_part_an[j])
                #print('N_dot_outside', N_dot_outside)
            else:
                gradC_inside = pars.C_max_an * (SV[block][(ptr.X_Li_ptr+j) + 1] - SV[block][(ptr.X_Li_ptr+j)]) * pars.del_r_an_inv
                #print('gradC_inside', gradC_inside)
                gradphi_inside = (VoltEquibAnode(SV[block][(ptr.X_Li_ptr+j) + 1]) - VoltEquibAnode(SV[block][(ptr.X_Li_ptr+j)])) * pars.del_r_an_inv
                #print('gradphi_inside', gradphi_inside)
                C_intf_inside = pars.C_max_an * (SV[block][(ptr.X_Li_ptr+j) + 1] + SV[block][(ptr.X_Li_ptr+j)]) * 0.5
                N_dot_inside = -((- pars.D_ss_an * gradC_inside ) * pars.SA_part_an[j+1])
                #print('N_dot_inside', N_dot_inside)

                gradC_outside = pars.C_max_an * (SV[block][(ptr.X_Li_ptr+j) - 1] - SV[block][(ptr.X_Li_ptr+j)]) * pars.del_r_an_inv
                #print('gradC_outside', gradC_outside)
                gradphi_outside = (VoltEquibAnode(SV[block][(ptr.X_Li_ptr+j) - 1]) - VoltEquibAnode(SV[block][(ptr.X_Li_ptr+j)])) * pars.del_r_an_inv
                #print('gradphi_outside', gradphi_outside)
                C_intf_outside = pars.C_max_an * (SV[block][(ptr.X_Li_ptr+j) - 1] + SV[block][(ptr.X_Li_ptr+j)]) * 0.5
                N_dot_outside = -((- pars.D_ss_an * gradC_outside ) * pars.SA_part_an[j])
                #print('N_dot_outside', N_dot_outside)

            dSV_dt[block][ptr.X_Li_ptr + j] = (N_dot_inside + N_dot_outside)/(pars.C_max_an * pars.V_part_an[j])

        #   - Temperature
        dSV_dt[block][ptr.T_ptr] = 0 # Temperature, Currently for all blocks temperature doesn't change

    # Separator
    for i in range(pars.N_sep):
        block = i + pars.N_an

        # Mass Transport
        #   - Li^+ diffusion coming in from the left block
        D_left = (1 / (pars.del_z_vec[block - 1] + pars.del_z_vec[block])) * (pars.del_z_vec[block - 1] * pars.D_eff_Li_ion_vec[block] + pars.del_z_vec[block] *pars.D_eff_Li_ion_vec[block])
        gradC_left = (SV[block - 1][ptr.X_Li_ion_ptr] - SV[block][ptr.X_Li_ion_ptr]) / np.absolute(pars.pos_vec[block - 1] - pars.pos_vec[block])
        C_intf_left = SV[block - 1][ptr.X_Li_ion_ptr] + gradC_left * pars.del_z_vec[block - 1] * 0.5
        phi_elyte_left = VoltEquibAnode(SV[block - 1][ptr.X_Li_ptr]) - SV[block - 1][ptr.delta_phi_ptr]
        #gradphi_left = (phi_elyte_left - phi_elyte) / np.absolute(pars.pos_vec[block - 1] - pars.pos_vec[block])
        gradphi_left = 0
        N_dot_left = -(- D_left * gradC_left - D_left * C_intf_left * F * RTinv * gradphi_left)

        #   - Li^+ diffusion coming in from the right block
        D_right = (1 / (pars.del_z_vec[block + 1] + pars.del_z_vec[block])) * (pars.del_z_vec[block + 1] * pars.D_eff_Li_ion_vec[block] + pars.del_z_vec[block] *pars.D_eff_Li_ion_vec[block])
        gradC_right = (SV[block + 1][ptr.X_Li_ion_ptr] - SV[block][ptr.X_Li_ion_ptr]) / np.absolute(pars.pos_vec[block + 1] - pars.pos_vec[block])
        C_intf_right = SV[block][ptr.X_Li_ion_ptr] + gradC_right * pars.del_z_vec[block] * 0.5
        phi_elyte_right = VoltEquibAnode(SV[block + 1][ptr.X_Li_ptr]) - SV[block + 1][ptr.delta_phi_ptr]
        #gradphi_right = (phi_elyte_right - phi_elyte) / np.absolute(pars.pos_vec[block + 1] - pars.pos_vec[block])
        gradphi_right = 0
        N_dot_right = -(- D_right * gradC_right - D_right * C_intf_right * F * RTinv * gradphi_right)

        N_dot_Li_ion = N_dot_left + N_dot_right

        # Governing Equations
        #   - Li^+ mole fraction
        dSV_dt[block][ptr.X_Li_ion_ptr] = (1 / (pars.C_Li_ion_ref * (1 - pars.eps_sep) * pars.del_z_vec[block])) * (N_dot_Li_ion)

        #   - Double Layer Voltage
        dSV_dt[block][ptr.delta_phi_ptr] = 0

        #   - Li mole fraction
        for j in range(pars.N_part):
            dSV_dt[block][ptr.X_Li_ptr + j] = 0

        #   - Temperature
        dSV_dt[block][ptr.T_ptr] = 0 # Temperature, Currently for all blocks temperature doesn't change

    # Cathode
    for i in range(pars.N_ca):
        block = i + pars.N_an + pars.N_sep

        # Charge Transfer
        phi_ed_equil = VoltEquibCathode(np.dot(SV[block][ptr.X_Li_ptr:ptr.X_Li_ptr + pars.N_part], pars.V_part_ca) / np.sum(pars.V_part_ca))
        #print(np.dot(SV[block][ptr.X_Li_ptr:ptr.X_Li_ptr + pars.N_part], pars.V_part_ca))
        #print('phi_ed_equil',phi_ed_equil)
        phi_ed_surf = VoltEquibCathode(SV[block][ptr.X_Li_ptr])
        phi_elyte = phi_ed_surf - SV[block][ptr.delta_phi_ptr]
        eta = (phi_ed_surf - phi_elyte) - (phi_ed_equil - pars.phi_elyte_equil)
        print('eta_ca', eta)
        RTinv = 1 / R / SV[block][ptr.T_ptr]
        i_Far = pars.i_o_ca * (exp(-pars.n_ca * F * pars.beta_ca * eta * RTinv) - exp(pars.n_ca * F * (1 - pars.beta_ca) * eta * RTinv))
        #print('i_Far_ca', i_Far)
        s_dot_Li = -pars.nu_Li * i_Far / pars.n_ca / F
        print('s_dot_Li_ca', s_dot_Li)
        s_dot_Li_ion = -pars.nu_Li_ion * i_Far / pars.n_ca / F
        print('s_dot_Li_ca', s_dot_Li)
        i_dl =  -pars.i_ext * (pars.A_geo_min / pars.A_surf_ca) - i_Far #%%%%%%%%%%%% maybe this is -i_ext
        # A_surf[block] use this if start making particles different sizes throughout the electrode

        # Mass Transport
        #   - Li^+ diffusion coming in from the left block
        D_left = (1 / (pars.del_z_vec[block - 1] + pars.del_z_vec[block])) * (pars.del_z_vec[block - 1] * pars.D_eff_Li_ion_vec[block] + pars.del_z_vec[block] * pars.D_eff_Li_ion_vec[block])
        gradC_left = (SV[block - 1][ptr.X_Li_ion_ptr] - SV[block][ptr.X_Li_ion_ptr]) / np.absolute(pars.pos_vec[block - 1] - pars.pos_vec[block])
        C_intf_left = SV[block - 1][ptr.X_Li_ion_ptr] + gradC_left * pars.del_z_vec[block - 1] * 0.5
        phi_elyte_left = VoltEquibCathode(SV[block - 1][ptr.X_Li_ptr]) - SV[block - 1][ptr.delta_phi_ptr]
        #gradphi_left = (phi_elyte_left - phi_elyte) / np.absolute(pars.pos_vec[block - 1] - pars.pos_vec[block])
        gradphi_left = 0
        N_dot_left = -(- D_left * gradC_left - D_left * C_intf_left * F * RTinv * gradphi_left)

        #   - Li^+ diffusion coming in from the right block
        if block == pars.N_tot-1:  # If current block is the first anode block (boundary condition)
            print('in if statement')
            N_dot_right = 0
        else:
            D_right = (1 / (pars.del_z_vec[block + 1] + pars.del_z_vec[block])) * (pars.del_z_vec[block + 1] * pars.D_eff_Li_ion_vec[block] + pars.del_z_vec[block] * pars.D_eff_Li_ion_vec[block])
            gradC_right = (SV[block + 1][ptr.X_Li_ion_ptr] - SV[block][ptr.X_Li_ion_ptr]) / np.absolute(pars.pos_vec[block + 1] - pars.pos_vec[block])
            C_intf_right = SV[block][ptr.X_Li_ion_ptr] + gradC_right * pars.del_z_vec[block] * 0.5
            phi_elyte_right = VoltEquibCathode(SV[block + 1][ptr.X_Li_ptr]) - SV[block + 1][ptr.delta_phi_ptr]
            #gradphi_right = (phi_elyte_right - phi_elyte) / np.absolute(pars.pos_vec[block + 1] - pars.pos_vec[block])
            gradphi_right = 0
            N_dot_right = -(- D_right * gradC_right - D_right * C_intf_right * F * RTinv * gradphi_right)

        N_dot_Li_ion = N_dot_left + N_dot_right

        # Governing Equations
        #   - Li^+ mole fraction
        dSV_dt[block][ptr.X_Li_ion_ptr] = (1 / (pars.C_Li_ion_ref * (1 - pars.eps_ca) * pars.del_z_vec[block])) * ((pars.A_surf_ca / pars.A_geo_min) + N_dot_Li_ion)

        #   - Double layer voltage
        dSV_dt[block][ptr.delta_phi_ptr] = i_dl * pars.C_dl_an_inv

        #   - Li mole fraction (Solid state diffusion)
        for j in range(pars.N_part):
            # Surface BC
            if j == 0:
                gradC_inside = pars.C_max_ca * (SV[block][ptr.X_Li_ptr + 1] - SV[block][ptr.X_Li_ptr]) * pars.del_r_ca_inv
                gradphi_inside = (VoltEquibCathode(SV[block][ptr.X_Li_ptr + 1]) - VoltEquibCathode(SV[block][ptr.X_Li_ptr])) * pars.del_r_ca_inv
                C_intf_inside = pars.C_max_ca * (SV[block][ptr.X_Li_ptr + 1] + SV[block][ptr.X_Li_ptr]) * 0.5
                N_dot_inside = -((- pars.D_ss_ca * gradC_inside ) * pars.SA_part_ca[j + 1])

                N_dot_outside = s_dot_Li * pars.A_surf_ca
            # Center BC
            elif j == (pars.N_part - 1):
                N_dot_inside = 0

                gradC_outside = pars.C_max_ca * (SV[block][(ptr.X_Li_ptr+j) - 1] - SV[block][(ptr.X_Li_ptr+j)]) * pars.del_r_ca_inv
                gradphi_outside = (VoltEquibCathode(SV[block][(ptr.X_Li_ptr+j) - 1]) - VoltEquibCathode(SV[block][(ptr.X_Li_ptr+j)])) * pars.del_r_ca_inv
                C_intf_outside = pars.C_max_ca * (SV[block][(ptr.X_Li_ptr+j) - 1] + SV[block][(ptr.X_Li_ptr+j)]) * 0.5
                N_dot_outside = -((- pars.D_ss_ca * gradC_outside ) * pars.SA_part_ca[j])
            else:
                gradC_inside = pars.C_max_ca * (SV[block][(ptr.X_Li_ptr+j) + 1] - SV[block][(ptr.X_Li_ptr+j)]) * pars.del_r_ca_inv
                gradphi_inside = (VoltEquibCathode(SV[block][(ptr.X_Li_ptr+j) + 1]) - VoltEquibCathode(SV[block][(ptr.X_Li_ptr+j)])) * pars.del_r_ca_inv
                C_intf_inside = pars.C_max_ca * (SV[block][(ptr.X_Li_ptr+j) + 1] + SV[block][(ptr.X_Li_ptr+j)]) * 0.5
                N_dot_inside = -((- pars.D_ss_ca * gradC_inside ) * pars.SA_part_ca[j + 1])

                gradC_outside = pars.C_max_ca * (SV[block][(ptr.X_Li_ptr+j) - 1] - SV[block][(ptr.X_Li_ptr+j)]) * pars.del_r_ca_inv
                gradphi_outside = (VoltEquibCathode(SV[block][(ptr.X_Li_ptr+j) - 1]) - VoltEquibCathode(SV[block][(ptr.X_Li_ptr+j)])) * pars.del_r_ca_inv
                C_intf_outside = pars.C_max_ca * (SV[block][(ptr.X_Li_ptr+j) - 1] + SV[block][(ptr.X_Li_ptr+j)]) * 0.5
                N_dot_outside = -((- pars.D_ss_ca * gradC_outside) * pars.SA_part_ca[j])

            dSV_dt[block][ptr.X_Li_ptr + j] = (N_dot_inside + N_dot_outside) / (pars.C_max_ca * pars.V_part_ca[j])

        #   - Temperature
        dSV_dt[block][ptr.T_ptr] = 0  # Temperature, Currently for all blocks temperature doesn't change

    print('dSV_dt', dSV_dt)
    dSV_dt = dSV_dt.flatten()
    #pbreak
    return dSV_dt


def DiffCoeffAnode(X):
    # Change to be a function of concentration
    # Currently Graphite
    D_o = 1e-16
    #print('Inside Fn', D_o)
    return D_o


def DiffCoeffCathode(X):
    # Change to be a function of concentration
    # Currently LiFePO_4
    D_o = 1e-14
    return D_o


def ElyteConductivity(T):
    # Change to be a function of temperature
    # Current salt is LiPF
    sigma = 0.1
    return sigma


def ExchangeCurrentDensityAnode(T):
    # Change to be a function of temperature
    # Currently Graphite
    i_o = 0.2
    return i_o


def ExchangeCurrentDensityCathode(T):
    # Change to be a function of temperature
    # Currently LiFePO_4
    i_o = 0.2
    return i_o


def VoltEquibCathode(X):
    # Returns the equilibrium voltage value that corresponses to lithiation
    # Currently LiFePO_4, using Li ref electrode
    # Concentration
    xp = [0, 2.22E-15, 1.00E-06, 0.00802, 0.016, 0.024, 0.0321, 0.0401, 0.0481, 0.0561, 0.0641, 0.0721, 0.0802, 0.0882,
          0.0962, 0.104, 0.112, 0.120, 0.128, 0.136, 0.144, 0.152, 0.160, 0.168, 0.176, 0.184, 0.192, 0.200, 0.208,
          0.216,  0.224, 0.232, 0.240, 0.248, 0.257, 0.265, 0.273, 0.281, 0.289, 0.297, 0.305, 0.313, 0.321, 0.329,
          0.337,  0.345, 0.353, 0.361, 0.369, 0.377, 0.385, 0.393, 0.401, 0.409, 0.417, 0.425, 0.433, 0.441, 0.449,
          0.457,  0.465, 0.473, 0.481, 0.489, 0.497, 0.505, 0.513, 0.521, 0.529, 0.537, 0.545, 0.553, 0.561, 0.569,
          0.577,  0.585, 0.593, 0.601, 0.609, 0.617, 0.625, 0.633, 0.641, 0.649, 0.657, 0.665, 0.673, 0.681, 0.689,
          0.697,  0.705, 0.713, 0.721, 0.729, 0.737, 0.745, 0.754, 0.762, 0.770, 0.778, 0.786, 0.794, 0.802, 0.810,
          0.818,  0.826, 0.834, 0.842, 0.850, 0.858, 0.866, 0.874, 0.882, 0.890, 0.898, 0.906, 0.914, 0.922, 0.930,
          0.938,  0.946, 0.954, 0.962, 0.970, 0.978, 0.986, 0.994, 1]
    # Voltage
    fp = [6.00, 3.75, 3.75, 3.48, 3.42, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41,
          3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41,
          3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40,
          3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40,
          3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.39, 3.39, 3.39, 3.39, 3.39, 3.39, 3.39, 3.39,
          3.39, 3.39, 3.39, 3.39, 3.39, 3.39, 3.38, 3.38, 3.38, 3.38, 3.38, 3.38, 3.38, 3.37, 3.37, 3.37, 3.37, 3.37,
          3.36, 3.36, 3.36, 3.35, 3.35, 3.34, 3.34, 3.33, 3.32, 3.31, 3.29, 3.27, 3.24, 3.18, 3.04, 2.90, 2.78, 2.61,
          2.38, 2.38]

    phi = np.interp(X, xp, fp)
    return phi


def VoltEquibAnode(X):
    # Returns the equilibrium voltage value that corresponses to lithiation
    # Currently Graphite, using Li ref electrode
    # Concentration
    xp = [0.0000000000, 0.0001014930, 0.0002346440, 0.0003677950, 0.0005009460, 0.0006368900, 0.0008436000,
          0.0009804760, 0.0011145600, 0.0012477100, 0.0013827200, 0.0015158700, 0.0016518200, 0.0018576000,
          0.0019944700, 0.0022617100, 0.0025289400, 0.0028399400, 0.0030745800, 0.0033501900, 0.0036546700,
          0.0039265600, 0.0042608300, 0.0045699700, 0.0049089000, 0.0052459700, 0.0055839600, 0.0059852800,
          0.0064266300, 0.0067646300, 0.0071342900, 0.0075719200, 0.0080449300, 0.0084527600, 0.0088205600,
          0.0092926400, 0.0097284100, 0.0102061000, 0.0107396000, 0.0112797000, 0.0117489000, 0.0122909000,
          0.0126298000, 0.0132741000, 0.0137751000, 0.0143170000, 0.0148580000, 0.0154306000, 0.0160061000,
          0.0166848000, 0.0172929000, 0.0179009000, 0.0184391000, 0.0189484000, 0.0194950000, 0.0202017000,
          0.0207492000, 0.0214271000, 0.0222697000, 0.0229066000, 0.0238787000, 0.0246264000, 0.0255026000,
          0.0263173000, 0.0273956000, 0.0280762000, 0.0289180000, 0.0301648000, 0.0312449000, 0.0321909000,
          0.0333408000, 0.0345979000, 0.0356091000, 0.0368530000, 0.0381399000, 0.0392246000, 0.0407740000,
          0.0421269000, 0.0433067000, 0.0444557000, 0.0460805000, 0.0479725000, 0.0495583000, 0.0515527000,
          0.0535137000, 0.0546906000, 0.0565855000, 0.0589878000, 0.0612821000, 0.0635819000, 0.0656705000,
          0.0681752000, 0.0694601000, 0.0720915000, 0.0751586000, 0.0781718000, 0.0816840000, 0.0854085000,
          0.0891330000, 0.0928575000, 0.0993754000, 0.1049620000, 0.1356890000, 0.1636230000, 0.1831770000,
          0.1869010000, 0.1971440000, 0.2129730000, 0.2288020000, 0.2437000000, 0.2651160000, 0.3936110000,
          0.4718260000, 0.4783440000, 0.5155890000, 0.5565580000, 0.7148500000, 0.8480010000, 0.8526560000,
          0.8982820000, 0.9271460000, 0.9513560000, 0.9597360000, 0.9718410000, 0.9876700000, 1.0000000000]
    # Voltage
    fp = [6.000000, 2.330000, 2.130000, 1.990000, 1.880000, 1.800000, 1.690000, 1.630000, 1.580000, 1.540000, 1.500000,
          1.460000, 1.430000, 1.400000, 1.380000, 1.370000, 1.350000, 1.330000, 1.320000, 1.300000, 1.290000, 1.270000,
          1.250000, 1.240000, 1.220000, 1.210000, 1.190000, 1.180000, 1.160000, 1.150000, 1.140000, 1.120000, 1.110000,
          1.090000, 1.080000, 1.070000, 1.060000, 1.040000, 1.020000, 1.010000, 0.995020, 0.982180, 0.967190, 0.951900,
          0.938440, 0.924380, 0.912140, 0.899300, 0.880030, 0.865050, 0.850980, 0.837520, 0.825590, 0.809080, 0.796850,
          0.784610, 0.770550, 0.758010, 0.745470, 0.729570, 0.714890, 0.703260, 0.690730, 0.678800, 0.665950, 0.652800,
          0.641180, 0.629560, 0.617630, 0.605710, 0.594080, 0.582160, 0.568700, 0.557690, 0.546380, 0.534750, 0.523740,
          0.512120, 0.498060, 0.488880, 0.477560, 0.465940, 0.454020, 0.442700, 0.431690, 0.419460, 0.407840, 0.396220,
          0.384900, 0.372360, 0.361050, 0.349730, 0.340250, 0.328020, 0.316390, 0.305380, 0.293460, 0.279390, 0.264400,
          0.251250, 0.234740, 0.223120, 0.214550, 0.203540, 0.192540, 0.184580, 0.175100, 0.163790, 0.152780, 0.141460,
          0.128920, 0.118530, 0.107520, 0.106900, 0.097730, 0.086720, 0.077240, 0.067450, 0.066530, 0.061030, 0.050020,
          0.038700, 0.028000, 0.016380, 0.005060, 0.005060]
    phi = np.interp(X, xp, fp)
    return phi


def VoltElyte(X, T):
    # Returns the voltage of the electrolyte based on the Li^+ concentration and temperature
    # Currently salt is LiPF_6
    phi = 2.5
    return phi


def SOCtoX(SOC):
    # Returns the mole fraction for both the anode and cathode at a given SOC
    # Max and min values are calculated from voltage limits and formation of the battery %%%%%%%%%%%%% But how?
    X_max_an = 1.0
    X_min_an = 0.063
    delta_X_an = X_max_an - X_min_an
    X_an = X_min_an + (delta_X_an * (SOC / 100))

    X_max_ca = 0.985
    X_min_ca = 0.106
    delta_X_ca = X_max_ca - X_min_ca
    X_ca = X_max_ca - (delta_X_ca * (SOC / 100))

    return X_an, X_ca