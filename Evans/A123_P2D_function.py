import numpy as np
from math import exp

# Constants
F = 96485
R = 8.3145


def residual(t, SV, pars, ptr):
    # Old equations
    #RTinv = 1 / R / pars.T
    #eta_an = SV[0] - pars.dPhi_eq_an
    #i_Far_an = pars.i_o_an * (exp(-pars.n_an * F * pars.beta_an * eta_an * RTinv)
    #                          - exp(pars.n_an * F * (1 - pars.beta_an) * eta_an * RTinv))
    #i_dl_an = pars.i_ext * pars.A_fac_an - i_Far_an
    #eta_ca = SV[1] - pars.dPhi_eq_ca
    #i_Far_ca = pars.i_o_ca * (exp(-pars.n_ca * F * pars.beta_ca * eta_ca * RTinv)
    #                          - exp(pars.n_ca * F * (1 - pars.beta_ca) * eta_ca * RTinv))
    #i_dl_ca = -pars.i_ext * pars.A_fac_ca - i_Far_ca
    # dSV_dt[i] = i_dl_an * pars.C_dl_an_inv
    # dSV_dt[1] = i_dl_ca*pars.C_dl_ca_inv

    # Testing the other functions
    #D_o = DiffCoeffAnode(SV[ptr.X_an_ptr])
    #print('Outside Fn',D_o)

    # Initialize dSV_dt
    dSV_dt = np.zeros_like(SV)
    #print(dSV_dt)

    # Temperature
    #   - Anode
    #print('ptr.T_an_ptr',ptr.T_an_ptr)
    for i in ptr.T_an_ptr:
        #print(i)
        dSV_dt[i] = 0
    #   - Seperator
    #print('ptr.T_sep_ptr', ptr.T_sep_ptr)
    for i in ptr.T_sep_ptr:
        #print(i)
        dSV_dt[i] = 0
    #   - Cathode
    for i in ptr.T_ca_ptr:
        dSV_dt[i] = 0
    #print(dSV_dt)

    # Voltage (May not use this SV)
    # Currently combining all of the voltage SV
    Voltage_ptr = np.concatenate((ptr.phi_an_ptr, ptr.phi_elyte_an_ptr, ptr.phi_sep_ptr, ptr.phi_ca_ptr, ptr.phi_elyte_ca_ptr))
    for i in Voltage_ptr:
        dSV_dt[i] = 0
    #print(dSV_dt)

    # Concentration
    #   - Anode
    #print('ptr.X_an_ptr', ptr.X_an_ptr)
    for i in ptr.X_an_ptr:
        #print(i)
        # Outer node
        dSV_dt[i] = 0  # %%%%%%% Put in actual values
        if pars.N_part > 1:
            if pars.N_part == 2:
                # Treat as center node
                dSV_dt[i+1] = 0  # %%%%%%% Put in actual values
            else:
                # Middle nodes
                for j in range(pars.N_part - 2):
                    dSV_dt[i + (j+1)] = 0  # %%%%%%% Put in actual values
                # Ceneter Node
                dSV_dt[i+(pars.N_part -1 )] = 0  # %%%%%%% Put in actual values
    #print('Before Sep',dSV_dt)
    #   - Separator
    for i in ptr.X_sep_ptr:
        dSV_dt[i] = 0  # %%%%%%% Put in actual values
    #   - Cathode
    for i in ptr.X_ca_ptr:
        # Outer node
        dSV_dt[i] = 0  # %%%%%%% Put in actual values
        if pars.N_part > 1:
            if pars.N_part == 2:
                # Treat as center node
                dSV_dt[i + 1] = 0  # %%%%%%% Put in actual values
            else:
                # Middle nodes
                for j in range(pars.N_part - 2):
                    dSV_dt[i + (j+1)] = 0  # %%%%%%% Put in actual values
                # Ceneter Node
                dSV_dt[i + (pars.N_part - 1)] = 0  # %%%%%%% Put in actual values

    #print(dSV_dt)
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


def VoltEquibAnode(X):
    # Returns the equilibrium voltage value that corresponses to lithiation
    # Currently Graphite, using Li ref electrode
    phi = 2
    return phi


def VoltEquibCathode(X):
    # Returns the equilibrium voltage value that corresponses to lithiation
    # Currently LiFePO_4, using Li ref electrode
    phi = 2
    return phi


def VoltElyte(X, T):
    # Returns the voltage of the electrolyte based on the Li^+ concentration and temperature
    # Currently salt is LiPF_6
    phi = 2.5
    return phi
