import numpy as np
from math import exp

# Constants
F = 96485
R = 8.3145


def residual(t, SV, pars, ptr):
# %%%%%%%% Fix this
    RTinv = 1/R/pars.T
    dSV_dt = np.zeros_like(SV)
    
    eta_an = SV[0] - pars.dPhi_eq_an
    i_Far_an = pars.i_o_an*(exp(-pars.n_an*F*pars.beta_an*eta_an*RTinv)
                      - exp(pars.n_an*F*(1-pars.beta_an)*eta_an*RTinv))
    i_dl_an = pars.i_ext*pars.A_fac_an - i_Far_an
    dSV_dt[0] = i_dl_an*pars.C_dl_an_inv
    
    eta_ca = SV[1] - pars.dPhi_eq_ca
    i_Far_ca = pars.i_o_ca*(exp(-pars.n_ca*F*pars.beta_ca*eta_ca*RTinv)
                      - exp(pars.n_ca*F*(1-pars.beta_ca)*eta_ca*RTinv))
    i_dl_ca = -pars.i_ext*pars.A_fac_ca - i_Far_ca

    dSV_dt[1] = i_dl_ca*pars.C_dl_ca_inv
    
    return dSV_dt


def DiffCoeffAnode(X):
    # Change to be a function of concentration
    # Currently Graphite
    D_o = 1e-16
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
