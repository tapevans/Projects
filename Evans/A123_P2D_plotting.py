"""
    This file runs and executes a youyr model, calculating the cell/device/system properties of interest as a function of time. 

    The code structure at this level is meant to be very simple, and delegates most functions to lower-level modules, which can be updated with new capabilties, over time.  The code:

        1 - Reads inputs and initializes the model
        
        2 - Calls the residual function and integrates over the user-defined    
            time span.
        3 - The simulation then returns the solution vector at the end of the 
            integration, which can be processed as needed to plot or analyze any quantities of interest.
"""

# Import necessary modules:
import numpy as np
from A123_P2D_function import VoltEquibAnode, VoltEquibCathode
from A123_P2D_init import t_final, pars, ptr
from matplotlib import pyplot as plt

def plots(solution):
    # Convert mole fraction into cell voltage
    AnodeVoltage = solution.y

    # plt.show()

    # If EIS, convert voltage into phase lag, Z_Re, and Z_Im




#anode_temp = np.transpose(solution.y[ptr.T_an_ptr[0]][:])
#print(anode_temp)
#print('Shape of anode_temp =', np.shape(anode_temp))
#for var in anode_temp:
    #print('Shape soln.t', np.shape(solution.t))
    #print('Shape var', np.shape(var))
    #pbreak
#plt.plot(solution.t, anode_temp)

for var in solution.y:
    plt.plot(solution.t, var)

#plt.legend(['Anode double layer', 'Cathode double layer'])

