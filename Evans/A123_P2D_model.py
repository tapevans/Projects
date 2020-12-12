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
from scipy.integrate import solve_ivp #integration function for ODE system.
import numpy as np
from A123_P2D_function import residual
from A123_P2D_init import SV_0, t_final, pars, ptr


time_span = np.array([0, t_final])
#print('Shape of SV_0 =', np.shape(SV_0))
solution = solve_ivp(lambda t, y: residual(t, y, pars, ptr), time_span, SV_0, rtol=1e-4, atol=1e-6)
#print('Shape of solution =', np.shape(solution))
#print('Shape of solution.t =', np.shape(solution.t))
#print('Shape of solution.y =', np.shape(solution.y))
#print('Solution.t =', solution.t)
#print('Solution.y =', solution.y)

# %%%%%%%%%%%% Add plotting stuff to another file
from matplotlib import pyplot as plt
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
plt.show()
