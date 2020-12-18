"""
    This file runs and executes a your model, calculating the cell/device/system properties of interest as a function of time.

    The code structure at this level is meant to be very simple, and delegates most functions to lower-level modules, which can be updated with new capabilties, over time.  The code:

        1 - Reads inputs and initializes the model
        
        2 - Calls the residual function and integrates over the user-defined    
            time span.
        3 - The simulation then returns the solution vector at the end of the 
            integration, which can be processed as needed to plot or analyze any quantities of interest.
"""

# Import necessary modules:
from scipy.integrate import solve_ivp  # integration function for ODE system.
#from assimulo.problem import Implicit_Problem
#from assimulo.solvers import IDA
import numpy as np
#from A123_P2D_function import residual
from A123_P2D_function import dSVdt
from A123_P2D_init import SV_0, t_final, pars, ptr  # dSVdt_0
from matplotlib import pyplot as plt

time_span = np.array([0, t_final])
t0 = time_span[0]
#print('Shape of SV_0 =', np.shape(SV_0))

#model = Implicit_Problem(residual, SV_0, dSVdt_0, t0)
#sim = IDA(model)
#ncp = 500  # Number of communication points (number of return points)
#t, y, yd = sim.simulate(t_final, ncp)

solution = solve_ivp(lambda t, y: dSVdt(t, y, pars, ptr), time_span, SV_0, 'BDF', rtol=1e-4, atol=1e-6)  # BDF  Radau
#print('Shape of solution =', np.shape(solution))
#print('Shape of solution.t =', np.shape(solution.t))
#print('Shape of solution.y =', np.shape(solution.y))
#print('Solution.t =', solution.t)
print('Solution.y =', solution.y)
file = open("solutiont.txt", "w")
file.write(str(solution.t))
file.close()

file = open("solutiony.txt", "w")
file.write(str(solution.y))
file.close()
# %%%%%%%%%%%% Add plotting stuff to another file

anode_temp = np.transpose(solution.y[ptr.T_ptr][:])
delta_phi_an = np.transpose(solution.y[ptr.delta_phi_ptr][:])
Li_ion_an = np.transpose(solution.y[ptr.X_Li_ptr][:])
print(Li_ion_an)
#print(solution.t)
#print('Shape of time=', np.shape(solution.t))
#print(anode_temp)
#print('Shape of anode_temp =', np.shape(anode_temp))

#for var in anode_temp:
    #print('Shape soln.t', np.shape(solution.t))
    #print('Shape var', np.shape(var))
    #pbreak
#plt.plot(solution.t, anode_temp)
#plt.show()

#plt.plot(solution.t, delta_phi_an)
#plt.show()

plt.plot(solution.t, Li_ion_an)
plt.show()

for var in solution.y:
    plt.plot(solution.t, var)
#plt.show()
#plt.legend(['Anode double layer', 'Cathode double layer'])

