import numpy as np
import matplotlib.pyplot as plt
import energy
import montecarlo
import thermodynamics as td
import os
import time

# We are setting k=J=1 for the simulations #

''' How to make this faster? for T~10^6, N~10^(1-2) and Temp ~ 10^1 
-> There are 10^{8-9} MC cycles that are completed
One way is to generate data catalogued into various sizes and temperatures. 
then use the stored data to calculate Thermodynamics properties '''

## MONTE CARLO PARAMETERS ##
t_eq = 10**3 # equilibriation time
t_a = 10**2 # Autocorrelation time (time between collection of samples)
N = [30] #[10,20,30, # Some chosen values for number of spins along a row/column
n = 100 # number of samples collected
T = t_eq + n*t_a # Total monte carlo time

dt = 0.1
Temp = np.arange(start=0.1,stop = 5.1, step = dt, dtype=float) # temperatures


store_disk = True # For storing samples to disk
#snap = True # If set to False, script does not collect snapshots and stop at equilibrium time
#snapshots=[] # for collecting snaps

# thermodynamic quantities
# spon_mag_per_spin = np.zeros_like(Temp) # spontaneous magnetisation per spin
# E = np.zeros_like(Temp) # Energy
# net_energy_change = 0

if store_disk:
    parent_dir = "C:\\Users\\archi\\Documents\\Statistical Physics\\Monte Carlo\\"
    folder = "snapshots\\"
    current_dir = os.path.join(parent_dir,folder)
    #os.mkdir(current_dir)
    for size in N:

        print("Storing for N = ",size)

        folder =  "N=" + str(size) + "\\"
        temp_dir = os.path.join(current_dir, folder)
        os.mkdir(temp_dir)

        for temp_index in range(len(Temp)):
            # note that for the next temperature, the starting spin lattice will be the last one obtained
            # in the simulation for same size and just lower temperature. Will this cause a problem?

            lattice = 2*np.random.randint(0,2,(size,size))-1

            folder = "temp=" + str(round(Temp[temp_index],1))  + "\\"
            temp_dir2 = os.path.join(temp_dir,folder)
            os.mkdir(temp_dir2)

            for t in range(T): 
                start = time.time()
                lattice,energy_change = montecarlo.mcCycle(lattice, size, Temp[temp_index]) # 1 MC cycle
                end = time.time()
                if t==0 and temp_index==0:
                    print("Time taken for one MC Cycle: ", start-end)
                    print('This will be repeated T = ',T)
                #if t==t_eq: 
                    #print('Equilibrium has been reached')
                #    if not snap:
                #        break
                if (t-t_eq)%t_a == 0:
                    r = (t-t_eq)/t_a
                    '''
                    #if r%1000 == 0: print("1000 samples collected")
                    spon_mag_per_spin[temp_index] += td.magnetisation(lattice)
                    #E[temp_index] += E[temp_index] + net_energy_change # Is this update operation correct?
                    net_energy_change = 0 #refreshing temp variable for next set of collections
                    #snapshots.append(lattice) # Collecting the lattice snapshot
                    '''
                    if r>0:
                        path =  temp_dir2 + "n=" + str(r)+".csv"
                        start = time.time()
                        np.savetxt(path, lattice, delimiter=' ')
                        end = time.time()
                        if r==0: 
                            print('time for saving - ',(start-end))
            '''
            spon_mag_per_spin[temp_index] /= n# Obtaining ensemble average (/n) 
            spon_mag_per_spin[temp_index] = np.abs(spon_mag_per_spin[temp_index])
            #E[temp_index] /= n
            '''
            print("Done for T = ", Temp[temp_index])

# Plotting spontaneous magnetisation
# plt.plot(Temp,spon_mag_per_spin)
# plt.show()