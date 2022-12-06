import energy
import numpy as np


def mcStep(lattice, flip, Temp):
    # 1 MC step
    i,j = flip
    del_E = energy.energy_update(lattice,J=1,flip=(i,j))
    energy_change = 0 # actual energy change at end of spin flip
    if del_E < 0:
        lattice[i,j] = - lattice[i,j]
        energy_change = del_E
    else :
        z = np.random.random(size= 1)
        if z < np.exp(-del_E/Temp) :
            lattice[i,j] = - lattice[i,j]
            energy_change = del_E

    return lattice, energy_change

def mcCycle(lattice, N, Temp):
    # 1 MC cycle
    energy_change_MC_cycle = 0 # net change in energy in 1 MC cycle
    for i in range(N):
        for j in range(N):
            lattice,del_E = mcStep(lattice, flip = (i,j), Temp = Temp) 
            energy_change_MC_cycle += del_E
    return lattice,energy_change_MC_cycle

