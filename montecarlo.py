import energy
import numpy as np


def mcStep(nn,N,Temp):
    # 1 MC step
    del_E = energy.energy_update(nn,J=1,N=N)
    flip = False
    energy_change = 0 # actual energy change at end of spin flip - only updated if spin is flipped, else not.
    if del_E < 0:
        flip = True
        energy_change = del_E 
    else :
        z = np.random.random(size= 1)
        if z < np.exp(-del_E/Temp) :
            flip = True
            energy_change = del_E

    return flip, energy_change

def mcCycle(lattice, N, Temp):
    # 1 MC cycle
    energy_change_MC_cycle = 0 # net change in energy in 1 MC cycle
    for i in range(N):
        for j in range(N):
            nn = np.array([lattice[i,j],lattice[(i+1)%N,j],lattice[i,(j+1)%N],lattice[(i-1)%N,j],lattice[i,(j-1)%N]])
            flip,del_E = mcStep(nn, N, Temp = Temp)  # The mcStep performs 1 mcStep and tells us if we 
            # should flip the i,j spin or not

            if flip:
                lattice[i,j] = - lattice[i,j]

            energy_change_MC_cycle += del_E
    return lattice,energy_change_MC_cycle

