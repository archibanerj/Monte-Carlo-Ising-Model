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

def Wolff(lattice,beta): 
    '''
    If we define this to be a part of a class, then we will not have to take the lattice as an input. 
    This may speed up the code
    '''

    size = np.size(lattice, axis=0)
    choose = np.random.randint(0,size),np.random.randint(0,size) # choosing random site

    stack = [choose]
    cluster = []
    cluster_spin = lattice[choose]

    while len(stack)>0 and len(stack)<size**2:

        #print("Stack = ", stack)
        #print("Cluster= ",cluster)
        
        i = stack[0][0]
        j = stack[0][1]
        if lattice[stack[0]] == cluster_spin:

            # Now we visit the nearest neighbours and add them to cluster with the prob. 1-exp(-2beta) if they are parallel
            nn = [((i+1)%size,j),(i,(j+1)%size),((i-1)%size,j),(i,(j-1)%size)]

            for pos in nn:
                if lattice[pos] == cluster_spin: 
                    pos_already_included = pos in cluster
                    if not pos_already_included:    
                        z = np.random.random(size= 1)
                        if z < 1 - np.exp(-2*beta):
                            stack.append(pos)

            # Now, I flip the initial spin I started from:
            lattice[stack[0]] = - lattice[stack[0]]

        # And I remove it from the array
        cluster.append(stack.pop(0))
        #print(len(stack))

    # Return the updated lattice, and the cluster used in the iteration.
    return lattice, cluster 
