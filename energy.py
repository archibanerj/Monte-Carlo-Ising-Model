import numpy as np

# In this script we calculate the energy of a given spin configuration of size
# N x N. We consider a simple square lattice in 2 dimensions
# This is done using nearest neighbour couplings between the spins
# We also consider a homogenous spin system with uniform J for all couplings

'''More efficient energy management if we use classes?'''

def p(i,N):
    # To implement periodic boundary conditions
    # This function send the index back to the 
    # other edge of the spin lattice if it reaches an edge while adding the energies.
    pos = i
    if i==-1:
        pos = i+N
    elif i==N:
        pos = i-N
    return pos
    

def hamilton(spins, J=1):
    N = np.size(spins,axis = 0) # Number of spins N
    E = 0
    for i in range(N): 
        for j in range(N):
            # We implement periodic boundary conditions here
            # using modulo
            E +=  spins[i][j]*spins[(i+1)%N][j] + spins[i][j]*spins[(i-1)%N][j]  +  spins[i][j]*spins[i][(j+1)%N]  +  spins[i][j]*spins[i][(j-1)%N]  
    
    return -J*E/2 # 2 for overcounting - each pair is counted twice.


def energy_update(NN, J, N):
    '''
    Function to calculate the energy difference between given 
    configuration and the resulting configuration with (i,j)-th spin flipped

    NN is the combined set of nearest neighbours of NN[1:] and the middle spin NN[0]
    N is the number of lattice points along an edge
    J is the coupling energy
    '''

    # Calculating energy contribution of the (i,j)th spin via nearest neighbour couplings
    del_E = -J * NN[0] * np.sum(NN[1:])

    # We flip the (i,j) th spin to get -del_E as the new energy contribution. 
    # Thus, the difference in energy between the two states (flipped - original) is
    diff =  -2*del_E

    return diff
