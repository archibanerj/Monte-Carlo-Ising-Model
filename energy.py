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
    

def hamilton(spins, J):
    N = np.size(spins,axis = 0) # Number of spins N
    E = 0
    for i in range(N): 
        for j in range(N):
            # We implement periodic boundary conditions here
            # using the function p
            # Problem here - overcounting. Fix this. 
            E +=  spins[i][j]*spins[p(i+1,N)][j] + spins[i][j]*spins[p(i-1,N)][j]  +  spins[i][j]*spins[i][p(j+1,N)]  +  spins[i][j]*spins[i][p(j-1,N)]  
    
    return -J*E/2 # Does the /2 fix the overcounting?


def energy_update(spins, J, flip):
    '''
    Function to calculate the energy difference between given 
    configuration and the resulting configuration with (i,j)-th spin flipped

    Can be simplified - we dont need to know the entire spins data - just the 4 n.n.s would do
    '''
    N = np.size(spins,axis=0)
    i,j = flip

    # Calculating energy contribution of the (i,j)th spin via nearest neighbour couplings
    del_E = -J * (spins[i][j]*spins[p(i+1,N)][j] + spins[i][j]*spins[p(i-1,N)][j]  +  spins[i][j]*spins[i][p(j+1,N)]  +  spins[i][j]*spins[i][p(j-1,N)])

    # We flip the (i,j) th spin to get -del_E as the new energy contribution. 
    # Thus, the difference in energy between the two states (flipped - original) is
    diff =  -2*del_E

    return diff
