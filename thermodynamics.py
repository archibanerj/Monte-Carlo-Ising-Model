import numpy as np
import energy

'''
Script for calculating thermodynamic quantities -
Also add functions for calculating Entropy, partition function, etc etc.
Also, calculation of critical exponents.
'''


def magnetisation(lattice):
    net = np.sum(lattice == 1) - np.sum(lattice == -1)
    m_perspin = net/(np.size(lattice, axis =0)**2)
    return m_perspin

def mag_susc(mags,n,T):
    '''
    magnetic susceptibility
    mags - set of magnetisations of equilibrium configurations collected at a temperature t
    n - number of spins!
    t - temperature
    '''
    term1 = np.sum(mags**2)/np.size(mags)
    term2 = (np.sum(np.abs(mags))/np.size(mags))**2
    return (n/T)*(term1-term2)

def spec_heat_stat(energies, T):
    '''
    Takes in energies of different equilibrium configurations at temperature T 
    Calculates specific heat statistically - like the magsusc function
    '''
    term1 = np.sum(energies**2)/np.size(energies)
    term2 = (np.sum(np.abs(energies))/np.size(energies))**2
    return (1/(T**2))*(term1-term2)

def spec_heat_diff(energies,temp):
    '''
    Takes in ensemble averaged value of system energies at each temperature

    Calculation of specific heat as a function of temperature
    Also returns critical temperature
    '''
    C = np.array(energies)[1:len(energies)] - np.array(energies)[0:(len(energies)-1)]
    dt = temp[1]-temp[0]
    C /= dt
    temp = temp[1:len(energies)]

    # Another method of calculating 

    return C,temp,(np.argmax(C)*dt)