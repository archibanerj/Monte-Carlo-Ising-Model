import numpy as np
import energy
# Script for calculating thermodynamic quantities

def magnetisation(lattice):
    net = np.sum(lattice == 1) - np.sum(lattice == -1)
    m_perspin = net/(np.size(lattice, axis =0)**2)
    return m_perspin