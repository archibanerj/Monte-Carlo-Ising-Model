import numpy as np
import matplotlib.pyplot as plt
import energy
import montecarlo
import thermodynamics as td
import os
import time

# We are setting k=J=1 for the simulations #
''' 
2. write functions to perform  (b) Finite Size Scaling (d) the data generation step. Essentially make everything modular. 
3. Use more temperature points near the transition temperature.
4. Make this more efficient. This is terrible. 

'''

def init_func(eq_time, AC_time, no_of_samples, temperatures):
    global t_eq,t_a,N,n,Temp,mag,susc,Energies,spec_heat,below_Tc,above_Tc,mag_above,mag_below,E_above,E_below,parent_dir,folder

    ## MONTE CARLO PARAMETERS ##
    t_eq = eq_time #10**3 # equilibriation time
    t_a = AC_time #10**2 # Autocorrelation time (time between collection of samples)
    N = [] # redundant
    n = no_of_samples # number of samples collected

    '''
    Temp = np.array(range(1,51))/10 # temperatures in steps of 0.1
    Temp1 = np.array(range(200,240))/100 # more dense sampling of temperature around critical value
    Temp = np.unique(np.append(Temp,Temp1)) # final temperature array
    '''

    Temp = temperatures

    ## Thermodynamic Quantities ##
    mag = [] # Ensemble-Averaged Magnetisation at each temperature
    susc = [] # (ensemble averaged)magnetic susceptibility at each temperature
    Energies = [] # Ensemble-Averaged Energy at each temperature
    spec_heat = [] # ensemble averaged specific heat at each temperature. 

    below_Tc = 2.0 # Temperature below Tc for studying magnetisation and energy distributions
    above_Tc = 2.5 # Temperature above Tc for studying magnetisation and energy distributions

    mag_below = [] # Set of all magnetisations of every monte carlo sample collected at given temperature below Tc
    E_below = [] # Set of all Energues of every monte carlo sample collected at given temperature above Tc
    mag_above = [] # Set of all magnetisations of every monte carlo sample collected at given temperature below Tc
    E_above = [] # Set of all Energies of every monte carlo sample collected at given temperature above Tc

    parent_dir = "C:\\Users\\archi\\Documents\\Statistical Physics\\"  # Parent Directory for storing snaps
    folder = "snapshots\\" # folder for storing snaps

def store(size,eq_time,AC_time,no_of_samples):
    ''' 
    Suggested change - change the autocorrelation time to temperature dependant.
    Also, change the following - folder to which everything is saved. This is for the new temperature values
    '''

    global T

    init_func(eq_time,AC_time,no_of_samples)

    T = t_eq + n*t_a # Total monte carlo time
    current_dir = os.path.join(parent_dir,folder)
    print("Storing for N = ",size)

    folder1 =  "N=" + str(size) + "\\"
    temp_dir = os.path.join(current_dir, folder1)
    os.mkdir(temp_dir)

    for temp_index in range(len(Temp)):

        lattice = 2*np.random.randint(0,2,(size,size))-1  # Creating random lattice

        folder2 = "temp=" + str(round(Temp[temp_index],2))  + "\\"
        temp_dir2 = os.path.join(temp_dir,folder2)
        os.mkdir(temp_dir2)

        t_auto = t_a[temp_index] # Autocorrelation time for given temperature
        t = 0 # Monte carlo chain time for the loop below

        #for t in range(T[temp_index]): 
        while t < T[temp_index]:

            start = time.time()
            #lattice,energy_change = montecarlo.mcCycle(lattice, size, Temp[temp_index]) # 1 MC cycle
            lattice,cluster = montecarlo.Wolff(lattice, 1/Temp[temp_index])
            end = time.time()

            '''
            if t==0 and temp_index==0:
                print("Time taken for one MC Cycle: ", end-start)
                print('This will be repeated T = ',T)
            
            '''

            if (t-t_eq)%t_auto == 0:
                r = (t-t_eq)/t_auto
                if r>0:
                    path =  temp_dir2 + "n=" + str(r)+".csv"
                    start = time.time()
                    np.savetxt(path, lattice, delimiter=' ') # Faster alternative?
                    end = time.time()
                    #if r==1.0: 
                        #print('time for saving = ',(end-start))
                        #print('This will be repeated = ',n)
            t += 1
        
        print("Done for T = ", Temp[temp_index])

def plot_stuff(mag,Energies,susc,spec_heat,temperature):

    plt.scatter(temperature, np.array(mag)) 
    plt.title("Magnetisations")
    plt.show()

    plt.plot(temperature, np.abs(np.array(mag)))
    plt.title("Magnetisation - Absolute")
    plt.show()

    plt.plot(temperature, np.array(Energies))
    plt.title("Internal Energy")
    plt.show()

    plt.plot(temperature,np.array(susc))
    plt.title("Mag. Susceptibility")
    plt.show()

    plt.plot(temperature,spec_heat,label='stat')
    plt.legend()
    plt.title('Specific Heat')
    plt.show()

def averaging(size,no_of_samples):
    ''' 
    Calculate all thermodynamic quantities for a given lattice size

    (We shall add more thermodynamic averages here as we add functions for more such 
    TD potentials in the td script)
    '''
    init_func(0,0,no_of_samples) # 0 and 0 are dummy variables 
    # because these are not required for calculation of averages
    fn = os.path.join(os.path.join(parent_dir,folder),"N=" + str(size) + "\\")
    
    for t in Temp:
        path = fn + "temp=" + str(round(t,2)) + "\\"
        m=0
        E=0
        mags_list =[]
        E_list = []
        for i in range(1,n):
            fname = path + "n="+str(round(float(i),1))+".csv"
            lattice = np.loadtxt(fname)
            lattice_mag = td.magnetisation(lattice)
            lattice_E = energy.hamilton(lattice, J=1)
            m += np.abs(lattice_mag) 
            # We are forced to do np.abs here otherwise the magnetisations +1 and -1 mix.
            # To avoid this issue, choose a better value of autocorrelation time. 
            # That would be a more scientific way to treat this problem
            #m += lattice_mag
            E += lattice_E
            mags_list.append(lattice_mag)
            E_list.append(lattice_E)
            '''
            if round(t,1) == below_Tc :
                mag_below.append(lattice_mag)
                E_below.append(lattice_E)
            elif round(t,1) == above_Tc :
                mag_above.append(lattice_mag)
                E_above.append(lattice_E)
            '''

        m /= n
        E /= n
        
        susc.append(td.mag_susc(np.array(mags_list),n=(size**2),T=t))
        spec_heat.append(td.spec_heat_stat(np.array(E_list),T=t))
        mag.append(m)
        Energies.append(E)
        print("Done for t=",t)
    
    return mag,Energies,susc,spec_heat



