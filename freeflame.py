    
import cantera as ct
import os
import csv
import matplotlib.pyplot as plt
import numpy as np


loglevel = 1  # amount of diagnostic output (0 to 5)

ratio = 2.00  # additional points will be added if the ratio of the spacing on either side of a grid point exceeds this value
slope = 0.05 # maximum difference in value between two adjacent points, scaled by the maximum difference in the profile (0.0 < slope < 1.0). Adds points in regions of high slope.
curve = 0.05  # maximum difference in slope between two adjacent intervals, scaled by the maximum difference in the profile (0.0 < curve < 1.0). Adds points in regions of high curvature.
prune = 0.03 # if the slope or curve criteria are satisfied to the level of ‘prune’, the grid point is assumed not to be needed and is removed. Set prune significantly smaller than ‘slope’ and ‘curve’. Set to zero to disable pruning the grid.



M_O2= 32/1000   
M_N2= 28.02/1000*3.76
    

def calc(name,tburner,width,pressure,sim_type,formula,equivalence_ratio,composition,mole_frac,stocking_values,resolution):
    gas = ct.Solution(formula)   #choosing the chemistry formulas of gri30
    gas.TP = tburner,pressure            #setting the gaz temperature and pressure
    dict_composant={}
    dict_M_num=[None] * len(composition)  
    dict_M_denum=[None] * len(composition) 
    resolution=resolution
    
                              #   STARTING THE BIG BIG LOOP
    n=1    
         

   
    for i in range(len(composition)):
        
        # composition_string={composition[composant]:mole_frac[composant]}
        dict_composant[composition[i]]= mole_frac[i]
        if composition[i] == "H2":
            dict_M_num[i]=2.016/1000*mole_frac[i]/equivalence_ratio
            dict_M_denum[i]=mole_frac[i]/equivalence_ratio
        if composition[i] == "CH4":
            dict_M_num[i]=16/1000*mole_frac[i]/equivalence_ratio
            dict_M_denum[i]=mole_frac[i]/equivalence_ratio
        if composition[i] == "NH3":
            dict_M_num[i]=17.031/1000*mole_frac[i]/equivalence_ratio
            dict_M_denum[i]=mole_frac[i]/equivalence_ratio
 
    gas.set_equivalence_ratio(equivalence_ratio,dict_composant, {'O2':1.0, 'N2':3.76})  #setting the molecular equation
    M_tot=(sum(dict_M_num)+M_O2+M_N2)/(sum(dict_M_denum)+1+3.76)

   
                        # CREATING THE STAGNATION FLOW OBJECT
    
    f = ct.FreeFlame(gas, width=width)
    
    f.set_grid_min(resolution) 
    f.set_refine_criteria(ratio=ratio, slope=slope, curve=curve)
       
    f.show_solution()
    
            # Solve with mixture-averaged transport model                      
                        #SOLVING THE SIMULATION
    try:
        f.solve(loglevel=loglevel, auto=True)
    except ct.CanteraError as e:
        print('Error: Did not make a freeflame =', n, e)
    # sim.solve(loglevel,refine_grid=True, auto=True)
    

                        # set the surface state
    f.show_solution()

                    # Solve with multi-component transport properties
    f.transport_model = 'Multi'
    f.solve(loglevel) # don't use 'auto' on subsequent solves
    f.show_solution()
    
    
    print(" \n\n Tmax is %.2f K , maxHeat is %.2f W/m^3  with the conditions: p=%s pa, tburner=%s K ,width= %s m , eq= %s" % (max(f.T),max(f.heat_release_rate),pressure,tburner,width,equivalence_ratio))
    
  
    #                              ## FINDING VELOCITY FOR STRAIN RATE
    rho_f = pressure / (8.314 / M_tot * tburner)     # fuel inlet density                        # fuel inlet mass flow rate (kg/m^2/s)
    # vel_f = mdotf / rho_f                  # fuel inlet velocity
     #compute flame thickness
    z= f.flame.grid
    T = f.T
    size = np.size(z)-1
    grad = np.zeros(size)
    for i in range(size):
      grad[i] = (T[i+1]-T[i])/(z[i+1]-z[i])
    thickness = (max(T) -min(T)) / max(grad)
    t_inv=f.velocity[0]/thickness
    vel=t_inv*width
    mdot_calculated= rho_f*vel

    if stocking_values=='y':
        outfile = '%s.h5' %name
        if os.path.exists(outfile):
            os.remove(outfile)
        try:
             f.write_hdf('%s.h5' %name, group='laminar flame',mode='a', 
                             description='solution with mixture-averaged transport')
        except ImportError:
              f.save('adiabatic_flame.xml', 'mix',
                     'solution with mixture-averaged transport')      

    return mdot_calculated,t_inv
    


   
