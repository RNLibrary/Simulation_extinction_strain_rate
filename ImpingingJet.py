"""
A detached flat flame stabilized at a stagnation point

This script simulates a mix of hydrogen and methane with a variable array of the dilution
ratio (H2) and of the flame speed (mdot).

Other info about the simulation:
The flame is stabilized in a strained flowfield at an axisymmetric stagnation point
on a non-reacting surface. The solution begins with a flame attached to the inlet 
(burner), and the mass flow rate is progressively increased, causing the flame to 
detach and move closer to the surface.
"""


import cantera as ct
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import h5py

loglevel = 1  # amount of diagnostic output (0 to 5)


ratio = 2.00  # additional points will be added if the ratio of the spacing on either side of a grid point exceeds this value
slope = 0.05 # maximum difference in value between two adjacent points, scaled by the maximum difference in the profile (0.0 < slope < 1.0). Adds points in regions of high slope.
curve = 0.05  # maximum difference in slope between two adjacent intervals, scaled by the maximum difference in the profile (0.0 < curve < 1.0). Adds points in regions of high curvature.
prune = 0.03 # if the slope or curve criteria are satisfied to the level of ‘prune’, the grid point is assumed not to be needed and is removed. Set prune significantly smaller than ‘slope’ and ‘curve’. Set to zero to disable pruning the grid.
resolution=1e-4


M_O2= 32/1000   
M_N2= 28.02/1000*3.76
    

def calc(name,tburner,tsurf,width,mdot,pressure,sim_type,formula,equivalence_ratio,composition,mole_frac,strain):
    gas = ct.Solution(formula)   #choosing the chemistry formulas of gri30
    gas.TP = tburner,pressure            #setting the gaz temperature and pressure
    dict_composant={}
    dict_M_num=[None] * len(composition)  
    dict_M_denum=[None] * len(composition) 
    matrice_grid=[]
    matrice_hrr=[]
    matrice_T=[]
    raffinement=10
    mdot_variable=np.linspace(0,mdot,raffinement)
    test_raf=0
    test_conv=0
    double_converge=0
                              #   STARTING THE BIG BIG LOOP
    n=1    
    y=1
         
                                # WRITING THE BELOW INFORMATION IN VALUES2.XML
    row1=['T','maxHeat','strainRate','','p','tburner','tsurf','mdot','witdh','Resolution','eq']
    with open(r'VALUES2_%s.csv' %name,"a",newline="") as csvfile:
                csvwriter= csv.writer(csvfile)
                csvwriter.writerow(row1)
           
    while True:
        
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
        
        sim = ct.ImpingingJet(gas=gas, width=width)
        
                              # set the mass flow rate at the inlet
        sim.inlet.mdot = mdot_variable[y]
        
                            # set the surface state
        sim.surface.T = tsurf
        
                             #SETTING THE RESOLUTIONS   
        
        sim.set_grid_min(resolution)       
        sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)      
        
        sim.set_initial_guess(products='equil')  # assume adiabatic equilibrium products
        sim.show_solution()
        
                       
                            #SOLVING THE SIMULATION
        try:
            sim.solve(loglevel,refine_grid=True, auto=True)
        except ct.CanteraError as e:
            print('Error: Did not converge at attempt n =', n, e)
         
                
                #  SAVING MAX T AND MAX HRR     
        maxHeat=max(sim.heat_release_rate)
        maxT=max(sim.T)
        
        
                                     #NORMALIZING FOR FLAME STRUCTURE
        normalized_grid = (sim.grid - min(sim.grid)) / (max(sim.grid) - min(sim.grid))
        normalized_hrr = (sim.heat_release_rate - min(sim.heat_release_rate)) / (max(sim.heat_release_rate) - min(sim.heat_release_rate))
        normalized_T = (sim.T - min(sim.T)) / (max(sim.T) - min(sim.T))
        matrice_T.append(normalized_T)   
        matrice_hrr.append(normalized_hrr)  
        matrice_grid.append(normalized_grid)  
        
                                     ## FINDING VELOCITY FOR STRAIN RATE
        rho_f = pressure / (8.314 / M_tot * tburner)     # fuel inlet density    
        mdotf= mdot_variable[y]                     # fuel inlet mass flow rate (kg/m^2/s)
        vel_f = mdotf / rho_f                  # fuel inlet velocity
   
    
    
    
                                    # # WRITING THE BELOW INFORMATION IN VALUES2.XML
        row=[str(maxT),str(maxHeat),str(vel_f/width),"",str(pressure),str(tburner),str(tsurf),str( mdot_variable[y]),str(width),str(resolution),str(equivalence_ratio)]
        print('results of',name)
        print(row)
        with open(r'VALUES2_%s.csv' %name,"a",newline="") as csvfile:
                csvwriter= csv.writer(csvfile)
                csvwriter.writerow(row)  
        
        if n==1:
            # f = h5py.File('%s.h5' %name)
            # f['mydataset'][:] = 0
            try:
                sim.write_hdf('%s.h5' %name, group='step #%s'% n,mode='a', 
                                description='solution with mixture-averaged transport')
            except ImportError:
                 sim.save('adiabatic_flame.xml', 'mix',
                        'solution with mixture-averaged transport')  
        
        else:
                
            try:
                # save to HDF container file if h5py is installed
                sim.write_hdf('%s.h5' %name, group='step #%s'% n,
                                    description='solution with mixture-averaged transport')
            except ImportError:
                sim.save('adiabatic_flame.xml', 'mix',
                       'solution with mixture-averaged transport')  
       
                    
                            # CONDITIONSSS
        if not np.isclose(max(sim.T), tsurf):
            n +=1
            y +=1
            if y==10 :
                print('index exceeded expectations')
                new_mdot=mdot_variable[y-1]*2
                old_mdot=mdot_variable[y-1]
                mdot_variable=np.linspace(mdot,new_mdot,raffinement)
                y=1 
            print ('la nouvelle vitesse est de',  mdot_variable[y])
            sim.set_initial_guess()
            
        elif np.isclose(max(sim.T), tsurf) and test_conv==0 and n!=1:
             print("la flamme a convergé pour la première fois")
             n +=1
             test_conv=1
             old_mdot=mdot_variable[y-1]
             max_mdot=mdot_variable[y]

             mdot_variable=np.linspace(old_mdot,max_mdot,raffinement)
             y=1
             print ('la nouvelle vitesse est de',  mdot_variable[y])
             sim.set_initial_guess()
    
        elif np.isclose(max(sim.T), tsurf) and test_conv==0 and n==1:
            print("la flamme s'est éteinte à la première simulation")

            max_mdot=mdot_variable[y]

            mdot_variable=np.linspace(0,max_mdot,raffinement)
            y=1
            n+=1
            test_conv=1
            print ('la nouvelle vitesse est de',  mdot_variable[y])
            sim.set_initial_guess()
            
        elif np.isclose(max(sim.T), tsurf) and test_conv==1: 
            if  np.isclose(matrice_T[n-1][-1],matrice_T[n-2][-1]) :
                print('La simulation a convergé 2 fois de suite')
                double_converge=1
            print ('la vitesse est  de',  mdot_variable[y])
            print ('convergence trouvé')
            break
    
            
            

                           #CREATING FLAME STRUCTURE
    fig=plt.figure()
    plt.rc('font', family='STIXGeneral',size=22)
    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    plt.figure(figsize=(6, 5))
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'STIXGeneral'
    plt.rcParams['mathtext.it'] = 'STIXGeneral:italic'
    plt.rcParams['mathtext.bf'] = 'STIXGeneral:bold'
    if double_converge==1:
        plt.plot(matrice_grid[n-3], matrice_hrr[n-3],'-',c='k',label='HRR',linewidth=2.5)
        plt.plot(matrice_grid[n-3],matrice_T[n-3],'-.',c='r',label='TEMPERATURE',linewidth=1)
    else:
        plt.plot(matrice_grid[n-2], matrice_hrr[n-2],'-',c='k',label='HRR',linewidth=2.5)
        plt.plot(matrice_grid[n-2],matrice_T[n-2],'-.',c='r',label='TEMPERATURE',linewidth=1)
    # plt.plot(normalized_grid, normalized_O2,'-',c='g',label='[O2]',linewidth=1)
    plt.ylabel(r'NORMALIZED DATA')
    plt.xlabel(r'NORMALIZED DISTANCE (m)')
    plt.title('%s'  %name)
    # plt.title('AMMONIA AND HYDROGEN LAMINAR FLAME #%s \n WITH [H2] AT %s (%%) AND MDOT AT %s (kg/m^2/s)'% (i,hh,md),size=19)
    plt.legend(bbox_to_anchor=(0, 0.2), loc='lower left',prop={'size': 12})
    plt.xlim(0.0,1.0)
    plt.savefig(r' %s.png' %name, bbox_inches='tight', dpi=150)
    # plt.savefig(r' HH[%s].png' %i, bbox_inches='tight', dpi=150)
    # plt.show()
    mdotf= mdot_variable[y-1]                     # fuel inlet mass flow rate (kg/m^2/s)
    vel_f = mdotf / rho_f                  # fuel inlet velocity
    strain=vel_f/width
    print('this is run',name)
    print("Strain rate of laminar : " ,strain)
    print('thats the strain rate of free flame',mdot/rho_f/width)
    
        
    return fig
        
    
