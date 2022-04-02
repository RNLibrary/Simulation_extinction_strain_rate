# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 15:04:09 2022

@author: Rani Naaman
"""
""" Main code for Rani's 2022 project """
import json
import sys
import freeflame
import ImpingingJet
import time
import concurrent.futures
from tqdm import tqdm
import multiprocessing  
import psutil
from concurrent.futures import ProcessPoolExecutor


def validate_key_config(key, config_dict):
    if key not in config_dict:
        print(f"Error the {key} key is required in {config_dict}")
        sys.exit()


with open("./data.json", 'r', encoding='utf-8') as f:
    config = json.load(f)

validate_key_config('simulations', config)

simulations = config['simulations']


def run_simulation(simulation, delay):
    validate_key_config('name', simulation)
    validate_key_config("tburner", simulation)
    validate_key_config("tsurf", simulation)
    validate_key_config("width", simulation)
    validate_key_config('pressure', simulation)
    validate_key_config('type', simulation)
    validate_key_config('formula', simulation)
    validate_key_config('composition', simulation)
    validate_key_config('equivalence ratio', simulation)
    validate_key_config('mole_frac', simulation)
    validate_key_config('stocking_values', simulation)

    name = simulation['name']    
    tsurf = simulation['tsurf']
    tburner = simulation['tburner']
    width = simulation['width']
    stocking_values=simulation['stocking_values']
    pressure = simulation['pressure']
    sim_type = simulation['type']
    formula = simulation['formula']
    equivalence_ratio = simulation['equivalence ratio']
    composition = simulation['composition']
    mole_frac = simulation['mole_frac']
    
    freeflame_values=freeflame.calc(name,tburner,width,pressure,sim_type,formula,equivalence_ratio,composition,mole_frac,stocking_values)
    # print(f"here is {freeflame_values[1]}")
    mdot_calculated=freeflame_values[0]
    strain=freeflame_values[1]
    
    print(f"Running {simulation['name']}")
    # if str.lower(sim_type) == "freeflame":
    #     freeflame.calc()
    if str.lower(sim_type) == "impingingjet":
        ImpingingJet.calc(name,tburner,tsurf,width,mdot_calculated,pressure,sim_type,formula,equivalence_ratio,composition,mole_frac,strain,stocking_values)       

    time.sleep(2 + delay)

                        # RUNNING FREEFLAME

n_cores=psutil.cpu_count(logical=False)
# def main():
#     with ProcessPoolExecutor(max_workers=n_cores) as executor:
    
#             for simulation in simulations:
    
#                 res = executor.map(run_simulation, simulation)  
#     print(list(res))

def main():
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_cores) as executor:
      futures = []
      delay = 1
      with tqdm(total=len(simulations)) as pbar:
          for simulation in simulations:
              futures.append(executor.submit(run_simulation, simulation, delay))
              delay += 2
          for future in concurrent.futures.as_completed(futures):
              pbar.update(1)
              print(future.result())   
      
if __name__ == '__main__':
    main()    

# with concurrent.futures.ThreadPoolExecutor(max_workers=n_cores) as executor:
#     futures = []
#     delay = 1
#     with tqdm(total=len(simulations)) as pbar:
#         for simulation in simulations:
#             futures.append(executor.submit(run_simulation, simulation, delay))
#             delay += 2
#         for future in concurrent.futures.as_completed(futures):
#             pbar.update(1)
#             print(future.result())


    

