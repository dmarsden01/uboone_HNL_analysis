#Script for calculating limit from .json
import os, sys, string, time
import json
import matplotlib.pyplot as plt
import numpy as np
import csv
import pyhf

def Load_json(loc):
    """Load a single .json file that contains pyhf model parameters."""
    full_loc = f"./jsons/{loc}"
    if os.path.isfile(full_loc) != True: 
        print(f"File {loc} does not exist, exiting.")
        return 1
    
    f = open(full_loc) #'f'./jsons/{HNL_mass}_{name_type}_Test.json'

    model = json.load(f)
    f.close()

    return model

def Load_data(loc):
    """Load a file that contains the data histogram and auxdata."""
    full_loc = f"./observed_data/{loc}"
    if os.path.isfile(full_loc) != True: 
        print(f"File {full_loc} does not exist, exiting.")
        return 1
    
    f = open(full_loc) #'f'./jsons/{HNL_mass}_{name_type}_Test.json'

    data = json.load(f)
    f.close()

    return data

def Calc_limit(model, data, poi_values, CL_alpha=0.1):
    """Given a pyhf model and data information, returns the CLs limit."""

    print("CL is " + str(100*(1-CL_alpha)) + "%.")

    obs_limit_single, exp_limits_single, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(data, 
                                                                                                        model, poi_values, 
                                                                                                        level=CL_alpha, return_results=True)
    
    return obs_limit_single, exp_limits_single

def Save_limit(mass, obs_limit, exp_limit, decay_type):
    """Given obs and exp limits, saves them as a csv."""

    fields = ["Mass", "mu"]
    obs_row = [mass, obs_limit]

    fields_exp = ["Mass", "2_sigma_down_mu", "1_sigma_down_mu", "mu", "1_sigma_up_mu", "2_sigma_up_mu"]
    exp_row = [mass]+exp_limit

    savename_obs = f'./limits/observed_{mass}_{decay_type}.csv'
    with open(savename_obs, "w") as s:
        w = csv.writer(s)

        w.writerow(fields)
        w.writerow(obs_row)

    savename_exp = f'./limits/expected_{mass}_{decay_type}.csv'
    with open(savename_exp, "w") as s:
        w = csv.writer(s)
        
        w.writerow(fields_exp)
        w.writerow(exp_row)
    
    

def Single_limit_calc(json_loc, data_loc, mass, decay_type):
    """Single limit calculator."""

    model = Load_json(json_loc)
    full_model = pyhf.Model(model)
    data = Load_data(data_loc)
    poi_values = np.linspace(0.001, 2, 50) #Should have this loaded in. 

    obs, exp = Calc_limit(full_model, data, poi_values, CL_alpha=0.1)

    print("obs:")
    print(obs)
    print("exp:")
    print(exp)

    Save_limit(mass, obs, exp, decay_type)


def Get_json_data_lists():
    """Uses files in the directory to make a list of locations."""
    jsons_path = "./jsons/"
    data_path = "./observed_data/"

    json_files = [f for f in os.listdir(jsons_path) if os.path.isfile(os.path.join(jsons_path, f))]
    data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f))]

    return json_files, data_files

def Get_name_dicts(json_files):
    """Given file string list, return dicts of the masses and decay types."""

    HNL_masses_dict, decay_type_dict = {}, {}

    for json_file in json_files:
        name_split = json_file.split("_")
        HNL_mass = name_split[0]
        decay_type = name_split[1]
        HNL_masses_dict[json_file]=HNL_mass
        decay_type_dict[json_file]=decay_type

    return HNL_masses_dict, decay_type_dict

def Link_json_to_data(json_masses_dict, json_decay_dict, data_masses_dict, data_decay_dict):
    """Given dictionaries with json and data file locations, creates a dict that links them."""
    Linked_dict = {}

    for json_file in json_masses_dict:
        HNL_mass = json_masses_dict[json_file]
        decay_type = json_decay_dict[json_file]

        for data_file in data_masses_dict: #Looping over data dict
            if (data_masses_dict[data_file] == HNL_mass) and (data_decay_dict[data_file] == decay_type):
                Linked_dict[json_file]=data_file
        # data_files = data_masses_dict.keys()[list(data_masses_dict.values()).index(HNL_mass)] #data files with HNL mass
        # if data_files 
        #     data_decay_files = data_files[decay_type]
    
    return Linked_dict

def Main():
    """Main function."""
    json_files, data_files=Get_json_data_lists()
    json_HNL_masses_dict, json_decay_type_dict=Get_name_dicts(json_files)
    data_HNL_masses_dict, data_decay_type_dict=Get_name_dicts(data_files)

    Linked_dict = Link_json_to_data(json_HNL_masses_dict, json_decay_type_dict, data_HNL_masses_dict, data_decay_type_dict)

    for json_name in json_HNL_masses_dict:
        Single_limit_calc(json_name, Linked_dict[json_name], json_HNL_masses_dict[json_name], json_decay_type_dict[json_name])

# Main()

