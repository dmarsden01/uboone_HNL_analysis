import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
import uproot
import uproot3
import functools
import xgboost
import pickle

import Utilities.Constants as Constants
import Utilities.Variables_list as Variables

#Making pickle files
def safe_div(x,y):
    if y == 0.0:
        return 0
    return x / y

def Shannon_entropy_M(Total_entries, Max_entries_in_bin):
    """
    Returns the "Goldilocks" statistic M, which should be between 2 and 3. 
    If smaller it is overbinned, if larger it is underbinned.
    """
    M = np.log2(Total_entries)/(np.log2(Total_entries/Max_entries_in_bin)+1)
    return M

def Get_chi_squared(Obs_hist, Exp_hist, errors_hist):
    """
    Returns the chi squared value for an input observed, expected and error on expected hist. 
    """
    if len(Obs_hist) != len(Exp_hist) or len(Obs_hist)!= len( errors_hist):
        print("Histograms of different lenghts! Exiting.")
        return 0
    Total = 0.0
    for i, obs in enumerate(Obs_hist):
        if Exp_hist[i]==0: continue
        calc = (Obs_hist[i]-Exp_hist[i])**2/errors_hist[i]**2
        Total += calc
    return Total

def Get_resolution(axis, Num_pixels):
    """
    Given the number of pixels in one dimension, returns the number which should be in the other,
    so that the resolution is in the ratio 1920x1080.
    """
    if (axis!="x") and (axis!="y"): 
        print("First argument should be \"x\" or \"y\". Exiting.")
        return 0
    SF = 1920/1080
    if axis=="x":
        print("Given x pixels, calculating y pixels.")
        result = round(Num_pixels/SF)
    if axis=="y":
        print("Given y pixels, calculating x pixels.")
        result = round(Num_pixels*SF)
    print(result)
    return result

def new_create_sample_list(Params):
    if Params["FLATTEN"] == True: Params["Flat_state"] = "flattened"
    else: Params["Flat_state"] = "unflattened"
    if Params["Only_keep_common_DetVar_evs"] == True: Params["Reduced_state"] = "reduced_evs"
    else: Params["Reduced_state"] = "all_evs"

    if Params["only_presel"]: #Only loading variables for preselection plots
        Params["variables_string"] = "Presel_vars"
        Params["variables"] = Variables.Preselection_vars_CRT + Variables.event_vars
        Params["variables_MC"] = Variables.Preselection_vars_CRT_MC + Variables.event_vars
    elif Params["Load_truth_vars"]:
        Params["variables_string"] = "Truth_vars"
        Params["variables"] = Variables.event_vars + Variables.Truth_vars + Variables.multiplicity_vars + ["NeutrinoEnergy2", 'shr_theta_v', 'shr_phi_v', 'trk_theta_v', 'trk_phi_v']
        Params["variables_MC"] = Variables.event_vars + Variables.Truth_vars + Variables.multiplicity_vars + ["NeutrinoEnergy2", 'shr_theta_v', 'shr_phi_v', 'trk_theta_v', 'trk_phi_v']
    else:
        Params["variables_string"] = "my_vars" #Standard variables I use
        Params["variables"] = Variables.Final_variable_list
        Params["variables_MC"] = Variables.Final_variable_list_MC

    if (Params["Run"] == "run1") or (Params["Run"] == "run2a"): Params["current"] = "FHC"
    elif (Params["Run"] == "run3") or (Params["Run"] == "run2b"): Params["current"] = "RHC"
    else: print("Need to choose either \"run1\", \"run3\", \"run2a\" or \"run2b\"")
    
    samples = []
    if Params["Load_standard_bkgs"] == True: samples.extend(["overlay","dirtoverlay","beamoff"])
    if Params["Load_data"] == True: samples.extend(["beamgood"])
    if Params["Load_DetVars"] == True: samples.extend(Constants.Detector_variations)
    if Params["Load_single_file"] == True: samples = [Params["single_file"]]
    #-----Signal-----#
    if Params["Load_lepton_signal"] == True:
        for mass in Constants.HNL_mass_samples:
            samples +=[str(mass) + "_ee"]
    if Params["Load_pi0_signal"] == True:
        for mass in Constants.HNL_mass_pi0_samples:
            samples +=[str(mass) + "_pi0"]
    if Params["Load_lepton_dirac"] == True:
        for mass in Constants.HNL_ee_dirac_mass_samples:
            samples +=[str(mass) + "_ee_dirac"]
    if Params["Load_pi0_dirac"] == True:
        for mass in Constants.HNL_pi0_dirac_mass_samples:
            samples +=[str(mass) + "_pi0_dirac"]
    #-----Signal DetVars-----#
    if Params["Load_Signal_DetVars"] == True:
        # for HNL_mass in Constants.HNL_mass_samples: #For when all detvar samples are made
        if Params["Run"] == "run1":
            for HNL_mass in [50, 100, 150]:
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_ee_"+DetVar]
        if Params["Run"] == "run3":
            for HNL_mass in [2, 10, 20, 50, 100]: #Don't have 150MeV sample yet
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_ee_"+DetVar]
    if Params["Load_pi0_signal_DetVars"] == True:
        if Params["Run"] == "run1": print("There are no pi0 detvar samples for run1.")
        if Params["Run"] == "run3":
            for HNL_mass in [150, 180, 220, 240,245]: #200 is broken for some reason
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_pi0_"+DetVar]
                    
    print(f"Loading these "+Params["Run"]+" samples: " + "\n" + str(samples))
    
    return Params, samples


def Get_all_sample_locs(Params):
    Run = Params["Run"]
    if (Run == "run1")  or (Run == "run2a"): current = "FHC"
    elif (Run == "run3") or (Run == "run2b"): current = "RHC"
    
    if (Run == "run1") or (Run == "run3"):
        sample_loc_dict = {"overlay":f'../NuMI_MC/SLIMMED_neutrinoselection_filt_{Run}_overlay.root',
                           "dirtoverlay":f'../NuMI_MC/neutrinoselection_filt_{Run}_dirt_overlay.root',
                           "beamoff":f'../NuMI_data/neutrinoselection_filt_{Run}_beamoff.root',
                           "beamgood":f'../NuMI_data/neutrinoselection_filt_{Run}_beamon_beamgood.root'}
        
    if (Run == "run2a") or (Run == "run2b"):
        if Run == "run2a": edited_run = "run1"
        if Run == "run2b": edited_run = "run3"
        sample_loc_dict = {"overlay":f'../NuMI_MC/SLIMMED_neutrinoselection_filt_{Run}_overlay.root',
                           "dirtoverlay":f'../NuMI_MC/neutrinoselection_filt_{edited_run}_dirt_overlay.root',
                           "beamoff":f'../NuMI_data/neutrinoselection_filt_{Run}_beamoff.root',
                           "beamgood":f'../NuMI_data/neutrinoselection_filt_{Run}_beamon.root'}
    
    signal_start_str = f'../NuMI_signal/KDAR_dump/sfnues/'
    #Majorana ee
    for mass in Constants.HNL_mass_samples:
        sample_loc_dict.update({str(mass)+"_ee":signal_start_str+f'sfnues_KDAR_dump_{mass}_ee_Umu4_majorana_{current}.root'})
    #Majorana pi0
    for mass in Constants.HNL_mass_pi0_samples:
        sample_loc_dict.update({str(mass)+"_pi0":signal_start_str+f'pi0/sfnues_KDAR_dump_{mass}_pi0_Umu4_majorana_{current}.root'})
    #Dirac ee
    for mass in Constants.HNL_ee_dirac_mass_samples:
        sample_loc_dict.update({str(mass)+"_ee_dirac":signal_start_str+f'sfnues_KDAR_dump_{mass}_ee_Umu4_dirac_{current}.root'})
    #Dirac pi0
    for mass in Constants.HNL_pi0_dirac_mass_samples:
        sample_loc_dict.update({str(mass)+"_pi0_dirac":signal_start_str+f'pi0/sfnues_KDAR_dump_{mass}_pi0_Umu4_dirac_{current}.root'})
    
    #Overlay DetVars
    for Var in Constants.Detector_variations:
        sample_loc_dict.update({Var:f"../NuMI_MC/DetVars/neutrinoselection_filt_{Run}_overlay_{Var}.root"})
    #Signal DetVars
    if Run=="run1": 
        DetVar_ee_masses = Constants.Run1_ee_DetVar_samples
        DetVar_pi0_masses = Constants.Run1_pi0_DetVar_samples
    if Run=="run3":
        DetVar_ee_masses = Constants.Run3_ee_DetVar_samples
        DetVar_pi0_masses = Constants.Run3_pi0_DetVar_samples
    else: 
        DetVar_ee_masses = []
        DetVar_pi0_masses = []
    #ee DetVars
    for mass in DetVar_ee_masses:       
        for DetVar in Constants.Detector_variations:
            sample_loc_dict.update({str(mass)+"_ee_"+DetVar:f"../NuMI_signal/KDAR_dump/sfnues/DetVars/{mass}_{DetVar}_{Run}.root"})
    #pi0 DetVars
    for mass in DetVar_pi0_masses:       
        for DetVar in Constants.Detector_variations:
            sample_loc_dict.update({str(mass)+"_pi0_"+DetVar:f"../NuMI_signal/KDAR_dump/sfnues/pi0/DetVars/{mass}_{DetVar}_{Run}.root"})
        
    return sample_loc_dict

def Edit_Weight_Tune(df_to_Tune): #This is taken from Aditya's code, Owen also has the same in his for overlay and dirt, there is the same block in PELEE code
    df_to_Tune.loc[ df_to_Tune['weightSplineTimesTune'] <= 0, 'weightSplineTimesTune' ] = 1.
    df_to_Tune.loc[ df_to_Tune['weightSplineTimesTune'] == np.inf, 'weightSplineTimesTune' ] = 1.
    df_to_Tune.loc[ df_to_Tune['weightSplineTimesTune'] > 50, 'weightSplineTimesTune' ] = 1.
    df_to_Tune.loc[ np.isnan(df_to_Tune['weightSplineTimesTune']) == True, 'weightSplineTimesTune' ] = 1.
    return df_to_Tune

def MC_weight_branch(df_MC): #Writes a new branch called "weight" including, ppfx, weightSplineTimesTune AND if pi0 are present, scales by pi0 factor
    df_MC["weight"] = df_MC["ppfx_cv"]*df_MC["weightSplineTimesTune"] 
    df_MC.loc[df_MC["npi0"]>0,"weight"] = df_MC["weight"][df_MC["npi0"]>0]*Constants.pi0_scaling_factor #If MC event contains pi0, need to scale down, derived from BNB data
    
def Make_fiducial_vars(df):
    n_pfps = df["n_pfps"].groupby(level="entry").apply(max) #OLD way, gave an error with one sample
    df_placeholder = df.query("trk_sce_start_x_v>-1e4").copy()
    # n_pfps = df_placeholder["n_pfps"].groupby(level="entry").apply(max) #New way as of 31 March
    
    min_x=df_placeholder[["trk_sce_start_x_v","trk_sce_end_x_v"]].min(axis=1).groupby(level="entry").apply(min)
    max_x=df_placeholder[["trk_sce_start_x_v","trk_sce_end_x_v"]].max(axis=1).groupby(level="entry").apply(max)
    min_y=df_placeholder[["trk_sce_start_y_v","trk_sce_end_y_v"]].min(axis=1).groupby(level="entry").apply(min)
    max_y=df_placeholder[["trk_sce_start_y_v","trk_sce_end_y_v"]].max(axis=1).groupby(level="entry").apply(max)
    min_z=df_placeholder[["trk_sce_start_z_v","trk_sce_end_z_v"]].min(axis=1).groupby(level="entry").apply(min)
    max_z=df_placeholder[["trk_sce_start_z_v","trk_sce_end_z_v"]].max(axis=1).groupby(level="entry").apply(max)
    
    del df_placeholder
    df2 = df.copy()
    
    df2["min_x"]=np.array(np.repeat(min_x, np.array(n_pfps)))
    df2["max_x"]=np.array(np.repeat(max_x, np.array(n_pfps)))
    df2["min_y"]=np.array(np.repeat(min_y, np.array(n_pfps)))
    df2["max_y"]=np.array(np.repeat(max_y, np.array(n_pfps)))
    df2["min_z"]=np.array(np.repeat(min_z, np.array(n_pfps)))
    df2["max_z"]=np.array(np.repeat(max_z, np.array(n_pfps)))
    return df2
    
def check_is_truth(df_MC, df_EXT, var): #This is for UNFLATTENED dataframes
    vals_MC = []
    vals_EXT = []
    for i in range(100):
        if isinstance(df_MC[var][i],np.ndarray):
            #print("variable is an array.")
            length_arr_MC = len(df_MC[var][i]) 
            if length_arr_MC == 0:
                vals_MC.append(None)
            else:
                for j in range(len(df_MC[var][i])):
                    vals_MC.append(df_MC[var][i][j])
        if isinstance(df_EXT[var][i],np.ndarray):
            length_arr_EXT = len(df_EXT[var][i]) 
            if length_arr_EXT == 0:
                vals_EXT.append(None)
            else:
                for j in range(len(df_EXT[var][i])):
                    vals_EXT.append(df_EXT[var][i][j])
        else:
            vals_MC.append(df_MC[var][i])
            vals_EXT.append(df_EXT[var][i])
    unique_MC = set(vals_MC)
    unique_EXT = set(vals_EXT)
    print("For the variable \"" + str(var) + "\"" + "\n")

    print("There are " + str(len(unique_MC)) + " unique values in MC.")
    print("There are " + str(len(unique_EXT))+ " unique values in EXT." + "\n")

    print_vals = input("Do you want to see the values of vars? y/n ")
    if print_vals == "y":
        print("MC values are : " + str(unique_MC))
        print("EXT values are: " + str(unique_EXT))
    return 0

def make_unique_ev_id(df): #df must have 'run', 'sub' and 'evt' branches
    if pd.Series(['run', 'sub', 'evt']).isin(df.columns).all():
        rse_list = []
        for entry in df.index: #Looping over all events in the dataframe
            rse = str(df['run'][entry]) + "_" + str(df['sub'][entry]) + "_" + str(df['evt'][entry])
            rse_list.append(rse)
        df['rse_id'] = rse_list #Writing a new branch with the unique event id
        return df.copy()
    else:
        print("Dataframe needs \"run\", \"sub\" and \"evt\" columns.")
        return 0
    
def make_common_evs_df(df_list): #Need to rewrite this so I don't get duplicate columns.
    overlapping_df = functools.reduce(lambda left,right: pd.merge(left, right, on=['rse_id'], how='inner'), df_list)
    print("Length is of common events list is " + str(len(overlapping_df)))
    return overlapping_df


def get_vars(sample, Params):
    MC_samples = ["overlay","dirtoverlay"] + Constants.Detector_variations
    if sample in MC_samples: variables = Params["variables_MC"]
    else: variables = Params["variables"]
    
    return variables
    
def Make_new_vars(df, sample, Params):
    weight_samples = ["overlay","dirtoverlay"] + Constants.Detector_variations 
    #Making fiducial variables
    if sample in Constants.Detector_variations and Params["Run"]=="run3":
        print("skipping fiducal vars, run3 overlay detvar has thetaXZ broken")
        final_file = df.copy()
        del(df)
    elif Params["FLATTEN"] == False:
        final_file = df.copy()
        del(df)
    else:
        final_file = Make_fiducial_vars(df)
        del(df)
    #Making and editing weight, only for MC samples
    if sample in weight_samples:
        Edit_Weight_Tune(final_file)
        MC_weight_branch(final_file)
    #Making "rse_id" column for all files
    make_unique_ev_id(final_file)
    
    return final_file

def get_pkl_savename(sample, loc_pkls, Params):
    save_str = loc_pkls
    Run=Params["Run"]
    split_end = sample.split("_")[-1]
    split_start = sample.split("_")[0]
    if sample in Constants.Detector_variations: save_str += "DetVars/overlay_"
    if Params["Load_Signal_DetVars"] == True: save_str += "Signal_DetVars/"
    if Params["Load_pi0_signal_DetVars"] == True: save_str += "Signal_DetVars/pi0/"
    save_str += f"{sample}_{Run}_"+Params["Flat_state"] 
    if (sample in Constants.Detector_variations) or (split_end in Constants.Detector_variations): save_str += "_"+Params["Reduced_state"]
    
    return save_str
    
def make_filtered_events(df, sample, common_evs, Params):
    if (sample in Constants.Detector_variations) or (Params["Load_Signal_DetVars"] == True) or (Params["Load_pi0_signal_DetVars"] == True):
        # if (common_evs == None) or (len(common_evs)==0): print("No common events loaded!")
        if len(common_evs)==0: print("No common events loaded!")
        if isinstance(common_evs, dict): 
            sample_str = sample.split('_')[0]
            sample_str += "_"
            sample_str += sample.split('_')[1]
            common_evs_list = common_evs[sample_str]
        else: common_evs_list = common_evs
        
        if Params["Only_keep_common_DetVar_evs"] == True:
            filtered = df.loc[(df['rse_id'].isin(common_evs_list['rse_id']))]
            final_file = filtered.copy()
            del(df)
            del(filtered)
        else: 
            final_file = df.copy()
            del(df)
    else: 
        final_file = df.copy()
        del(df)
        
    return final_file

def Check_fiducial_problems(Params, sample):
    """
    For identifying whether certain events would cause issues with creating fiducial vars.
    Returns True if they would.
    """
    if (Params["Run"]=="run2a") and (sample == "beamoff"):
        return True
    elif (Params["Run"]=="run3") and (sample == "WireModThetaXZ"):
        return True
    else: return False

def Check_swtrig_problems(Params, sample):
    """
    For identifying whether certain dataframes are missing \"swtrig_pre\" and \"swtrig_post\".
    Returns True if they would.
    """
    if (Params["Run"]=="run2a" or Params["Run"]=="run2b") and (sample == "beamoff" or sample == "beamgood"):
        return True

    else: return False

def Remove_fiducial_problem_row(df):
    """
    Pass dataframe with n_pfps and \"trk_sce_start_x_v\",\"trk_sce_end_x_v\".
    Returns a dataframe with rows removed which cause issues with fiducial var creation.
    """
    n_pfps = df["n_pfps"].groupby(level="entry").apply(max)
    df_placeholder = df.query("trk_sce_start_x_v>-1e4").copy()
    
    min_x=df_placeholder[["trk_sce_start_x_v","trk_sce_end_x_v"]].min(axis=1).groupby(level="entry").apply(min)
    
    n_pfps_entries = n_pfps.index
    min_x_entries = min_x.index
    
    problem_indices = []
    for index in n_pfps_entries:
        if index not in min_x_entries: problem_indices.append(index)
        
    if len(problem_indices) < 1:
        print("No problem rows")
        return df
    
    for index in problem_indices:
        cleaned = df.drop(labels=374939, axis=0)
    return cleaned

def Create_extra_swtrig_vars(df):
    """
    Copies the values of \"swtrig\" to create \"swtrig_pre\" and \"swtrig_post\".
    Returns the new dataframe.
    """
    df['swtrig_pre'] = df['swtrig'].copy()
    df['swtrig_post'] = df['swtrig'].copy()
    
    df.drop(columns=['swtrig'])
    
    return df

def New_load_and_pkl(samples, sample_loc, loc_pkls, common_evs, Params, save_str=""):
    for sample in samples:
        variables = get_vars(sample, Params)
        root_file = uproot3.open(sample_loc[sample])['nuselection/NeutrinoSelectionFilter']
        if(Check_swtrig_problems(Params, sample)):
            print("swtrig_pre not in vars, creating")
            variables_placeholder=variables.copy()
            variables_placeholder.remove('swtrig_pre')
            variables_placeholder.remove('swtrig_post')
            variables_placeholder+=['swtrig']
        else: variables_placeholder=variables.copy()
        df_file = root_file.pandas.df(variables_placeholder, flatten=Params["FLATTEN"])
        if(Check_swtrig_problems(Params, sample)):
            df_file=Create_extra_swtrig_vars(df_file)
        if(Check_fiducial_problems(Params, sample)):
            df_file=Remove_fiducial_problem_row(df_file)
        new_file = Make_new_vars(df_file, sample, Params) #Edits weight tune, creates weights, event ids and fiducial variables
        final_file = make_filtered_events(new_file, sample, common_evs, Params) #For getting rid of events which aren't in all Detvar samples     

        print("Pickling "+Params["Run"]+f" {sample} file")
        filename = get_pkl_savename(sample, loc_pkls, Params)
        final_file.to_pickle(filename+save_str+".pkl")
        del(final_file)
    print("\n"+"Finished all!")
    
def Make_common_evs(samples, sample_loc, Params):
    if Params["Only_keep_common_DetVar_evs"]==False:
        return None
    if Params["Load_DetVars"] == True: #ONLY loading overlay DetVars
        df_rse_list = []
        for sample in samples:
            root_file = uproot3.open(sample_loc[sample])['nuselection/NeutrinoSelectionFilter']
            df_file = root_file.pandas.df(['run','sub','evt'], flatten=False)
            make_unique_ev_id(df_file)
            df_rse = df_file[['rse_id']]
            df_rse_list.append(df_rse)

        common_evs = make_common_evs_df(df_rse_list)
        return common_evs
    
    elif (Params["Load_Signal_DetVars"] == True) or (Params["Load_pi0_signal_DetVars"] == True): #Loading e+e- OR pi0 DetVars
        masses = []
        for sample in samples:
            sample_split = sample.split("_")
            if len(sample_split) != 3: print("Some samples not signal DetVars")
            name = sample_split[0] + "_" + sample_split[1]
            if name not in masses: masses.append(name)
        print("masses are " + str(masses)) #Only for testing
        common_evs_dict = {}
        for mass in masses:
            df_rse_list = []
            for DetVar in Constants.Detector_variations:
                sample = mass+"_"+DetVar
                root_file = uproot3.open(sample_loc[sample])['nuselection/NeutrinoSelectionFilter']
                df_file = root_file.pandas.df(['run','sub','evt'], flatten=False)
                make_unique_ev_id(df_file)
                df_rse = df_file[['rse_id']]
                df_rse_list.append(df_rse)
            common_evs_dict[mass] = make_common_evs_df(df_rse_list)
        return common_evs_dict
    else: return None


#-----Loading pkls------#
def Load_initial_pkls(samples, Params, loc_pkls, filename):
    """
    For loading initial pkl files. 
    """
    sample_dict = {}

    reduced_str = "_"+Params["Reduced_state"]
    if Params["Load_DetVars"] == True: loc_pkls += "DetVars/overlay_"
    elif Params["Load_Signal_DetVars"] == True: loc_pkls += "Signal_DetVars/"
    elif Params['Load_pi0_signal_DetVars'] == True: loc_pkls += "Signal_DetVars/pi0/"
    else: reduced_str = "" #No reduced string for non-DetVar samples
    
    Run, flat = Params["Run"], Params["Flat_state"]
    
    for sample in samples:
        sample_dict[sample] = pd.read_pickle(loc_pkls+f"{sample}_{Run}_{flat}{reduced_str}{filename}.pkl")
    
    return sample_dict
                
#POT counting
def POT_counter(file): #Takes uproot file
    Total_POT = file["pot"].array().sum()
    return Total_POT


#------------------------#
#------Preselection------#
#------------------------#

def make_unique_events_df(df):
    """
    Input dataframe. Return copy of dataframe with no duplicate events.
    """
    placeholder=df.drop_duplicates(subset=["run","evt","sub"]).copy()
    return placeholder

def count_unique_events(df):
    """
    Input dataframe. Return number of unique events in that dataframe.
    """
    placeholder=df.drop_duplicates(subset=["run","evt","sub"]).copy()
    unique_evs = len(placeholder)
    del placeholder
    return unique_evs

def Preselection_weighted_efficiency(samples, cut_dict): #Need to account for weigthing in overlay and dirt samples
    """
    Input dict of unflattened dataframes.
    Returns the dict of efficiencies and pre-selected samples.
    """
    
    Efficiency_dict, Preselected = {}, {}
    for sample in samples:
        if sample == "overlay" or sample == "dirtoverlay" or sample in Constants.Detector_variations:
            weight = samples[sample]["weight"]
            NumEvs = sum(weight)
        else:
            NumEvs = len(samples[sample]) #Signal, beamoff and beamon don't have weights.
        
        effic_list = [1.0]
        Preselected[sample]=samples[sample].copy()
        for cut in cut_dict.keys():
            # samples[sample]=samples[sample].query(cut_dict[cut])
            Preselected[sample]=Preselected[sample].query(cut_dict[cut])
            if sample == "overlay" or sample == "dirtoverlay" or sample in Constants.Detector_variations:
                weight = Preselected[sample]["weight"]
                Num_selected = sum(weight)
            else:
                Num_selected = len(Preselected[sample])
            effic_list.append(Num_selected/NumEvs)
        Efficiency_dict[sample]=effic_list
        
    return Efficiency_dict, Preselected
        
def Flattened_Preselection_weighted_efficiency(samples, cut_dict, Run): #Need to account for weigthing in overlay and dirt samples
    """
    Input dict of flattened dataframes, dict of cuts and Run.
    Returns the dict of efficiencies and pre-selected samples.
    """
    
    Efficiency_dict, Preselected = {}, {}
    for sample in samples:
        if sample == "overlay" or sample == "dirtoverlay" or sample in Constants.Detector_variations:
            if Run == "run1":NumEvs = Constants.run1_sum_weights[sample] #This is the total BEFORE any preselection
            if Run == "run3":NumEvs = Constants.run3_sum_weights[sample]
        else:
            if Run == "run1":NumEvs = Constants.run1_event_numbers[sample]
            if Run == "run3":NumEvs = Constants.run3_event_numbers[sample]
        
        effic_list = [1.0]
        Preselected[sample]=samples[sample].copy()
        for cut in cut_dict.keys():
            Preselected[sample]=Preselected[sample].query(cut_dict[cut])
            if sample == "overlay" or sample == "dirtoverlay" or sample in Constants.Detector_variations:
                unique_placeholder = make_unique_events_df(Preselected[sample])
                weight = unique_placeholder["weight"]
                Num_selected = sum(weight)
            else:
                unique_placeholder = make_unique_events_df(Preselected[sample])
                Num_selected = len(unique_placeholder)
            effic_list.append(Num_selected/NumEvs)
        Efficiency_dict[sample]=effic_list
        
    return Efficiency_dict, Preselected
    
def Preselection_DetVars(samples, cut_dict): #Not making efficiency plots for DetVars
    """
    Input dict of flattened dataframes and dict of cuts.
    Returns dict of the pre-selected samples.
    """
    Preselected = {}
    for sample in samples:
        Preselected[sample]=samples[sample].copy()
        for cut in cut_dict.keys():
            Preselected[sample]=Preselected[sample].query(cut_dict[cut])
    
    return Preselected


def Get_signal_efficiency_range(Params, Preselection_dict, Efficiency_dict):
    """
    Input the Params and pre-selection efficiency dict.
    Returns the max and min signal efficiency dicts with values for each mass point.
    """
    Preselection_signal_min, Preselection_signal_max = [], []
    min_presel_effic, max_presel_effic = 1.0, 0.0
    
    if Params["Load_lepton_signal"] == True: HNL_masses = Constants.HNL_ee_samples_names
    if Params["Load_pi0_signal"] == True: HNL_masses = Constants.HNL_mass_pi0_samples_names
    if (Params["Load_pi0_signal"] == True) and (Params["Load_lepton_signal"] == True): 
        HNL_masses = Constants.HNL_ee_samples_names+Constants.HNL_mass_pi0_samples_names

    if Params["Load_lepton_dirac"] == True: HNL_masses = Constants.HNL_ee_dirac_names
    if Params["Load_pi0_dirac"] == True: HNL_masses = Constants.HNL_pi0_dirac_names

    for i in range(len(Preselection_dict)+1):
        for HNL_mass in HNL_masses: 
            if Efficiency_dict[HNL_mass][i] > max_presel_effic:
                max_presel_effic = Efficiency_dict[HNL_mass][i]
            if Efficiency_dict[HNL_mass][i] < min_presel_effic:
                min_presel_effic = Efficiency_dict[HNL_mass][i]
        Preselection_signal_max.append(max_presel_effic)
        Preselection_signal_min.append(min_presel_effic)
        max_presel_effic = 0.0
        
    return Preselection_signal_min, Preselection_signal_max


def Get_ev_nums_weights_POT(Run):
    """
    Input the Run "run1" or "run3".
    Returns the total event numbers, sum of event weights and POT normalisation dicts.
    """
    if (Run != "run1") and (Run != "run3"):
        print("Need to select \"run1\" or \"run3\"")
        return 0
    if Run == "run1":
        return Constants.run1_event_numbers, Constants.run1_sum_weights, Constants.run1_POT_scaling_dict
    if Run == "run3":
        return Constants.run3_event_numbers, Constants.run3_sum_weights, Constants.run3_POT_scaling_dict
    
def Get_significance_after_cut(Params, samples, cut):
    """
    Input Params, samples and cut (string used for .query()).
    Returns single number which is signal/sqrt(Total_bkg)
    """
    Preselected = {}
    Num_dict = {}
    
    ev_numbers, ev_sum_weights, POT_norm = Get_ev_nums_weights_POT(Params["Run"])
    
    if Params["Load_lepton_signal"] == True: HNL_masses = Constants.HNL_ee_samples_names
    if Params["Load_pi0_signal"] == True: HNL_masses = Constants.HNL_mass_pi0_samples_names
    if (Params["Load_pi0_signal"] == True) and (Params["Load_lepton_signal"] == True): 
        HNL_masses = Constants.HNL_ee_samples_names+Constants.HNL_mass_pi0_samples_names
    if Params["Load_lepton_dirac"] == True: HNL_masses = Constants.HNL_ee_dirac_names
    if Params["Load_pi0_dirac"] == True: HNL_masses = Constants.HNL_pi0_dirac_names
    
    for sample in samples:
        
        Preselected[sample]=samples[sample].copy()

        Preselected[sample]=Preselected[sample].query(cut)
        if sample == "overlay" or sample == "dirtoverlay" or sample in Constants.Detector_variations:
            unique_placeholder = make_unique_events_df(Preselected[sample]) #Get number of events, not objects
            weight = unique_placeholder["weight"]
            Num_selected = sum(weight)*POT_norm[sample]
        elif sample in HNL_masses:
            unique_placeholder = make_unique_events_df(Preselected[sample]) #Get number of events, not objects
            Num_selected = len(unique_placeholder)
        else:
            unique_placeholder = make_unique_events_df(Preselected[sample]) #Get number of events, not objects
            Num_selected = len(unique_placeholder)*POT_norm[sample]

        Num_dict[sample]=Num_selected
    
    Tot_bkg = Num_dict["overlay"]+Num_dict["dirtoverlay"]+Num_dict["beamoff"]
    
    max_Num = 0
    min_Num = 1e7
    for HNL_mass in HNL_masses:
        if Num_dict[HNL_mass] > max_Num: max_Num = Num_dict[HNL_mass]
        if Num_dict[HNL_mass] < min_Num: min_Num = Num_dict[HNL_mass]
    
    Significance_max = max_Num/np.sqrt(Tot_bkg)
    Significance_min = min_Num/np.sqrt(Tot_bkg)
    
    return Significance_max, Significance_min

def Significance_scan(Params, samples_dict, cut_string_start, scan_start, scan_end, numsteps=20):
    """
    Scans through cut values to find the maximum significance.
    """
    Significance_dict_max, Significance_dict_min = {}, {}
    var_scan = np.linspace(scan_start, scan_end, numsteps)
    cut_scan = []
    for var in var_scan:
        cut_scan.append(cut_string_start + str(var))
        
    for i, cut in enumerate(cut_scan):
        var = var_scan[i]
        Significance_max, Significance_min = Get_significance_after_cut(Params, samples_dict, cut)
        Significance_dict_max[var] = Significance_max
        Significance_dict_min[var] = Significance_min
        
    var_sig_max = max(Significance_dict_max, key=Significance_dict_max.get)
    var_sig_min = max(Significance_dict_min, key=Significance_dict_min.get)
    
    print(f"Maximum signal max significance is at {var_sig_max}")
    print(f"Minimum signal max significance is at {var_sig_min}")
    
    return Significance_dict_max, Significance_dict_min, var_sig_max, var_sig_min


def Print_efficiency_numbers(Params, Preselected_dict, Efficiency_dict):
    """
    Input the Params and efficiency dict.
    Prints the efficiency and event numbers for putting into a table.
    """
    print(Params["Run"])
    for sample in Efficiency_dict:
        print(f"{sample} efficiency is " + str(Efficiency_dict[sample][-1]*100) + "%")
    
    ev_numbers, ev_sum_weights, POT_norm = Get_ev_nums_weights_POT(Params["Run"])
    
    Num_selected_dict = {}
    for sample in Preselected_dict:
        if sample == "overlay" or sample == "dirtoverlay" or sample in Constants.Detector_variations:
            unique_placeholder = make_unique_events_df(Preselected_dict[sample])
            weight = unique_placeholder["weight"]
            Num_selected = sum(weight)*POT_norm[sample]
            Num_selected_dict[sample]=Num_selected
        else:
            unique_placeholder = make_unique_events_df(Preselected_dict[sample])
            Num_selected = len(unique_placeholder)*POT_norm[sample]
            Num_selected_dict[sample]=Num_selected
        print(f"{sample} " + str(Num_selected))
        
    selected_bkg_sum = Num_selected_dict["overlay"]+Num_selected_dict["dirtoverlay"]+Num_selected_dict["beamoff"]
    initial_bkg_sum  = ev_sum_weights["overlay"]*POT_norm["overlay"]+ev_sum_weights["dirtoverlay"]*POT_norm["dirtoverlay"]+ev_numbers["beamoff"]*POT_norm["beamoff"]
    print("Sum of bkgs: " + str(selected_bkg_sum))
    print("Sum of bkgs effic: " + str((selected_bkg_sum/initial_bkg_sum)*100) + "%")

    if "beamgood" in Num_selected_dict.keys():
        print("Data/prediction: " + str(Num_selected_dict["beamgood"]/selected_bkg_sum))
        
def Get_effic_wrt_previous(Params, Preselection_dict, Efficiency_dict):
    """
    Input Params and dict of absolute efficiencies.
    Returns a dict of efficiencies wrt previous cut and list of signal min and max for that.
    """
    effic_wrt_prev = {}
    lowest_signal_wrt_prev, highest_signal_wrt_prev = {}, {}
    for sample in Efficiency_dict:
        effic_list = []
        for i, effic in enumerate(Efficiency_dict[sample]):
            if i==0: continue
            effic_list.append(Efficiency_dict[sample][i]/Efficiency_dict[sample][i-1])
        effic_wrt_prev[sample] = effic_list

    lowest_signal_wrt_prev, highest_signal_wrt_prev = [], []

    if Params["Load_lepton_signal"] == True: HNL_masses = Constants.HNL_ee_samples_names
    if Params["Load_pi0_signal"] == True: HNL_masses = Constants.HNL_mass_pi0_samples_names
    if (Params["Load_pi0_signal"] == True) and (Params["Load_lepton_signal"] == True): 
        HNL_masses = Constants.HNL_ee_samples_names+Constants.HNL_mass_pi0_samples_names

    if Params["Load_lepton_dirac"] == True: HNL_masses = Constants.HNL_ee_dirac_names
    if Params["Load_pi0_dirac"] == True: HNL_masses = Constants.HNL_pi0_dirac_names

    for i in range(len(Preselection_dict)):
        min_presel_effic, max_presel_effic = 1.0, 0.0
        for HNL_mass in HNL_masses: 
            if effic_wrt_prev[HNL_mass][i] > max_presel_effic:
                max_presel_effic = effic_wrt_prev[HNL_mass][i]
            if effic_wrt_prev[HNL_mass][i] < min_presel_effic:
                min_presel_effic = effic_wrt_prev[HNL_mass][i]
        highest_signal_wrt_prev.append(max_presel_effic)
        lowest_signal_wrt_prev.append(min_presel_effic)
    
    return effic_wrt_prev, lowest_signal_wrt_prev, highest_signal_wrt_prev
        
def Save_preselected_pkls(Prepared_dict, Params, loc_pkls, save_str):
    """
    Pass the \"Prepared dict\" to be saved as Preselected .pkl files.
    Saves the prselected .pkls in the location specified.
    """
    Run=Params["Run"]
    for sample in Prepared_dict:
        print("Saving "+Params["Run"]+f" Preselected {sample} .pkl")
        start_str = loc_pkls
        end_str = save_str
        split_end = sample.split("_")[-1]
        split_start = sample.split("_")[0]
        if sample in Constants.Detector_variations: #overlay DetVars
            start_str += f"DetVars/Preselected_overlay_"
            end_str = Params["Reduced_state"] + "_" + end_str
        elif Params["Load_Signal_DetVars"] == True: #e+e- DetVars
            start_str += f"Signal_DetVars/Preselected_"
            end_str = Params["Reduced_state"] + "_" + end_str
        elif Params["Load_pi0_signal_DetVars"] == True: #pi0 DetVars
            start_str += f"Signal_DetVars/pi0/Preselected_"
            end_str = Params["Reduced_state"] + "_" + end_str   
        # elif Params["Load_pi0_signal"] == True:  #pi0 samples
        elif split_end == "pi0":  #pi0 samples
            start_str += f"pi0_selection/Preselected_"
        else:
            start_str += f"Preselected_"
        
        # print("Saving as " + f"{start_str}{sample}_{Run}_"+Params["Flat_state"]+f"_{end_str}.pkl")
        Prepared_dict[sample].to_pickle(f"{start_str}{sample}_{Run}_"+Params["Flat_state"]+f"_{end_str}.pkl")

#------------------------#
#-------BDT training-----#
#------------------------#

def Remove_high_trk_score_objects(samples, threshold=0.97):
    """
    The pre-selection used to have a track score cut.
    Now I just do it after the standard pre-selection, before saving the .pkl files
    """
    score_cut_dict = {}
    print(f"Removing tracks with trk_score_v > {threshold}")
    for sample in samples:
        score_cut_dict[sample] = samples[sample].copy()
        score_cut_dict[sample] = score_cut_dict[sample].query(f"trk_score_v < {threshold}")
    
    return score_cut_dict

def Fixed_Prepare_dfs_for_xgb(samples): 
    """
    Takes a dictionary of dataframes.
    Returns a dictionary with non-reco vals changed to -9999.
    """
    cleaned_dict = {}
    print("Changing non-reco values to -9999")
    for sample in samples:
        cleaned_dict[sample] = samples[sample].copy()
        value = -1e15
        new_value = -9999
        first_entry = cleaned_dict[sample].index[0]
        for variable in cleaned_dict[sample].keys():
            if isinstance(cleaned_dict[sample][variable][first_entry], (int,float,np.int32,np.float32,np.uint32)):
                if(len(cleaned_dict[sample].loc[cleaned_dict[sample][variable] < value]) > 0):
                    cleaned_dict[sample].loc[(cleaned_dict[sample][variable] < value), variable] = new_value #Gets rid of hugely negative values
                if(len(cleaned_dict[sample].loc[cleaned_dict[sample][variable] == -1.0]) > 0):
                    cleaned_dict[sample].loc[(cleaned_dict[sample][variable] == -1.0), variable] = new_value #Gets rid of exactly -1.0
                if(len(cleaned_dict[sample].loc[cleaned_dict[sample][variable] == np.nan]) > 0):
                    cleaned_dict[sample].loc[(cleaned_dict[sample][variable] == np.nan), variable] = new_value #Gets rid of nans
                if(len(cleaned_dict[sample].loc[cleaned_dict[sample][variable] == np.inf]) > 0):
                    cleaned_dict[sample].loc[(cleaned_dict[sample][variable] == np.inf), variable] = new_value #Gets rid of infs
                if(len(cleaned_dict[sample].loc[cleaned_dict[sample][variable] > (-1*value)]) > 0):
                    cleaned_dict[sample].loc[(cleaned_dict[sample][variable] > (-1*value)), variable] = new_value #Gets rid of hugely positive values
                if variable in ["shr_energy_tot", "trk_energy", "trk_energy_hits_tot", "trk_energy_tot"]:
                    cleaned_dict[sample].loc[(cleaned_dict[sample][variable] == 0.0), variable] = new_value #Sets non-reco trk, shr energies to -9999
                if variable == "subcluster":
                    cleaned_dict[sample].loc[(cleaned_dict[sample][variable] > 2e8), variable] = new_value #Sets non-reco subcluster to -9999
 
    return cleaned_dict

# def Sophisticated_dfs_for_xgb(df): #Requires a minimum of 2 reconstructed objects
#     #Take highest E object
#     print("Write this")
#     variable = 'pfnplanehits_Y'
    # df.loc
    #Look at 2nd highest E object, if within x cm save these two
    #If not look for next highest E object and repeat. 
    #If exhausted of objects, remove event
    
    #Make event-wise variables to feed into BDT
    #Save opening angle
    #Save theta, phi of highest E object
    #Save theta, phi of lower E object
    #Save length of highest E object
    #Save length of lower E object
    #Save invariant mass of 2 objects
    #Save total E of both objects

def only_keep_highest_E(samples):
    """
    Takes a dictionary of dataframes.
    Returns a dictionary with only the highest energy objects remaining.
    """
    highest_E_samples = {}
    print("Only keeping object with highest associated pfnplanehits_Y")
    for sample in samples:
        highest_E_samples[sample] = samples[sample].copy()
        highest_E_samples[sample]["highest_E"]=highest_E_samples[sample]['pfnplanehits_Y'].groupby("entry").transform(max) == highest_E_samples[sample]['pfnplanehits_Y']
        highest_E_samples[sample] = highest_E_samples[sample].query("highest_E").copy()
    return highest_E_samples

def create_test_samples_list(Params): #Returns the list of samples to run over
    samples = [] #A list of all the samples which will be loaded and pickled
    if Params["Load_standard"] == True:
        samples.extend(["overlay","dirtoverlay","beamoff"])
    if Params["Load_lepton_signal"] == True:
        samples.extend(Constants.HNL_ee_samples_names)
    if Params["Load_pi0_signal"] == True:
        samples.extend(Constants.HNL_mass_pi0_samples_names)
        # for HNL_mass in Constants.HNL_mass_pi0_samples_names:
            # samples+=[str(HNL_mass)+"_pi0"]
    if Params["Load_lepton_dirac"] == True: samples.extend(Constants.HNL_ee_dirac_names)
    if Params["Load_pi0_dirac"] == True: samples.extend(Constants.HNL_pi0_dirac_names)
    
    if Params["Load_DetVars"] == True: #This is overlay DetVars
        samples.extend(Constants.Detector_variations)
    if Params["Load_Signal_DetVars"] == True:
        # for HNL_mass in Constants.HNL_mass_samples: #For when all detvar samples are made
        if Params["Run"] == "run1": masses = Constants.Run1_ee_DetVar_samples
        if Params["Run"] == "run3": masses = Constants.Run3_ee_DetVar_samples
        for HNL_mass in masses: #[50, 100, 150, 180, 200] #Now use names
            for DetVar in Constants.Detector_variations:
                samples+=[str(HNL_mass)+"_"+DetVar]
    if Params['Load_pi0_signal_DetVars'] == True:
        if Params["Run"] == "run1": masses = Constants.Run1_pi0_DetVar_samples
        if Params["Run"] == "run3": masses = Constants.Run3_pi0_DetVar_samples
        for HNL_mass in masses: #[50, 100, 150, 180, 200] #Now use names
            for DetVar in Constants.Detector_variations:
                samples+=[str(HNL_mass)+"_"+DetVar]
                
    if Params["Load_data"] == True:
        samples.extend(["beamgood"])
        
    if Params["Load_single_file"] == True:
        samples.extend([Params["single_file"]])
        
    print(f"Loading these "+Params["Run"]+" samples: " + "\n" + str(samples))
    
    return samples

def create_sig_detsys_samples_list(Params): #Returns the list of samples to run over
    samples = [] #A list of all the samples which will be loaded and pickled
    masses = []
    # for HNL_mass in Constants.HNL_mass_samples:
    if Params["Load_lepton_signal"] == True:
        if Params["Run"] == "run1":
            for HNL_mass in [150]: #While I only have 150MeV sample
                masses+=[HNL_mass]
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_ee_"+DetVar]
        if Params["Run"] == "run3":
            for HNL_mass in [2, 10, 20, 50, 100]: #Don't have 150MeV sample yet
                masses+=[HNL_mass]
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_ee_"+DetVar]
    if Params["Load_pi0_signal"] == True:
        if Params["Run"] == "run1":
            print("There are no run1 pi0 detector variation samples")
        if Params["Run"] == "run3":
            for HNL_mass in [150,180,220,240,245]: #200 broken for some reason
                masses+=[HNL_mass]
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_pi0_"+DetVar]
        
    print(f"Loading these "+Params["Run"]+" samples: " + "\n")
    print(samples)
    
    if Params["Use_logit"] == True:
        Params["logit_str"] = "logit"
    else: Params["logit_str"] = "standard"
    
    return samples, masses

def Split_samples(Presel_overlay, signal_samples_dict, Presel_EXT, Presel_dirt, 
                  overlay_train_frac=0.7, signal_train_frac=0.7, EXT_train_frac=0.3, dirt_train_frac=0.3,
                  print_vals=True):
    """
    Input pkl files and train_vs_test fractions.
    Return the split samples.
    """
    overlay_train = Presel_overlay[:int(len(Presel_overlay)*overlay_train_frac)]
    overlay_test = Presel_overlay[int(len(Presel_overlay)*overlay_train_frac):]
    
    EXT_train = Presel_EXT[:int(len(Presel_EXT)*EXT_train_frac)]
    EXT_test = Presel_EXT[int(len(Presel_EXT)*EXT_train_frac):]
    
    dirt_train = Presel_dirt[:int(len(Presel_dirt)*dirt_train_frac)]
    dirt_test = Presel_dirt[int(len(Presel_dirt)*dirt_train_frac):]
    
    signal_train_dict, signal_test_dict = {}, {}
    
    if(print_vals):
        print(f"Length overlay train {len(overlay_train)}")
        print(f"Length overlay test {len(overlay_test)}")
        print(f"Length EXT train {len(EXT_train)}")
        print(f"Length EXT test {len(EXT_test)}")
        print(f"Length dirt train {len(dirt_train)}")
        print(f"Length dirt test {len(dirt_test)}")
    
    for HNL_mass in signal_samples_dict:
        signal_train_dict[HNL_mass] = signal_samples_dict[HNL_mass][:int(len(signal_samples_dict[HNL_mass])*signal_train_frac)]
        signal_test_dict[HNL_mass] = signal_samples_dict[HNL_mass][int(len(signal_samples_dict[HNL_mass])*signal_train_frac):]
        if(dirt_train_frac):
            print(f"Length {HNL_mass} train {len(signal_train_dict[HNL_mass])}")
            print(f"Length {HNL_mass} test {len(signal_test_dict[HNL_mass])}")
    
    split_dict = {"overlay_train":overlay_train, "overlay_test":overlay_test,
                  "signal_train_dict":signal_train_dict, "signal_test_dict":signal_test_dict,
                  "EXT_train":EXT_train, "EXT_test":EXT_test,
                  "dirt_train":dirt_train, "dirt_test":dirt_test}
        
    return split_dict

def Save_test_pkls(Params, loc_pkls, save_str, overlay_test, signal_test_dict, EXT_test = [], dirt_test = []):
    """
    Input Params, save_str, overlay test df, signal_test_dict, and EXT test sample if using.
    Saves the dataframes as .pkl files.
    """
    start_str = loc_pkls
    if Params["Load_pi0_signal"] == False: start_str+="BDT_Test_dfs/"
    if Params["Load_pi0_signal"] == True: start_str+="BDT_Test_dfs/pi0_selection/"
    
    print(f"Pickling overlay test sample")
    overlay_test.to_pickle(start_str+"Test_overlay_"+Params["Run"]+f"_flattened{save_str}.pkl")
    
    for HNL_mass in signal_test_dict:
        print(f"Pickling {HNL_mass} test sample")
        signal_test_dict[HNL_mass].to_pickle(start_str+f"Test_{HNL_mass}_"+Params["Run"]+f"_flattened{save_str}.pkl")
        
    if Params["EXT_in_training"] == True: 
        print(f"Pickling beamoff test sample")
        EXT_test.to_pickle(start_str+"Test_beamoff_"+Params["Run"]+f"_flattened{save_str}.pkl")
        
    if Params["dirt_in_training"] == True: 
        print(f"Pickling dirt test sample")
        dirt_test.to_pickle(start_str+"Test_dirtoverlay_"+Params["Run"]+f"_flattened{save_str}.pkl")
        
    if Params["EXT_in_training"] == False: print("Not saving beamoff test, as not using in training.") 
    if Params["dirt_in_training"] == False: print("Not saving dirt test, as not using in training.")
        

def Make_train_labels_and_dicts(Params, bdt_vars, overlay_train, EXT_train, dirt_train, signal_train_dict):
    """
    Input Params, bdt variables, training samples.
    Return dicts of labels indicating if the event is signal (1) or bkg (0).
    """
    combined_dict, labels_dict = {}, {}
    
    for HNL_mass in signal_train_dict:
        # labels_dict[HNL_mass] = [1]*len(signal_train_dict[HNL_mass][bdt_vars]) + [0]*len(overlay_train[bdt_vars])
        combined_dict[HNL_mass] = pd.concat([signal_train_dict[HNL_mass][bdt_vars], overlay_train[bdt_vars]])
        labels_dict[HNL_mass] = [1]*len(signal_train_dict[HNL_mass]) + [0]*len(overlay_train)
        
    bkg_train = overlay_train.copy()
        
    if Params["EXT_in_training"] == True:
        for HNL_mass in signal_train_dict:
            # labels_dict[HNL_mass] = labels_dict[HNL_mass] + [0]*len(EXT_train[bdt_vars])
            combined_dict[HNL_mass] = pd.concat([combined_dict[HNL_mass][bdt_vars], EXT_train[bdt_vars]])
            labels_dict[HNL_mass] = labels_dict[HNL_mass] + [0]*len(EXT_train)
        bkg_train=pd.concat([overlay_train, EXT_train])
        
    if Params["dirt_in_training"] == True:
        for HNL_mass in signal_train_dict:
            # labels_dict[HNL_mass] = labels_dict[HNL_mass] + [0]*len(EXT_train[bdt_vars])
            combined_dict[HNL_mass] = pd.concat([combined_dict[HNL_mass][bdt_vars], dirt_train[bdt_vars]])
            labels_dict[HNL_mass] = labels_dict[HNL_mass] + [0]*len(dirt_train)
        bkg_train=pd.concat([bkg_train, dirt_train])
    
    return combined_dict, labels_dict, bkg_train

def Prepare_for_xgb(Params, bdt_vars, combined_train_dict, signal_train_dict, bkg_train,  signal_test_dict, bkg_test,
                    labels_train_dict, missing=-9999.0):
    """
    Input training dict, test dict, labels and xgboost parameters.
    Returns the DMatrix forms of the dataframes for training.
    """
    xgb_train_bkg = xgboost.DMatrix(bkg_train[bdt_vars],label=[0]*len(bkg_train[bdt_vars]), missing=missing, feature_names=bdt_vars)
    xgb_test_bkg = xgboost.DMatrix(bkg_test[bdt_vars], label=[0]*len(bkg_test[bdt_vars]), missing=missing, feature_names=bdt_vars)
    
    xgb_train_dict, xgb_sig_test_dict, xgb_sig_train_dict = {}, {}, {}
    
    for HNL_mass in signal_train_dict:
        xgb_train_dict[HNL_mass] = xgboost.DMatrix(combined_train_dict[HNL_mass][bdt_vars], label=labels_train_dict[HNL_mass], 
                                                   missing=missing, feature_names=bdt_vars)
        
        xgb_sig_test_dict[HNL_mass] = xgboost.DMatrix(signal_test_dict[HNL_mass][bdt_vars], label=[1]*len(signal_test_dict[HNL_mass][bdt_vars]), 
                                                      missing=missing, feature_names=bdt_vars)
        
        xgb_sig_train_dict[HNL_mass] = xgboost.DMatrix(signal_train_dict[HNL_mass][bdt_vars],label=[1]*len(signal_train_dict[HNL_mass][bdt_vars]),
                                                  missing=missing, feature_names=bdt_vars)
        
        
    return xgb_train_dict, xgb_sig_train_dict, xgb_sig_test_dict, xgb_train_bkg, xgb_test_bkg
        
def Train_BDTs(Params, bdt_vars, BDT_name, xgb_train_dict, xgb_sig_test_dict, xgb_test_bkg, xgb_param, xgb_combined_test_dict,
               progress, num_round = 150, early_stop=30, missing=-9999.0):
    """
    Input training dict, test dict, labels and xgboost parameters.
    Saves the BDT models as .json files.
    """
    for HNL_mass in xgb_train_dict:
        # watchlist = [(xgb_train_dict[HNL_mass], 'train'), (xgb_sig_test_dict[HNL_mass], 'test_sig'), (xgb_test_bkg,'test_bkg')] for logloss
        watchlist = [(xgb_train_dict[HNL_mass], 'train'), (xgb_combined_test_dict[HNL_mass], 'test_combined')]
        print(f"Training {HNL_mass} BDT" + "\n")
        # bdt = xgboost.train(xgb_param, xgb_train_dict[HNL_mass], num_round, watchlist, evals_result=progress, verbose_eval=False)
        individual_progress=dict()
        bdt = xgboost.train(params=xgb_param, dtrain=xgb_train_dict[HNL_mass], num_boost_round=num_round, evals=watchlist, 
                            early_stopping_rounds=early_stop, evals_result=individual_progress, verbose_eval=False)
        progress[HNL_mass] = individual_progress
        start_save_str = "bdts/"
        if Params["Load_pi0_signal"] == True: start_save_str += "pi0_selection/"
        
        bdt.save_model(start_save_str+Params["Run"]+f"_{HNL_mass}{BDT_name}.json") #Saving with native function for compatibility, but doesn't save params used properly
        pickle.dump(bdt, open(start_save_str+Params["Run"]+f"_{HNL_mass}{BDT_name}.pkl", "wb")) #DOES correctly save the params used in training
        
            
def Get_sample_norms(Params, sig_names, sig_train, overlay_train, EXT_train, dirt_train):
    """
    Input Params, signal names list and the training fractions used.
    Returns a dict of sample_norms NOT WEIGHTED by event weights. 
    """
    SF_test_sig = 1.0/(1-sig_train)
    SF_test_overlay = 1.0/(1-overlay_train)
    SF_EXT = 1.0/(1-EXT_train)
    SF_dirt = 1.0/(1-dirt_train)

    if Params["Run"] == "run1": POT_scale_dict = Constants.run1_POT_scaling_dict
    elif Params["Run"] == "run3": POT_scale_dict = Constants.run3_POT_scaling_dict
        
    overlay_scale = POT_scale_dict["overlay"]*SF_test_overlay
    beamoff_scale = POT_scale_dict["beamoff"]*SF_EXT
    dirtoverlay_scale = POT_scale_dict["dirtoverlay"]*SF_dirt
    
    sample_norms={'overlay':overlay_scale,
                  'dirtoverlay':dirtoverlay_scale,
                  'beamoff':beamoff_scale}
    
    for HNL_mass in sig_names:
        sample_norms[HNL_mass]=SF_test_sig
            
    return sample_norms

def Get_weighted_sample_norms(Params, sample_dict, sig_names, sig_train, overlay_train, EXT_train, dirt_train):
    """
    Input Params, dict of sample dataframes with weights (if applicable) and the training fractions used.
    Returns a dict of sample_norms. 
    """
    SF_test_sig = 1.0/(1-sig_train)
    SF_test_overlay = 1.0/(1-overlay_train)
    SF_EXT = 1.0/(1-EXT_train)
    SF_dirt = 1.0/(1-dirt_train)

    if Params["Run"] == "run1": POT_scale_dict = Constants.run1_POT_scaling_dict
    elif Params["Run"] == "run3": POT_scale_dict = Constants.run3_POT_scaling_dict
        
    overlay_scale = POT_scale_dict["overlay"]*SF_test_overlay
    beamoff_scale = POT_scale_dict["beamoff"]*SF_EXT
    dirtoverlay_scale = POT_scale_dict["dirtoverlay"]*SF_dirt #dirt not used in training
    
    sample_norms={'overlay_test':np.array(sample_dict['overlay']["weight"]*overlay_scale),
                  'dirtoverlay':np.array(sample_dict['dirtoverlay']["weight"]*dirtoverlay_scale),
                  'beamoff':np.ones(len(sample_dict['beamoff']['n_pfps']))*beamoff_scale}
    
    for HNL_mass in sig_names:
        signal_scale_list = np.ones(len(sample_dict[HNL_mass]['n_pfps']))*SF_test_sig
        sample_norms[HNL_mass]=signal_scale_list
            
    return sample_norms

def SaveToRoot(nbins,xlims,bkg_overlay,bkg_dirt,bkg_EXT,sig,data,fileName='test.root'):
    nBins = nbins
    binLimits = xlims
  ### Save files 
    rFile = ROOT.TFile(f'bdt_output/{fileName}','RECREATE')
    tData1 = ROOT.TH1F("Signal","Signal",nBins,binLimits[0],binLimits[1])
    for i in range(nBins):
        tData1.SetBinContent(i+1,sig['hist'][i])
        tData1.SetBinError(i+1,sig['err'][i])
    tData2 = ROOT.TH1F("bkg_overlay","bkg_overlay",nBins,binLimits[0],binLimits[1])
    for i in range(nBins):
        tData2.SetBinContent(i+1,bkg_overlay['hist'][i])
        tData2.SetBinError(i+1,bkg_overlay['err'][i])
    tData3 = ROOT.TH1F("bkg_dirt","bkg_dirt",nBins,binLimits[0],binLimits[1])
    for i in range(nBins):
        tData3.SetBinContent(i+1,bkg_dirt['hist'][i])
        tData3.SetBinError(i+1,bkg_dirt['err'][i])
    tData4 = ROOT.TH1F("bkg_EXT","bkg_EXT",nBins,binLimits[0],binLimits[1])
    for i in range(nBins):
        tData4.SetBinContent(i+1,bkg_EXT['hist'][i])
        tData4.SetBinError(i+1,bkg_EXT['err'][i])
    tData5 = ROOT.TH1F("Data","Data",nBins,binLimits[0],binLimits[1])
    for i in range(nBins):
        tData5.SetBinContent(i+1,data['hist'][i])
        tData5.SetBinError(i+1,data['err'][i])
    rFile.Write()
    rFile.Close()

#Logistic transformations for BDT scores
def logit(x): #The "logit" function, also known as the inverse of the logistic function
    return np.log(x/(1-x))

def invlogit(x): #The "logistic" function, also known as the inverse of the "logit" function
    return np.exp(x)/(1+np.exp(x))

#Systematic errors from reweighting functions

def All_reweight_err(df, Multisim, var_name, BINS, x_range, Norm):
    results_dict = {}
    n_bins = len(BINS)-1

    Nuniverse = Constants.Multisim_univs[Multisim]
    n_tot = np.empty([Nuniverse, n_bins])
    n_cv_tot = np.empty(n_bins)
    n_tot.fill(0)
    n_cv_tot.fill(0)

    variable = df[var_name] #The BDT output score
    syst_weights = df[Multisim] #An array of length of the number of events, each entry is an array of length Nunivs
    spline_fix_cv  = df["weight"]*Norm
    spline_fix_var = df["weight"]*Norm

    s = syst_weights
    df_weights = pd.DataFrame(s.values.tolist())
    n_cv, bins = np.histogram(variable, range=x_range, bins=BINS, weights=spline_fix_cv)
    n_cv_tot += n_cv

    if(Multisim == "weightsGenie"): #special treatment as ["weightSplineTimesTune"] is included in genie weights
        if not df_weights.empty:
            for i in range(Nuniverse):
                weight = df_weights[i].values / 1000.
                weight[weight == 1]= df["weightSplineTimesTune"].iloc[weight == 1]
                weight[np.isnan(weight)] = df["weightSplineTimesTune"].iloc[np.isnan(weight)]
                weight[weight > 50] = df["weightSplineTimesTune"].iloc[weight > 50] # why 30 not 50?
                weight[weight <= 0] = df["weightSplineTimesTune"].iloc[weight <= 0]
                weight[weight == np.inf] = df["weightSplineTimesTune"].iloc[weight == np.inf]

                n, bins = np.histogram(variable, 
                                       weights=np.nan_to_num(weight*spline_fix_var/df["weightSplineTimesTune"]), range=x_range, bins=BINS)
                n_tot[i] += n

    if(Multisim == "weightsPPFX"): #special treatment as ["PPFXPcv"] is included in ppfx weights
        if not df_weights.empty:
            for i in range(Nuniverse):
                weight = df_weights[i].values / 1000.
                weight[weight == 1]= df["ppfx_cv"].iloc[weight == 1]
                weight[np.isnan(weight)] = df["ppfx_cv"].iloc[np.isnan(weight)]
                weight[weight > 100] = df["ppfx_cv"].iloc[weight > 100]
                weight[weight < 0] = df["ppfx_cv"].iloc[weight < 0]
                weight[weight == np.inf] = df["ppfx_cv"].iloc[weight == np.inf]

                n, bins = np.histogram(variable, weights=weight*np.nan_to_num(spline_fix_var/df["ppfx_cv"]), range=x_range, bins=BINS)
                n_tot[i] += n

    if(Multisim == "weightsReint"):
        if not df_weights.empty:
            for i in range(Nuniverse):
                weight = df_weights[i].values / 1000.
                weight[np.isnan(weight)] = 1
                weight[weight > 100] = 1
                weight[weight < 0] = 1
                weight[weight == np.inf] = 1
                n, bins = np.histogram(variable, weights=weight*spline_fix_var, range=x_range, bins=BINS)
                n_tot[i] += n
    cov = np.empty([len(n_cv), len(n_cv)])
    cov.fill(0)

    for n in n_tot:
        for i in range(len(n_cv)):
            for j in range(len(n_cv)):
                cov[i][j] += (n[i] - n_cv_tot[i]) * (n[j] - n_cv_tot[j])

    cov /= Nuniverse
    results_dict[Multisim] = [cov,n_cv_tot,n_tot,bins]
    return results_dict

def make_stat_err(var, bins, weights_times_SF):
    """
    Takes a variable, histogram bins and the full weights (can include multiplied by normalisation).
    Returns the associated Poisson error on each bin.
    """
    hist_unweighted = np.histogram(var,bins=bins)[0]
    hist_weighted = np.histogram(var,bins=bins,weights=weights_times_SF)[0]
    Total_SF = hist_weighted/hist_unweighted
    stat_err=np.sqrt(hist_unweighted)*Total_SF
    return stat_err

def Get_signal_name_type(Params):
    """
    Given script parameters will returen "ee", "pi0", "ee_dirac" or "pi0_dirac".
    """
    if Params["Load_lepton_signal"] == True:
        return "ee"
    elif Params["Load_pi0_signal"] == True: 
        return "pi0"
    elif Params["Load_lepton_dirac"] == True: 
        return "ee_dirac"
    elif Params["Load_pi0_dirac"] == True: 
        return "pi0_dirac"
    else:
        print("Sample types not recognised, returning 1.")
        return 1
    
#Limit setting functions    
def pyhf_params(Params):
    if Params["Stats_only"] == True:
        print("Calculating stats-only limit.")
    elif Params["Use_flat_sys"] == True:
        print("Using FLAT systematic uncertainty on signal and background")
        perc_overlay = Params["Flat_bkg_overlay_frac"]*100
        perc_dirt = Params["Flat_bkg_dirt_frac"]*100
        print("With " + str(perc_overlay) + "% on overlay, and " + str(perc_dirt) + "% on dirt.")
        perc_signal_KDAR = Params["Signal_flux_error"]
        perc_signal_detvar = Params["Flat_sig_detvar"]
        total_sig_err = np.sqrt(perc_signal_KDAR**2 + perc_signal_detvar**2)
        print("With " + str(total_sig_err*100) + "% on all signal")
    else: 
        perc_dirt = Params["Flat_bkg_dirt_frac"]*100
        print(f"Using fully evaluated systematic uncertainty for background. Dirt error is {perc_dirt}%.")
        perc_flux = Params["Signal_flux_error"]*100
        print(f"Using fully evaluated systematic uncertainty for signal. Using {perc_flux}% flux error.")

def add_hists_vals(hist_list):
    Total_hist = np.zeros_like(hist_list[0].values())
    for hist in hist_list:
        Total_hist += hist.values()
    return Total_hist

def add_all_errors(err_list): #adds in quadrature, assuming all hists are same shape
    Total_hist = np.zeros_like(err_list[0])
    for i in range(len(err_list[0])): #Looping over the bins
        for errs in err_list: #Looping over the histograms
            Total_hist[i] += errs[i]**2 #Adding error from each hist in quadrature
        Total_hist[i] = np.sqrt(Total_hist[i])
    return Total_hist

def add_all_errors_dict(err_dict): #adds in quadrature, assuming all hists are same shape
    list_keys = list(err_dict.keys())
    Total_hist = np.zeros_like(err_dict[list_keys[0]])
    for i in range(len(err_dict[list_keys[0]])): #Looping over the bins
        for errs in err_dict.keys(): #Looping over the histograms
            Total_hist[i] += err_dict[errs][i]**2 #Adding error from each hist in quadrature
        Total_hist[i] = np.sqrt(Total_hist[i])
    return Total_hist

def append_r3_to_r1_Owen(r1_list, r3_list): #Must be in a list form already, also removes two 0.0 values at the end of this hist
    r1_list.remove(0.0)
    r1_list.remove(0.0)
    r3_list.remove(0.0)
    r3_list.remove(0.0)
    TOTAL = r1_list + r3_list
    return TOTAL

def append_lists(hist_list): #hist_list should be a list of lists
    for hist in hist_list:
        Total += hist
    return Total

def get_full_errors_nu(hist, ppfx, gen):
    Total_err = np.zeros_like(hist.values())
    stat_err = hist.errors()
    det_sys_err = 0.3*hist.values()  #30% flat det sys err
    #det_sys_err = 0.7*hist.values()  #30% flat det sys err
    ppfx_err = []
    for i in range(len(ppfx.values())):
        ppfx_err.append(hist.values()[i]*ppfx.values()[i])
    gen_err = []
    for j in range(len(gen.values())):
        gen_err.append(hist.values()[j]*gen.values()[j])
    for k in range(len(Total_err)):
        Total_err[k] = np.sqrt(stat_err[k]**2 + det_sys_err[k]**2 + ppfx_err[k]**2 + gen_err[k]**2)
    return Total_err

def get_full_errors_nu_FLAT_INPUTS(hist):
    Total_err = np.zeros_like(hist.values())
    stat_err = hist.errors()
    det_sys_err = 0.3*hist.values()  #30% flat det sys err
    ppfx_err = 0.3*hist.values() #Took flat 30% ppfx
    gen_err = 0.2*hist.values() #Took flat 20% ppfx
    #det_sys_err = 0.7*hist.values()  #30% flat det sys err
    for k in range(len(Total_err)):
        Total_err[k] = np.sqrt(stat_err[k]**2 + det_sys_err[k]**2 + ppfx_err[k]**2 + gen_err[k]**2)
    return Total_err

def get_full_errors_signal(hist):
    Total_err = np.zeros_like(hist.values())
    stat_err = hist.errors()
    det_sys_err = 0.15*hist.values() #15% flat det sys err
    flux_err = 0.3*hist.values() #30% flat flux err
    #flux_err = 0.7*hist.values() #30% flat flux err
    for k in range(len(Total_err)):
        Total_err[k] = np.sqrt(stat_err[k]**2 + det_sys_err[k]**2 + flux_err[k]**2)
    return Total_err

def get_full_errors_dirt(hist):
    Total_err = np.zeros_like(hist.values())
    stat_err = hist.errors()
    dirt_unconstrained_err = 1.0*hist.values() #100% unconstrained err
    for k in range(len(Total_err)):
        Total_err[k] = np.sqrt(stat_err[k]**2 + dirt_unconstrained_err[k]**2)
    return Total_err

def get_stat_errors(hist):
    Total_err = np.zeros_like(hist.values())
    stat_err = hist.errors()
    return stat_err

def remove_first_half_hist(hist_list):
    length = len(hist_list)
    slice_at = int(np.floor(length/2))
    sliced_hist = hist_list[slice_at:]
    return sliced_hist

def make_zero_bin_unc(hist_dict, SF_dict, Params):
    """
    Given dict of hists and scaling factors.
    Returns the Poisson errors due to zero counts.
    For signal it assumes there is no scaling before this step.
    """
    zero_bins_errors = {}
    bkg_sample_names = ['bkg_overlay','bkg_EXT','bkg_dirt']
    conversion_dict = {'bkg_overlay':'overlay','bkg_EXT':'beamoff','bkg_dirt':'dirtoverlay'}
    if Params["Load_lepton_hists"] == True: 
        for HNL_mass in hist_dict:
            conversion_dict.update({HNL_mass:f"{HNL_mass}_ee"})
    if Params["Load_pi0_hists"] == True: 
        for HNL_mass in hist_dict:
            conversion_dict.update({HNL_mass:f"{HNL_mass}_pi0"})
    
    for HNL_mass in hist_dict:
        zero_bins_per_mass = {}
        for bkg in bkg_sample_names:
            zero_bins_per_mass[bkg] = np.zeros_like(hist_dict[HNL_mass][bkg].values())
            check_hist = hist_dict[HNL_mass][bkg].values()
            for i, val in enumerate(check_hist):
                if check_hist[i] == 0: 
                    one_sigma = 1.149 #This is the 1 sigma error on a poisson distribution for k=0
                    zero_bins_per_mass[bkg][i] = one_sigma*SF_dict[conversion_dict[bkg]]
                    
        zero_bins_per_mass["signal"] = np.zeros_like(hist_dict[HNL_mass]["signal"].values())
        check_hist = hist_dict[HNL_mass]["signal"].values()
        for i, val in enumerate(check_hist):
            if check_hist[i] == 0: 
                one_sigma = 1.149 #This is the 1 sigma error on a poisson distribution for k=0
                zero_bins_per_mass["signal"][i] = one_sigma*SF_dict[conversion_dict[HNL_mass]]
                
        zero_bins_errors[HNL_mass] = zero_bins_per_mass
        
    return zero_bins_errors

def Calculate_total_uncertainty(Params, hist_dict): #Takes the dictionary of all root files
    BKG_ERR_dict, SIGNAL_ERR_dict = {}, {}
    bkg_sample_names = ['bkg_overlay','bkg_EXT','bkg_dirt']
    overlay_sys_names = ["ppfx_uncertainty","Genie_uncertainty","Reinteraction_uncertainty","overlay_DetVar_uncertainty"]
    for HNL_mass in hist_dict:
        bkg_stat_err_dict, bkg_sys_err_dict = {}, {} #Clean for each mass point
        for name in bkg_sample_names:
            bkg_stat_err_dict[name]=hist_dict[HNL_mass][name].errors() #Load in stat error from error saved in hist
        sig_stat_err = hist_dict[HNL_mass]['signal'].errors()
        if Params["Stats_only"] == True: #Set all systematic errors to zero
            for name in bkg_sample_names:
                bkg_sys_err_dict[name] = np.zeros_like(hist_dict[HNL_mass][name].errors())
            sig_sys_err =  np.zeros_like(hist_dict[HNL_mass]['signal'].errors())
        elif Params["Use_flat_sys"] == True:
            for name in bkg_sample_names:
                bkg_sys_err_dict[name] = hist_dict[HNL_mass][name].values()*Params["Flat_"+name+"_frac"]
            sig_flux_err = hist_dict[HNL_mass]['signal'].values()*Params["Signal_flux_error"]
            sig_detvar_err = hist_dict[HNL_mass]['signal'].values()*Params["Flat_sig_detvar"]
            sig_sys_err = np.sqrt(sig_flux_err**2 + sig_detvar_err**2)
        elif Params["Use_flat_sys"] == False:
            overlay_sys_dict = {}
            for sys in overlay_sys_names:
                overlay_sys_dict[sys] = hist_dict[HNL_mass][sys].values()
            bkg_sys_err_dict['bkg_overlay'] = add_all_errors_dict(overlay_sys_dict)
            bkg_sys_err_dict['bkg_EXT'] = np.zeros_like(hist_dict[HNL_mass]['bkg_EXT'].errors())
            bkg_sys_err_dict['bkg_dirt'] = hist_dict[HNL_mass]['bkg_dirt'].values()*Params["Flat_bkg_dirt_frac"]
            
            sig_detvar_err = hist_dict[HNL_mass]["signal_DetVar_uncertainty"].values()
            sig_flux_err = hist_dict[HNL_mass]['signal'].values()*Params["Signal_flux_error"]
            sig_sys_err = add_all_errors([sig_detvar_err,sig_flux_err])
            
        #Evaluating final stat+sys errors    
        bkg_stat_plus_sys_dict={}
        for name in bkg_sample_names:
            bkg_stat_plus_sys_dict[name]=add_all_errors([bkg_stat_err_dict[name],bkg_sys_err_dict[name]]) #WRONG
        
        total_bkg_err = add_all_errors_dict(bkg_stat_plus_sys_dict) #Now adding the errors of overlay, EXT and dirt in quadrature
        total_sig_err = add_all_errors([sig_stat_err,sig_sys_err])
        
        BKG_ERR_dict[HNL_mass] = total_bkg_err
        SIGNAL_ERR_dict[HNL_mass] = total_sig_err
    return BKG_ERR_dict, SIGNAL_ERR_dict

def Uncertainty_breakdown(Params, hist_dict, bkg_reweight_err_dict=None, bkg_detvar_dict=None, sig_detvar_dict=None): #Takes the dictionary of all root files
    BKG_ERR_dict, SIGNAL_ERR_dict = {}, {}
    for HNL_mass in hist_dict:
        bkg_stat_err_list = [hist_dict[HNL_mass]['bkg_overlay'].errors(), 
                             hist_dict[HNL_mass]['bkg_EXT'].errors(), 
                             hist_dict[HNL_mass]['bkg_dirt'].errors()]
        sig_stat_err = hist_dict[HNL_mass]['signal'].errors()
        print("Signal stat error:")
        print(sig_stat_err)
        if Params["Stats_only"] == True:
        #As default the errors saved in the files are stat errors, this will change once I properly calculate them
            bkg_err_list = bkg_stat_err_list
            sig_err = sig_stat_err
        elif Params["Use_flat_sys"] == True:
            zero_bins = []
            for i,val in enumerate(hist_dict[HNL_mass]['bkg_overlay'].values()):
                if val == 0:
                    zero_bins.append(i)
                    print(f"{HNL_mass} last bin 0, setting error to 2.0")
            if len(zero_bins) != 0:
                bkg_sys_err_list = [hist_dict[HNL_mass]['bkg_overlay'].values()*Params["Flat_bkg_overlay_frac"] + np.ones_like(hist_dict[HNL_mass]['bkg_overlay'].values())*2.0, #This is horrible need to rewrite 
                                    np.zeros_like(hist_dict[HNL_mass]['bkg_EXT'].errors()), #No systematic error on the EXT sample
                                    hist_dict[HNL_mass]['bkg_dirt'].values()*Params["Flat_bkg_dirt_frac"]]
            else:    
                bkg_sys_err_list = [hist_dict[HNL_mass]['bkg_overlay'].values()*Params["Flat_bkg_overlay_frac"], 
                                    np.zeros_like(hist_dict[HNL_mass]['bkg_EXT'].errors()), #No systematic error on the EXT sample
                                    hist_dict[HNL_mass]['bkg_dirt'].values()*Params["Flat_bkg_dirt_frac"]]
            bkg_err_list = [add_all_errors([bkg_stat_err_list[0],bkg_sys_err_list[0]]), #adding the sys and stat error in quadrature for each bkg type
                            add_all_errors([bkg_stat_err_list[1],bkg_sys_err_list[1]]),
                            add_all_errors([bkg_stat_err_list[2],bkg_sys_err_list[2]])]
        elif Params["Use_flat_sys"] == False:
            ppfx_unc = hist_dict[HNL_mass]["ppfx_uncertainty"].values()
            genie_unc = hist_dict[HNL_mass]["Genie_uncertainty"].values()
            reint_unc = hist_dict[HNL_mass]["Reinteraction_uncertainty"].values()
            # detvar_unc = bkg_detvar_dict[HNL_mass]["Total_DetVar_uncertainty"].values() #Don't know what this looks like yet, as I haven't made
            detvar_unc = hist_dict[HNL_mass]['bkg_overlay'].values()*Params["Overlay_detvar_frac"] #Just setting as flat. Too much variation in samples
            tot_overlay_sys = add_all_errors([ppfx_unc, genie_unc, reint_unc, detvar_unc])
            bkg_sys_err_list = [tot_overlay_sys, 
                                np.zeros_like(hist_dict[HNL_mass]['bkg_EXT'].errors()), #No systematic error on the EXT sample
                                hist_dict[HNL_mass]['bkg_dirt'].values()*Params["Flat_bkg_dirt_frac"]] #Don't have reweight or DetVar samples for dirt
            bkg_err_list = [add_all_errors([bkg_stat_err_list[0],bkg_sys_err_list[0]]), #adding the sys and stat error in quadrature for each bkg type
                            add_all_errors([bkg_stat_err_list[1],bkg_sys_err_list[1]]),
                            add_all_errors([bkg_stat_err_list[2],bkg_sys_err_list[2]])]
            bkg_stat_err_total = add_all_errors([bkg_stat_err_list[0],bkg_stat_err_list[1],bkg_stat_err_list[2]])
            print("bkg stat error:")
            print(bkg_stat_err_total)
            print("bkg flux error:")
            print(ppfx_unc)
            print("bkg genie error:")
            print(genie_unc)
            print("bkg reint error:")
            print(reint_unc)
            
        if (Params["Stats_only"] == False) and (Params["Use_flat_sys"] == True):
            zero_bins = []
            for i,val in enumerate(hist_dict[HNL_mass]['signal'].values()):
                if val == 0:
                    zero_bins.append(i)
                    print(f"{HNL_mass} signal last bin 0, setting error to 2.0")
            if len(zero_bins) != 0:
                sig_sys_err = hist_dict[HNL_mass]['signal'].values()*Params["Signal_flux_error"]+2.0
            else:
                sig_sys_err = hist_dict[HNL_mass]['signal'].values()*Params["Signal_flux_error"]
            sig_err = add_all_errors([sig_stat_err,sig_sys_err])
        if (Params["Stats_only"] == False) and (Params["Use_flat_sys"] == False):
            sig_detvar_err = hist_dict[HNL_mass]['signal_DetVar_uncertainty'].values()
            sig_flux_err = hist_dict[HNL_mass]['signal'].values()*Params["Signal_flux_error"]
            sig_err = add_all_errors([sig_stat_err,sig_detvar_err,sig_flux_err]) #Adding stat, detvar and flux errors in quadrature
        total_bkg_err = add_all_errors(bkg_err_list) #Now adding the errors of overlay, EXT and dirt in quadrature
        BKG_ERR_dict[HNL_mass] = total_bkg_err
        SIGNAL_ERR_dict[HNL_mass] = sig_err
        print("Total bkg error:")
        print(total_bkg_err)
        print("Total signal error:")
        print(sig_err)
    return BKG_ERR_dict, SIGNAL_ERR_dict

def make_overflow_bin(bins_dict, bins_cents_dict):
    """
    For making the final "overflow" bin the same size as the previous bins, i.e one integer in width.
    """
    bins_overflow, bins_cent_overflow = {}, {}
    for HNL_mass in bins_dict:
        overflow_bin = bins_cents_dict[HNL_mass][-2]+1 #Just adding one to the penultimate bin centre val. 
        bins_cent_overflow[HNL_mass] = bins_cents_dict[HNL_mass].copy()
        bins_cent_overflow[HNL_mass][-1] = overflow_bin
        bins_overflow[HNL_mass] = bins_dict[HNL_mass].copy()
        bins_overflow[HNL_mass][-1] = bins_dict[HNL_mass][-2]+1 #Just adding one to the penultimate bin end val. 
    return bins_overflow, bins_cent_overflow

def make_xlims_dict(bins_dict, lower = None):
    """
    Making a dict of xlims for plotting several mass points at once.
    Also returns a dict of xticks for the purpose of indicating the overflow.
    """
    xlims_adjusted, xticks_adjusted = {}, {}
    for HNL_mass in bins_dict:
        if isinstance(lower,(int, float)): lower_val = lower
        else: lower_val = bins_dict[HNL_mass][0]
        xlims_adjusted[HNL_mass] = [lower_val,bins_dict[HNL_mass][-1]]
        ticks = np.arange(bins_dict[HNL_mass][0], bins_dict[HNL_mass][-1], 1)
        ticks_strings = []
        for val in ticks:
            ticks_strings.append(str(int(val)))
        ticks_strings[-1] = str(ticks_strings[-1])+"+"
        xticks_adjusted[HNL_mass] = ticks_strings
        
    return xlims_adjusted, xticks_adjusted


def remove_part_hist(hist_list, numbins):
        length = len(hist_list)
        slice_at = length - int(numbins)
        if slice_at < 0:
            print("Trying to use greater number of bins than available, using full dist.")
            return hist_list
        else:
            sliced_hist = hist_list[slice_at:]
            return sliced_hist

def Make_into_lists(Params, BKG_dict, SIGNAL_dict, TOT_ERR_dict):
    
    BKG_dict_FINAL, SIGNAL_dict_FINAL= {}, {}
    ERR_dict_FINAL = {}
    for HNL_mass in BKG_dict:
        ERR_list_dict = {}
        BKG = np.ndarray.tolist(BKG_dict[HNL_mass])
        SIGNAL = np.ndarray.tolist(SIGNAL_dict[HNL_mass])
        for err_dict in TOT_ERR_dict:
            ERR_list_dict[err_dict]=np.ndarray.tolist(TOT_ERR_dict[err_dict][HNL_mass])
        if Params["Use_part_only"] == True:
            numbins = Params["Num_bins_for_calc"] #Number of bins in signal region to use for CLs calc
            BKG=remove_part_hist(BKG, numbins)
            SIGNAL=remove_part_hist(SIGNAL, numbins)
            for err_dict in ERR_list_dict:
                ERR_list_dict[err_dict]=remove_part_hist(ERR_list_dict[err_dict], numbins)
            
        BKG_dict_FINAL[HNL_mass] = BKG
        SIGNAL_dict_FINAL[HNL_mass] = SIGNAL
        ERR_dict_FINAL[HNL_mass] = ERR_list_dict

    output_dict = {"BKG_dict":BKG_dict_FINAL, "SIGNAL_dict":SIGNAL_dict_FINAL}
    for err_dict in TOT_ERR_dict:
        new_err_dict_placeholder = {}
        for HNL_mass in BKG_dict:
            new_err_dict_placeholder[HNL_mass] = ERR_dict_FINAL[HNL_mass][err_dict]
        
        output_dict.update({err_dict:new_err_dict_placeholder})
        
    return output_dict

def Create_final_appended_runs_dict(list_input_dicts):
    
    def append_list_of_lists(input_list):
        output_list = []
        for i in range(len(input_list)):
            output_list = output_list + input_list[i]
        return output_list

    Total_dict = {}
    all_keys = list(list_input_dicts[0].keys())
    first_key = all_keys[0]
    for HNL_mass in list_input_dicts[0][first_key]:
    # for HNL_mass in Constants.HNL_mass_samples:
        Appended_dict = {}
        for dict_type in list_input_dicts[0].keys():
            list_placeholder = []
            for input_dict in list_input_dicts: #This loops over the dicts for different runs
                list_placeholder.append(input_dict[dict_type][HNL_mass]) 
            Appended = append_list_of_lists(list_placeholder)
            Appended_dict[dict_type] = Appended
        Total_dict[HNL_mass] = Appended_dict
    print(Total_dict.keys())
    return Total_dict
    
def check_duplicate_events(df):
    rse_list = df['rse_id'].to_list()

    seen = set()
    dupes = []

    for x in rse_list:
        if x in seen:
            dupes.append(x)
        else:
            seen.add(x)
    print("Number of duplicates is " + str(len(dupes)))
    print("Number of unique events is " + str(len(seen)))
    
def Load_pyhf_files(filenames, Params_pyhf, location='Uncertainties/', HNL_masses=Constants.HNL_mass_samples):
    loc_hists = location

    hist_dict_run1, hist_dict_run3, theta_dict = {}, {}, {}

    #Loading in the .root files
    if Params_pyhf["Load_lepton_hists"] == True:
        for HNL_mass in HNL_masses:
            hist_dict_run1[HNL_mass] = uproot.open(loc_hists+f'run1_{HNL_mass}MeV_' + filenames)
            if Params_pyhf["Load_single_r1_file"] == True: hist_dict_run3[HNL_mass] = hist_dict_run1[HNL_mass] 
            else:
                hist_dict_run3[HNL_mass] = uproot.open(loc_hists+f'run3_{HNL_mass}MeV_' + filenames)
            theta_dict[HNL_mass] = hist_dict_run1[HNL_mass]["theta"].values()[0] #assuming scaled theta is the same for all runs, only 1 value saved

    elif Params_pyhf["Load_pi0_hists"] == True:
        pi0_dict_run1, pi0_dict_run3 = {}, {}
        for HNL_mass in Constants.HNL_mass_pi0_samples:
            hist_dict_run1[HNL_mass] = uproot.open(loc_hists+f'pi0/run1_{HNL_mass}MeV_' + filenames)
            if Params_pyhf["Load_single_r1_file"] == True: hist_dict_run3[HNL_mass] = hist_dict_run1[HNL_mass]
            else:
                hist_dict_run3[HNL_mass] = uproot.open(loc_hists+f'pi0/run3_{HNL_mass}MeV_' + filenames)
            theta_dict[HNL_mass] = hist_dict_run1[HNL_mass]["theta"].values()[0]
            
    else:
        print("Loading Dirac samples")
        for HNL_mass in HNL_masses:
            hist_dict_run1[HNL_mass] = uproot.open(loc_hists+f'run1_{HNL_mass}MeV_' + filenames)
            if Params_pyhf["Load_single_r1_file"] == True: hist_dict_run3[HNL_mass] = hist_dict_run1[HNL_mass] 
            else:
                hist_dict_run3[HNL_mass] = uproot.open(loc_hists+f'run3_{HNL_mass}MeV_' + filenames)
            theta_dict[HNL_mass] = hist_dict_run1[HNL_mass]["theta"].values()[0]

    #list_of_dicts = [hist_dict_run1, hist_dict_run3] #Add runs 2, 4 when available, not using yet

    theta_squared = Constants.theta_mu_4*Constants.theta_mu_4

    all_hists_list = ['bkg_overlay;1', 'bkg_dirt;1', 'bkg_EXT;1', 'signal;1', 'data;1', 'theta;1', 
                      'ppfx_uncertainty;1', 'Genie_uncertainty;1', 'Reinteraction_uncertainty;1', 
                      'ppfx_uncertainty_frac;1', 'Genie_uncertainty_frac;1', 'Reinteraction_uncertainty_frac;1', 
                      'overlay_DetVar_uncertainty;1', 'overlay_DetVar_uncertainty_frac;1', 'signal_DetVar_uncertainty;1', 'signal_DetVar_uncertainty_frac;1']
    missing_hists = []
    for hist_name in all_hists_list:
        if hist_name not in hist_dict_run1[HNL_mass].keys(): missing_hists.append(hist_name)
    if len(missing_hists) == 0: print("No missing histograms in Run1")
    else:
        print("Missing hists for Run1 are: ")
        print(missing_hists)
    print("thetas are:")
    print(theta_dict)
    print("Done")
    
    return hist_dict_run1, hist_dict_run3, theta_dict

def New_Load_pyhf_files(filenames, Params_pyhf, location='Uncertainties/', HNL_masses=Constants.HNL_mass_samples):
    loc_hists = location

    hist_dict_run1, hist_dict_run3, theta_dict = {}, {}, {}

    #Loading in the .root files
    if Params_pyhf["Load_lepton_hists"] == True:
        for HNL_mass in HNL_masses:
            hist_dict_run1[HNL_mass] = uproot.open(loc_hists+f'run1_{HNL_mass}_' + filenames)
            if Params_pyhf["Load_single_r1_file"] == True: hist_dict_run3[HNL_mass] = hist_dict_run1[HNL_mass] 
            else:
                hist_dict_run3[HNL_mass] = uproot.open(loc_hists+f'run3_{HNL_mass}_' + filenames)
            theta_dict[HNL_mass] = hist_dict_run1[HNL_mass]["theta"].values()[0] #assuming scaled theta is the same for all runs, only 1 value saved

    elif Params_pyhf["Load_pi0_hists"] == True:
        pi0_dict_run1, pi0_dict_run3 = {}, {}
        for HNL_mass in Constants.HNL_mass_pi0_samples:
            hist_dict_run1[HNL_mass] = uproot.open(loc_hists+f'pi0/run1_{HNL_mass}_' + filenames)
            if Params_pyhf["Load_single_r1_file"] == True: hist_dict_run3[HNL_mass] = hist_dict_run1[HNL_mass]
            else:
                hist_dict_run3[HNL_mass] = uproot.open(loc_hists+f'pi0/run3_{HNL_mass}_' + filenames)
            theta_dict[HNL_mass] = hist_dict_run1[HNL_mass]["theta"].values()[0]
            
    else:
        print("Loading Dirac samples")
        for HNL_mass in HNL_masses:
            hist_dict_run1[HNL_mass] = uproot.open(loc_hists+f'run1_{HNL_mass}_' + filenames)
            if Params_pyhf["Load_single_r1_file"] == True: hist_dict_run3[HNL_mass] = hist_dict_run1[HNL_mass] 
            else:
                hist_dict_run3[HNL_mass] = uproot.open(loc_hists+f'run3_{HNL_mass}_' + filenames)
            theta_dict[HNL_mass] = hist_dict_run1[HNL_mass]["theta"].values()[0]

    #list_of_dicts = [hist_dict_run1, hist_dict_run3] #Add runs 2, 4 when available, not using yet

    theta_squared = Constants.theta_mu_4*Constants.theta_mu_4

    all_hists_list = ['bkg_overlay;1', 'bkg_dirt;1', 'bkg_EXT;1', 'signal;1', 'data;1', 'theta;1', 
                      'ppfx_uncertainty;1', 'Genie_uncertainty;1', 'Reinteraction_uncertainty;1', 
                      'ppfx_uncertainty_frac;1', 'Genie_uncertainty_frac;1', 'Reinteraction_uncertainty_frac;1', 
                      'overlay_DetVar_uncertainty;1', 'overlay_DetVar_uncertainty_frac;1', 'signal_DetVar_uncertainty;1', 'signal_DetVar_uncertainty_frac;1']
    missing_hists = []
    for hist_name in all_hists_list:
        if hist_name not in hist_dict_run1[HNL_mass].keys(): missing_hists.append(hist_name)
    if len(missing_hists) == 0: print("No missing histograms in Run1")
    else:
        print("Missing hists for Run1 are: ")
        print(missing_hists)
    print("thetas are:")
    print(theta_dict)
    print("Done")
    
    return hist_dict_run1, hist_dict_run3, theta_dict

#Limit plot
def Pandafy(path):
    cols = ['Mass','Value']
    df = pd.read_csv(path,names=cols)
    firstLine = pd.DataFrame([[df['Mass'][0],1.]],columns=cols)
    lastLine = pd.DataFrame([[df['Mass'][-1:].values[0],1.]],columns=cols)
    df = pd.concat([firstLine,df])
    df = pd.concat([df,lastLine])
    return df

def Pandafy_new(path):
    cols = ['Mass','Value']
    df = pd.read_csv(path,names=cols)
    return df