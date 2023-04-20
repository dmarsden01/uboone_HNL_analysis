import matplotlib.pyplot as plt
import awkward as awk
import matplotlib as mpl
import numpy as np
import math
import pandas as pd
import uproot
import uproot3
import functools

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
        

def create_sample_list(Params): #Returns an extended parameter dict and a the list of samples to run over
    if Params["FLATTEN"] == True: Params["Flat_state"] = "flattened"
    else: Params["Flat_state"] = "unflattened"
    if Params["Only_keep_common_DetVar_evs"] == True: Params["Reduced_state"] = "reduced_evs"
    else: Params["Reduced_state"] = "all_evs"

    if Params["only_presel"]:
        Params["variables_string"] = "Presel_vars"
        Params["variables"] = Variables.Preselection_vars_CRT + Variables.event_vars
        Params["variables_MC"] = Variables.Preselection_vars_CRT_MC + Variables.event_vars
    elif Params["Load_truth_vars"]:
        Params["variables_string"] = "Truth_vars"
        Params["variables"] = Variables.event_vars + Variables.event_vars
        Params["variables_MC"] = Variables.event_vars + Variables.event_vars
    else:
        Params["variables_string"] = "my_vars"
        Params["variables"] = Variables.Final_variable_list
        Params["variables_MC"] = Variables.Final_variable_list_MC

    if Params["Run"] == "run1": Params["current"] = "FHC"
    elif Params["Run"] == "run3": Params["current"] = "RHC"
    else: print("Need to choose either \"run1\" or \"run3\"")

    samples = [] #A list of all the samples which will be loaded and pickled
    if Params["Load_lepton_signal"] == True: samples.extend(["signal"])
    if Params["Load_standard_bkgs"] == True: samples.extend(["overlay","dirtoverlay","beamoff"])
    if Params["Load_pi0_signal"] == True: samples.extend(["pi0_signal"])
    if Params["Load_DetVars"] == True: samples.extend(Constants.Detector_variations)
    if Params["Load_Signal_DetVars"] == True:
        # for HNL_mass in Constants.HNL_mass_samples: #For when all detvar samples are made
        if Params["Run"] == "run1":
            for HNL_mass in [50, 100, 150]:
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_"+DetVar]
        if Params["Run"] == "run3":
            for HNL_mass in [2, 10, 20, 50, 100]: #Don't have 150MeV sample yet
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_"+DetVar]
    if Params["Load_pi0_signal_DetVars"] == True:
        if Params["Run"] == "run1": print("There are no pi0 detvar samples for run1.")
        if Params["Run"] == "run3":
            for HNL_mass in [150, 180, 220, 240,245]: #200 is broken for some reason
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_"+DetVar]
    if Params["Load_data"] == True: samples.extend(["beamgood"])
    if Params["Load_single_file"] == True: samples = [Params["single_file"]]
        
    print(f"Loading these "+Params["Run"]+" samples: " + "\n" + str(samples))
    
    return Params, samples

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
    
def make_common_evs_df(df_list):
    overlapping_df = functools.reduce(lambda left,right: pd.merge(left, right, on=['rse_id'], how='inner'), df_list)
    print("Length is of common events list is " + str(len(overlapping_df)))
    return overlapping_df

# def make_common_indices_df(df_list):
#     overlapping_df = functools.reduce(lambda left,right: pd.merge(left, right, on=['rse_id'], how='inner'), df_list)
#     print("Length is of common events list is " + str(len(overlapping_df)))
#     return overlapping_df

def Load_and_pkl_samples(samples, sample_loc, loc_pkls, common_evs, Params, save_str=""):
    for sample in samples: #Looping over all samples, should make a function for this.
        if sample in Constants.Detector_variations: #Checks if it is an overlay Detector Variation sample
            print(f"Loading overlay {sample} "+Params["Run"]+" with uproot")
            NuMI_MC_overlay=uproot3.open(f"../NuMI_MC/DetVars/neutrinoselection_filt_"+Params["Run"]+f"_overlay_{sample}.root")[Constants.root_dir+"/"+Constants.main_tree]
            df_overlay = NuMI_MC_overlay.pandas.df(Params["variables_MC"], flatten=Params["FLATTEN"])
            file = df_overlay
            new_file = file.copy()
            del(file)
            if Params["Run"]=="run3":
                print("skipping fiducal vars")
                final_file = new_file.copy() # Not making fiducial variables because of broken thetaXZ file
            else:
                final_file = Make_fiducial_vars(new_file)
            del(new_file)
            Edit_Weight_Tune(final_file)
            MC_weight_branch(final_file)
            make_unique_ev_id(final_file) #This creates "rse_id" branch
            if Params["Only_keep_common_DetVar_evs"] == True:
                filtered = final_file.loc[(final_file['rse_id'].isin(common_evs['rse_id']))]
                new_overlay = filtered.copy()
                del(final_file)
                del(filtered)
            else:
                new_overlay = final_file.copy()
                del(final_file)
            print("Pickling "+Params["Run"]+f" overlay {sample} file")
            new_overlay.to_pickle(loc_pkls+"DetVars/overlay_"+Params["Run"]+"_"+Params["variables_string"]+f"_{sample}_"+Params["Flat_state"]+"_"+Params["Reduced_state"]+save_str+".pkl")
            del(new_overlay)
        elif Params["Load_Signal_DetVars"] == True: #I should ONLY load these samples in this case.
            first_str = sample.split("_")[0]
            HNL_mass = int(first_str)
            NuMI_MC_signal=uproot3.open("../NuMI_signal/KDAR_dump/sfnues/DetVars/"+f"{sample}_"+Params["Run"]+".root")[Constants.root_dir+"/"+Constants.main_tree]
            df_signal = NuMI_MC_signal.pandas.df(Params["variables"], flatten=Params["FLATTEN"])
            file = df_signal
            new_file = file.copy()
            del(file)
            final_file = Make_fiducial_vars(new_file)
            del(new_file)
            final_file = make_unique_ev_id(final_file) #This creates "rse_id" branch
            if Params["Only_keep_common_DetVar_evs"] == True:
                filtered = final_file.loc[(final_file['rse_id'].isin(common_evs[HNL_mass]['rse_id']))]
                new_overlay = filtered.copy()
                del(final_file)
                del(filtered)
            else:
                new_overlay = final_file.copy()
                del(final_file)
            new_overlay.to_pickle(loc_pkls+"Signal_DetVars/"+Params["Run"]+f"_{sample}_"+Params["Reduced_state"]+save_str+".pkl")
            del(new_overlay)
                
        elif Params["Load_pi0_signal_DetVars"] == True: #I should ONLY load these samples in this case.
            first_str = sample.split("_")[0]
            HNL_mass = int(first_str)
            NuMI_MC_signal=uproot3.open("../NuMI_signal/KDAR_dump/sfnues/pi0/DetVars/"+f"{sample}_"+Params["Run"]+".root")[Constants.root_dir+"/"+Constants.main_tree]
            df_signal = NuMI_MC_signal.pandas.df(Params["variables"], flatten=Params["FLATTEN"])
            file = df_signal
            new_file = file.copy()
            del(file)
            final_file = Make_fiducial_vars(new_file)
            del(new_file)
            final_file = make_unique_ev_id(final_file) #This creates "rse_id" branch
            if Params["Only_keep_common_DetVar_evs"] == True:
                filtered = final_file.loc[(final_file['rse_id'].isin(common_evs[HNL_mass]['rse_id']))]
                new_overlay = filtered.copy()
                del(final_file)
                del(filtered)
            else:
                new_overlay = final_file.copy()
                del(final_file)
            new_overlay.to_pickle(loc_pkls+"Signal_DetVars/pi0/"+Params["Run"]+f"_{sample}_"+Params["Reduced_state"]+save_str+".pkl")
            del(new_overlay)
            
        else: #Standard sample types
            print(f"Loading {sample} "+Params["Run"]+" file(s) with uproot")
            # if Constants.sample_type[sample] == "MC_signal": #Meant to be only e+e- signal samples
            if sample == "signal": #Meant to be only e+e- signal samples
                for HNL_mass in Constants.HNL_mass_samples:
                    file_loc = sample_loc[sample]+f"{HNL_mass}_ee_Umu4_majorana_"+Params["current"]+".root"
                    uproot_file = uproot3.open(file_loc)[Constants.root_dir+"/"+Constants.main_tree]
                    if Params["Load_truth_vars"] == True:
                        df_signal = uproot_file.pandas.df(Params["variables"] + Variables.Truth_vars, flatten=Params["FLATTEN"])
                    else:
                        df_signal = uproot_file.pandas.df(Params["variables"], flatten=Params["FLATTEN"])
                    file = df_signal
                    new_signal = file.copy()
                    del(file)
                    if (Params["Load_truth_vars"] == False) and (Params["FLATTEN"] == True):
                        final_signal = Make_fiducial_vars(new_signal) #New since last save
                    else:
                        final_signal = new_signal.copy()
                    del(new_signal)
                    print("Pickling "+Params["Run"]+ f" {HNL_mass}MeV file")
                    final_signal.to_pickle(loc_pkls+f"signal_{HNL_mass}MeV_"+Params["Run"]+"_"+Params["variables_string"]+"_"+Params["Flat_state"]+save_str+".pkl")
                    del(final_signal)
            elif sample == "pi0_signal":
                for HNL_mass in Constants.HNL_mass_pi0_samples:
                    file_loc = sample_loc[sample]+f"{HNL_mass}_pi0_Umu4_majorana_"+Params["current"]+".root"
                    uproot_file = uproot3.open(file_loc)[Constants.root_dir+"/"+Constants.main_tree]
                    if Params["Load_truth_vars"] == True:
                        df_signal = uproot_file.pandas.df(Params["variables"] + Variables.Truth_vars, flatten=Params["FLATTEN"])
                    else:
                        df_signal = uproot_file.pandas.df(Params["variables"], flatten=Params["FLATTEN"])
                    file = df_signal
                    new_signal = file.copy()
                    del(file)
                    if Params["Load_truth_vars"] == False:
                        final_signal = Make_fiducial_vars(new_signal)
                    if Params["Load_truth_vars"] == True:
                        final_signal = new_signal.copy()
                    del(new_signal)
                    print("Pickling "+Params["Run"]+ f" pi0 {HNL_mass}MeV file")
                    final_signal.to_pickle(loc_pkls+f"pi0_signal_{HNL_mass}MeV_"+Params["Run"]+"_"+Params["variables_string"]+"_"+Params["Flat_state"]+save_str+".pkl")
                    del(final_signal)
                    
            elif (Params["Load_single_file"] == True) and (isinstance(sample,int)):
                HNL_mass = sample
                file_loc = sample_loc["signal"]+f"{HNL_mass}_Umu4_majorana_numi_"+Params["current"]+".root"
                uproot_file = uproot3.open(file_loc)[Constants.root_dir+"/"+Constants.main_tree]
                if Params["Load_truth_vars"] == True:
                    df_signal = uproot_file.pandas.df(Params["variables"] + Variables.Truth_vars, flatten=Params["FLATTEN"])
                else:
                    df_signal = uproot_file.pandas.df(Params["variables"], flatten=Params["FLATTEN"])
                file = df_signal
                new_signal = file.copy()
                del(file)
                final_signal = Make_fiducial_vars(new_signal)
                del(new_signal)
                print("Pickling "+Params["Run"]+ f" {HNL_mass}MeV file")
                final_signal.to_pickle(loc_pkls+f"signal_{HNL_mass}MeV_"+Params["Run"]+"_"+Params["variables_string"]+"_"+Params["Flat_state"]+save_str+".pkl")
                del(final_signal)
                    
            else:
                uproot_file = uproot3.open(sample_loc[sample])[Constants.root_dir+'/'+Constants.main_tree]
                if Constants.sample_type[sample] == "MC":
                    df = uproot_file.pandas.df(Params["variables_MC"], flatten=Params["FLATTEN"])
                    file = df
                    Edit_Weight_Tune(file)
                    MC_weight_branch(file)
                else:
                    df = uproot_file.pandas.df(Params["variables"], flatten=Params["FLATTEN"])
                    file = df
                new_file = file.copy()
                del(file)
                if (Params["Load_truth_vars"] == False) and (Params["FLATTEN"] == True):
                    final_file = Make_fiducial_vars(new_file)
                else:
                    final_file = new_file.copy()
                del(new_file)
                print("Pickling "+Params["Run"] +f" {sample} file")
                final_file.to_pickle(loc_pkls+f"{sample}_"+Params["Run"]+"_"+Params["variables_string"]+"_"+Params["Flat_state"]+save_str+".pkl")
                del(final_file)

#POT counting
def POT_counter(file): #Takes uproot file
    Total_POT = file["pot"].array().sum()
    return Total_POT

#Preselection
def Preselection_weighted_efficiency(samples, cut_dict, Efficiency_dict, Preselected): #Need to account for weigthing in overlay and dirt samples
    for sample in samples:
        if sample == "overlay" or sample == "dirtoverlay":
            weight = samples[sample]["weight"]
            NumEvs = sum(weight)
        else:
            NumEvs = len(samples[sample])
        
        effic_list = [1.0]
        for cut in cut_dict.keys():
            samples[sample]=samples[sample].query(cut_dict[cut])
            if sample == "overlay" or sample == "dirtoverlay":
                weight = samples[sample]["weight"]
                Num_selected = sum(weight)
            else:
                Num_selected = len(samples[sample])
            effic_list.append(Num_selected/NumEvs)
        Efficiency_dict[sample]=effic_list
        #samples.update()
        Selected = samples[sample].copy()
        placeholder_dict = {sample:Selected}
        Preselected.update(placeholder_dict) 

#Doing BDT training
def Prepare_dfs_for_xgb(df): #The default value for missing data in XGB is 0. So this changes those very large negative values to -9999.
    value = -1e15
    new_value = -9999
    first_entry = df.index[0]
    for variable in df.keys():
        if isinstance(df[variable][first_entry], (int,float,np.int32,np.float32,np.uint32)):
        # if variable=='rse_id': #Should come up with a better way of checking the "type" of variable, in case it is not int or float.
        #     continue           #But don't know how to access the first extant row of a dataframe (since some have been removed). 
        # else: 
            if(len(df.loc[df[variable] < value]) > 0):
                df.loc[(df[variable] < value), variable] = new_value #Sets the new value
            if(len(df.loc[df[variable] == -1.0]) > 0):
                df.loc[(df[variable] == -1.0), variable] = new_value #Sets the new value
            if(len(df.loc[df[variable] == np.nan]) > 0):
                df.loc[(df[variable] == np.nan), variable] = new_value #Sets the new value
            if(len(df.loc[df[variable] == np.inf]) > 0):
                df.loc[(df[variable] == np.inf), variable] = new_value #Sets the new value
        # else:
        #     print(variable)
            
    df_edited = df.copy() 
    return df_edited

def only_keep_highest_E(df):
    df["highest_E"]=df['pfnplanehits_Y'].groupby("entry").transform(max) == df['pfnplanehits_Y']
    df_new = df.query("highest_E").copy()
    return df_new

def create_test_samples_list(Params): #Returns the list of samples to run over
    samples = [] #A list of all the samples which will be loaded and pickled
    if Params["Load_standard"] == True:
        samples.extend(["overlay","dirtoverlay","beamoff"])
    if Params["Load_lepton_signal"] == True:
        samples.extend(Constants.HNL_mass_samples)
    if Params["Load_pi0_signal"] == True:
        for HNL_mass in Constants.HNL_mass_pi0_samples:
            samples+=[str(HNL_mass)+"_pi0"]
    if Params["Load_DetVars"] == True: #This is overlay DetVars
        samples.extend(Constants.Detector_variations)
    if Params["Load_Signal_DetVars"] == True:
        # for HNL_mass in Constants.HNL_mass_samples: #For when all detvar samples are made
        for HNL_mass in [50, 100, 150, 180, 200]:
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
                    samples+=[str(HNL_mass)+"_"+DetVar]
        if Params["Run"] == "run3":
            for HNL_mass in [2, 10, 20, 50, 100]: #Don't have 150MeV sample yet
                masses+=[HNL_mass]
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_"+DetVar]
    if Params["Load_pi0_signal"] == True:
        if Params["Run"] == "run1":
            print("There are no run1 pi0 detector variation samples")
        if Params["Run"] == "run3":
            for HNL_mass in [150,180,220,240,245]: #200 broken for some reason
                masses+=[HNL_mass]
                for DetVar in Constants.Detector_variations:
                    samples+=[str(HNL_mass)+"_"+DetVar]
        
    print(f"Loading these "+Params["Run"]+" samples: " + "\n")
    print(samples)
    
    if Params["Use_logit"] == True:
        Params["logit_str"] = "logit"
    else: Params["logit_str"] = "standard"
    
    return samples, masses

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
        print("Using fully evaluated systematic uncertainty for background. Dirt will still be 100%.")
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

def make_zero_bin_unc(hist_dict, SF_dict):
    zero_bins_errors = {}
    bkg_sample_names = ['bkg_overlay','bkg_EXT','bkg_dirt']
    conversion_dict = {'bkg_overlay':'overlay','bkg_EXT':'beamoff','bkg_dirt':'dirtoverlay'}
    for HNL_mass in hist_dict:
        zero_bins_per_mass = {}
        for bkg in bkg_sample_names:
            zero_bins_per_mass[bkg] = np.zeros_like(hist_dict[HNL_mass][bkg].values())
            check_hist = hist_dict[HNL_mass][bkg].values()
            for i, val in enumerate(check_hist):
                if check_hist[i] == 0: 
                    one_sigma = 1.149 #This is the 1 sigma error on a poisson distribution for k=0
                    zero_bins_per_mass[bkg][i] = one_sigma*SF_dict[conversion_dict[bkg]]
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

def Calculate_total_uncertainty_OLD(Params, hist_dict, bkg_reweight_err_dict=None, bkg_detvar_dict=None, sig_detvar_dict=None): #Takes the dictionary of all root files
    BKG_ERR_dict, SIGNAL_ERR_dict = {}, {}
    for HNL_mass in Constants.HNL_mass_samples:
        bkg_stat_err_list = [hist_dict[HNL_mass]['bkg_overlay'].errors(), 
                             hist_dict[HNL_mass]['bkg_EXT'].errors(), 
                             hist_dict[HNL_mass]['bkg_dirt'].errors()]
        sig_stat_err = hist_dict[HNL_mass]['Signal'].errors()
        if Params["Stats_only"] == True:
        #As default the errors saved in the files are stat errors, this will change once I properly calculate them
            bkg_err_list = bkg_stat_err_list
            sig_err = sig_stat_err
        elif Params["Use_flat_sys_bkg"] == True:
            bkg_sys_err_list = [hist_dict[HNL_mass]['bkg_overlay'].values()*Params["Flat_overlay_bkg_frac"], 
                                np.zeros_like(hist_dict[HNL_mass]['bkg_EXT'].errors()), #No systematic error on the EXT sample
                                hist_dict[HNL_mass]['bkg_dirt'].values()*Params["Flat_dirt_bkg_frac"]]
            bkg_err_list = [add_all_errors([bkg_stat_err_list[0],bkg_sys_err_list[0]]), #adding the sys and stat error in quadrature for each bkg type
                            add_all_errors([bkg_stat_err_list[1],bkg_sys_err_list[1]]),
                            add_all_errors([bkg_stat_err_list[2],bkg_sys_err_list[2]])]
        elif Params["Use_flat_sys_bkg"] == False:
            ppfx_unc = bkg_reweight_err_dict[HNL_mass]["ppfx_uncertainty"].values()
            genie_unc = bkg_reweight_err_dict[HNL_mass]["Genie_uncertainty"].values()
            reint_unc = bkg_reweight_err_dict[HNL_mass]["Reinteraction_uncertainty"].values()
            detvar_unc = bkg_detvar_dict[HNL_mass]["Total_DetVar_uncertainty"].values() #Don't know what this looks like yet, as I haven't made
            tot_overlay_sys = add_all_errors([ppfx_unc, genie_unc, reint_unc, detvar_unc])
            bkg_sys_err_list = [tot_overlay_sys, 
                                0, #No systematic error on the EXT sample
                                hist_dict[HNL_mass]['bkg_dirt'].values()*Params["Flat_dirt_bkg_frac"]] #Don't have reweight or DetVar samples for dirt
            bkg_err_list = [add_all_errors([bkg_stat_err_list[0],bkg_sys_err_list[0]]), #adding the sys and stat error in quadrature for each bkg type
                            add_all_errors([bkg_stat_err_list[1],bkg_sys_err_list[1]]),
                            add_all_errors([bkg_stat_err_list[2],bkg_sys_err_list[2]])]
        if (Params["Stats_only"] == False) and (Params["Use_flat_sys_signal"] == True):
            sig_sys_err = hist_dict[HNL_mass]['Signal'].values()*Params["Flat_sig_frac"]
            sig_err = add_all_errors([sig_stat_err,sig_sys_err])
        if (Params["Stats_only"] == False) and (Params["Use_flat_sys_signal"] == False):
            sig_detvar_err = sig_detvar_dict[HNL_mass]["Total_DetVar_uncertainty"].values()
            sig_flux_err = hist_dict[HNL_mass]['Signal'].values()*Params["Flat_overlay_bkg_frac"]
            sig_err = add_all_errors([sig_stat_err,sig_detvar_err,sig_flux_err]) #Adding stat, detvar and flux errors in quadrature
        total_bkg_err = add_all_errors(bkg_err_list) #Now adding the errors of overlay, EXT and dirt in quadrature
        BKG_ERR_dict[HNL_mass] = total_bkg_err
        SIGNAL_ERR_dict[HNL_mass] = sig_err
    return BKG_ERR_dict, SIGNAL_ERR_dict

def Make_into_lists_OLD(Params, BKG_dict, SIGNAL_dict, BKG_ERR_dict, SIGNAL_ERR_dict):
    
    def remove_part_hist(hist_list, numbins):
        length = len(hist_list)
        slice_at = length - int(numbins)
        if slice_at < 0:
            print("Trying to use greater number of bins than available, using full dist.")
            slice_at = 0
        sliced_hist = hist_list[slice_at:]
        return sliced_hist
    
    BKG_dict_FINAL, BKG_ERR_dict_FINAL, SIGNAL_dict_FINAL, SIGNAL_ERR_dict_FINAL = {}, {}, {}, {}
    for HNL_mass in BKG_dict:
        BKG = np.ndarray.tolist(BKG_dict[HNL_mass])
        BKG_ERR = np.ndarray.tolist(BKG_ERR_dict[HNL_mass])
        SIGNAL = np.ndarray.tolist(SIGNAL_dict[HNL_mass])
        SIGNAL_ERR = np.ndarray.tolist(SIGNAL_ERR_dict[HNL_mass])
        if Params["Use_part_only"] == True:
            numbins = Params["Num_bins_for_calc"] #Number of bins in signal region to use for CLs calc
            BKG=remove_part_hist(BKG, numbins)
            BKG_ERR=remove_part_hist(BKG_ERR, numbins)
            SIGNAL=remove_part_hist(SIGNAL, numbins)
            SIGNAL_ERR=remove_part_hist(SIGNAL_ERR, numbins)
            
        BKG_dict_FINAL[HNL_mass] = BKG
        BKG_ERR_dict_FINAL[HNL_mass] = BKG_ERR
        SIGNAL_dict_FINAL[HNL_mass] = SIGNAL
        SIGNAL_ERR_dict_FINAL[HNL_mass] = SIGNAL_ERR
        
    output_dict = {"BKG_dict":BKG_dict_FINAL, "BKG_ERR_dict":BKG_ERR_dict_FINAL, 
                   "SIGNAL_dict":SIGNAL_dict_FINAL, "SIGNAL_ERR_dict":SIGNAL_ERR_dict_FINAL}
        
    return output_dict

def Make_into_lists(Params, BKG_dict, SIGNAL_dict, TOT_ERR_dict):
    
    def remove_part_hist(hist_list, numbins):
        length = len(hist_list)
        slice_at = length - int(numbins)
        if slice_at < 0:
            print("Trying to use greater number of bins than available, using full dist.")
            slice_at = 0
        sliced_hist = hist_list[slice_at:]
        return sliced_hist
    
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

    if Params_pyhf["Load_pi0_hists"] == True:
        pi0_dict_run1, pi0_dict_run3 = {}, {}
        for HNL_mass in Constants.HNL_mass_pi0_samples:
            hist_dict_run1[HNL_mass] = uproot.open(loc_hists+f'pi0/run1_{HNL_mass}MeV_' + filenames)
            if Params_pyhf["Load_single_r1_file"] == True: hist_dict_run3[HNL_mass] = hist_dict_run1[HNL_mass]
            else:
                hist_dict_run3[HNL_mass] = uproot.open(loc_hists+f'pi0/run3_{HNL_mass}MeV_' + filenames)
            theta_dict[HNL_mass] = hist_dict_run1[HNL_mass]["theta"].values()[0]

    #list_of_dicts = [hist_dict_run1, hist_dict_run3] #Add run2 when available, not using yet

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