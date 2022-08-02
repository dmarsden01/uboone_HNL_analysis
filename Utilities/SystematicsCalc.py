
import pandas as pd
import root_numpy as rn
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import awkward as awk
import matplotlib as mpl
import ROOT
import os
import pickle as pl
from tqdm import tqdm_notebook as tqdm
import uproot
# NuMI needs to add PPFX workaround for dirt

whichdf="vert_df"

def sys_err(sample,name, var_name, query, x_range, bins,NormScale):
    # how many universes?
    Nuniverse = 600 #100 #len(df)
#     print("Universes",Nuniverse)

    if(isinstance(bins, int)):n_bins=bins
    else: n_bins=len(bins)-1

    
    n_tot = np.empty([Nuniverse, n_bins])
    n_cv_tot = np.empty(n_bins)
    n_tot.fill(0)
    n_cv_tot.fill(0)

        # for pi0 fit only
        #if ((t in ["ncpi0","ccpi0"]) and (name == "weightsGenie") ):
        #    continue

    queried_sample = sample.query(query)
    variable = queried_sample[var_name]
    syst_weights = queried_sample[name]
    #print ('N universes is :',len(syst_weights))
    spline_fix_cv  = queried_sample["weight"]*NormScale
    spline_fix_var = queried_sample["weight"]*NormScale
#     if (name == "weightsGenie"): 
#         spline_fix_var = queried_sample["ppfx_cv"]*NormScale

    s = syst_weights
    df = pd.DataFrame(s.values.tolist())
    #print (df)
    #continue
    n_cv, bins = np.histogram(
        variable,
        range=x_range,
        bins=bins,
        weights=spline_fix_cv)
    n_cv_tot += n_cv
    
    if(name == "weightsGenie"): #special treatment as ["weightSplineTimesTune"] is included in genie weights
        if not df.empty:
            for i in range(Nuniverse):
                weight = df[i].values / 1000.
                weight[weight == 1]= queried_sample["weightSplineTimesTune"].iloc[weight == 1]
                weight[np.isnan(weight)] = queried_sample["weightSplineTimesTune"].iloc[np.isnan(weight)]
                weight[weight > 30] = queried_sample["weightSplineTimesTune"].iloc[weight > 30] #build up of events at 65 weight 
                weight[weight < 0] = queried_sample["weightSplineTimesTune"].iloc[weight < 0] 
                weight[weight == np.inf] = queried_sample["weightSplineTimesTune"].iloc[weight == np.inf]
                
                n, bins = np.histogram(
                    variable, weights=weight*np.nan_to_num(spline_fix_var/queried_sample["weightSplineTimesTune"]), range=x_range, bins=bins)
                n_tot[i] += n
                
    elif(name == "weightsPPFX"): #special treatment as ["PPFXPcv"] is included in ppfx weights
        if not df.empty:
            for i in range(Nuniverse):
                weight = df[i].values / 1000.
                weight[weight == 1]= queried_sample["ppfx_cv"].iloc[weight == 1]
                weight[np.isnan(weight)] = queried_sample["ppfx_cv"].iloc[np.isnan(weight)]
                weight[weight > 100] = queried_sample["ppfx_cv"].iloc[weight > 100]
                weight[weight < 0] = queried_sample["ppfx_cv"].iloc[weight < 0]
                weight[weight == np.inf] = queried_sample["ppfx_cv"].iloc[weight == np.inf]
                
                n, bins = np.histogram(
                    variable, weights=weight*np.nan_to_num(spline_fix_var/queried_sample["ppfx_cv"]), range=x_range, bins=bins)
                n_tot[i] += n
    else:       
        if not df.empty:
            for i in range(Nuniverse):
                weight = df[i].values / 1000.
                weight[np.isnan(weight)] = 1
                weight[weight > 100] = 1
                weight[weight < 0] = 1
                weight[weight == np.inf] = 1
                n, bins = np.histogram(
                    variable, weights=weight*spline_fix_var, range=x_range, bins=bins)
                n_tot[i] += n

    cov = np.empty([len(n_cv), len(n_cv)])
    cov.fill(0)

#     print()
#     bin1s=[]
    for n in n_tot:
#         print(n/2)
#         bin1s.append(n[2])
        for i in range(len(n_cv)):
            for j in range(len(n_cv)):
                cov[i][j] += (n[i] - n_cv_tot[i]) * (n[j] - n_cv_tot[j])


    cov /= Nuniverse
#     print(np.sqrt(np.diag(cov))/n_cv)
    return cov,n_cv_tot,n_tot,bins


def sys_err_old(sample,name, var_name, query, x_range, n_bins,NormScale):
    # how many universes?
    Nuniverse = 600 #100 #len(df)
#     print("Universes",Nuniverse)

    n_tot = np.empty([Nuniverse, n_bins])
    n_cv_tot = np.empty(n_bins)
    n_tot.fill(0)
    n_cv_tot.fill(0)

        # for pi0 fit only
        #if ((t in ["ncpi0","ccpi0"]) and (name == "weightsGenie") ):
        #    continue

    queried_sample = sample.query(query)
    variable = queried_sample[var_name]
    syst_weights = queried_sample[name]
    #print ('N universes is :',len(syst_weights))
    spline_fix_cv  = queried_sample["weight"]*NormScale
    spline_fix_var = queried_sample["weight"]*NormScale
#     if (name == "weightsGenie"): 
#         spline_fix_var = queried_sample["ppfx_cv"]*NormScale

    s = syst_weights
    df = pd.DataFrame(s.values.tolist())
    #print (df)
    #continue
    n_cv, bins = np.histogram(
        variable,
        range=x_range,
        bins=n_bins,
        weights=spline_fix_cv)
    n_cv_tot += n_cv
    
    if(name == "weightsGenie"): #special treatment as ["weightSplineTimesTune"] is included in genie weights
        if not df.empty:
            for i in range(Nuniverse):
                weight = df[i].values / 1000.
                weight[weight == 1]= queried_sample["weightSplineTimesTune"].iloc[weight == 1]
                weight[np.isnan(weight)] = queried_sample["weightSplineTimesTune"].iloc[np.isnan(weight)]
                weight[weight > 30] = queried_sample["weightSplineTimesTune"].iloc[weight > 30] #build up of events at 65 weight 
                weight[weight < 0] = queried_sample["weightSplineTimesTune"].iloc[weight < 0] 
                weight[weight == np.inf] = queried_sample["weightSplineTimesTune"].iloc[weight == np.inf]
                
                n, bins = np.histogram(
                    variable, weights=weight*np.nan_to_num(spline_fix_var/queried_sample["weightSplineTimesTune"]), range=x_range, bins=n_bins)
                n_tot[i] += n
                
    elif(name == "weightsPPFX"): #special treatment as ["PPFXPcv"] is included in ppfx weights
        if not df.empty:
            for i in range(Nuniverse):
                weight = df[i].values / 1000.
                weight[weight == 1]= queried_sample["ppfx_cv"].iloc[weight == 1]
                weight[np.isnan(weight)] = queried_sample["ppfx_cv"].iloc[np.isnan(weight)]
                weight[weight > 100] = queried_sample["ppfx_cv"].iloc[weight > 100]
                weight[weight < 0] = queried_sample["ppfx_cv"].iloc[weight < 0]
                weight[weight == np.inf] = queried_sample["ppfx_cv"].iloc[weight == np.inf]
                
                n, bins = np.histogram(
                    variable, weights=weight*np.nan_to_num(spline_fix_var/queried_sample["ppfx_cv"]), range=x_range, bins=n_bins)
                n_tot[i] += n
    else:       
        if not df.empty:
            for i in range(Nuniverse):
                weight = df[i].values / 1000.
                weight[np.isnan(weight)] = 1
                weight[weight > 100] = 1
                weight[weight < 0] = 1
                weight[weight == np.inf] = 1
                n, bins = np.histogram(
                    variable, weights=weight*spline_fix_var, range=x_range, bins=n_bins)
                n_tot[i] += n

    cov = np.empty([len(n_cv), len(n_cv)])
    cov.fill(0)

#     print()
#     bin1s=[]
    for n in n_tot:
#         print(n/2)
#         bin1s.append(n[2])
        for i in range(len(n_cv)):
            for j in range(len(n_cv)):
                cov[i][j] += (n[i] - n_cv_tot[i]) * (n[j] - n_cv_tot[j])


    cov /= Nuniverse
#     print(np.sqrt(np.diag(cov))/n_cv)
    return cov,n_cv_tot,n_tot,bins





def det_sys_err(sample,name, var_name, query, x_range, n_bins,NormScale):
    # how many universes?
    Nuniverse = 600 #100 #len(df)
#     print("Universes",Nuniverse)

    n_tot = np.empty([Nuniverse, n_bins])
    n_cv_tot = np.empty(n_bins)
    n_tot.fill(0)
    n_cv_tot.fill(0)

        # for pi0 fit only
        #if ((t in ["ncpi0","ccpi0"]) and (name == "weightsGenie") ):
        #    continue

    queried_sample = sample.query(query)
    variable = queried_sample[var_name]
    syst_weights = queried_sample[name]
    #print ('N universes is :',len(syst_weights))
    spline_fix_cv  = queried_sample["weight"]*NormScale
    spline_fix_var = queried_sample["weight"]*NormScale
#     if (name == "weightsGenie"): 
#         spline_fix_var = queried_sample["ppfx_cv"]*NormScale

    s = syst_weights
    df = pd.DataFrame(s.values.tolist())
    #print (df)
    #continue
    n_cv, bins = np.histogram(
        variable,
        range=x_range,
        bins=n_bins,
        weights=spline_fix_cv)
    n_cv_tot += n_cv
    
    if(name == "weightsGenie"): #special treatment as ["weightSplineTimesTune"] is included in genie weights
        if not df.empty:
            for i in range(Nuniverse):
                weight = df[i].values / 1000.
                weight[weight == 1]= queried_sample["weightSplineTimesTune"].iloc[weight == 1]
                weight[np.isnan(weight)] = queried_sample["weightSplineTimesTune"].iloc[np.isnan(weight)]
                weight[weight > 30] = queried_sample["weightSplineTimesTune"].iloc[weight > 30]
                weight[weight < 0] = queried_sample["weightSplineTimesTune"].iloc[weight < 0] 
                weight[weight == np.inf] = queried_sample["weightSplineTimesTune"].iloc[weight == np.inf]
                
                n, bins = np.histogram(
                    variable, weights=weight*np.nan_to_num(spline_fix_var/queried_sample["weightSplineTimesTune"]), range=x_range, bins=n_bins)
                n_tot[i] += n
    else:       
        if not df.empty:
            for i in range(Nuniverse):
                weight = df[i].values / 1000.
                weight[np.isnan(weight)] = 1
                weight[weight > 100] = 1
                weight[weight < 0] = 1
                weight[weight == np.inf] = 1
                n, bins = np.histogram(
                    variable, weights=weight*spline_fix_var, range=x_range, bins=n_bins)
                n_tot[i] += n

    cov = np.empty([len(n_cv), len(n_cv)])
    cov.fill(0)

#     print()
    bin1s=[]
    for n in n_tot:
#         print(n/2)
        bin1s.append(n[2])
        for i in range(len(n_cv)):
            for j in range(len(n_cv)):
                cov[i][j] += (n[i] - n_cv_tot[i]) * (n[j] - n_cv_tot[j])


    cov /= Nuniverse
#     print(np.sqrt(np.diag(cov))/n_cv)
    return cov,n_cv_tot,n_tot,bins