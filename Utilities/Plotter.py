import matplotlib.pyplot as plt
import awkward as awk
import matplotlib as mpl
import numpy as np
import math
import Utilities.SystematicsCalc as SC
import Utilities.Consts as C

import Utilities.Constants as Constants
import Utilities.Functions as Functions

plt.rcParams.update({'font.size': 22})

params = {'legend.fontsize': 15,
          'figure.figsize': (15, 5),
         'axes.labelsize': 'large',
         'axes.titlesize':'medium',
         'xtick.labelsize':'small',
         'ytick.labelsize':'small'}
plt.rcParams.update(params)



from math import floor, log10
# Define function for string formatting of scientific notation
def sci_notation(num, decimal_digits=0, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"{0:.{1}f}\cdot10^{{{{".format(coeff,precision)+"{0:d}".format(exponent)+"}}$"


def Plot_preselection_variable(variable, samples=[], sample_norms=[], xlabel=[],xlims=[0,0],bins=40,figsize=[10,10],dpi=100,MergeBins=False, discrete=False, HNL_mass = 0, HNLplotscale=100000,density=False,legloc="best",logy = "False", cutline = 0.0):
    
    if(samples==[]): raise Exception("Specify samples dict") 
    if(xlabel==[]): xlabel=variable
    
    #beamgood=samples["beamgood"] #I should loop through samples instead, so don't always need data
    beamoff=samples["beamoff"]
    overlay=samples["overlay"]
    dirtoverlay=samples["dirtoverlay"]
    signal=samples["signal"]
    
    #var_Data=beamgood[var]
    var_Offbeam=beamoff[variable]
    var_Overlay=overlay[variable]
    var_Dirt=dirtoverlay[variable]
    var_HNL=signal[variable]
    
    weight_Offbeam=np.ones(len(var_Offbeam))*sample_norms["beamoff"]
    weight_Overlay=overlay["weight"]*sample_norms["overlay"]
    weight_Dirt=dirtoverlay["weight"]*sample_norms["dirtoverlay"]
    weight_signal=np.ones(len(var_HNL))*sample_norms["signal"]*HNLplotscale
    
    if xlims[0] == 0 and xlims[1] == 0: xlims = [min(var_Overlay),max(var_Overlay)]
    
    if(isinstance(bins, int)):
        nbins=bins
        bins=np.linspace(xlims[0],xlims[1],nbins+1)
    else: nbins=len(bins)-1
    
    if(MergeBins): #remove bins with zero bkg prediction
        totbkg=np.histogram(var_Offbeam,bins=bins,range=xlims)[0]+np.histogram(var_Dirt,bins=bins,range=xlims)[0]+np.histogram(var_Overlay,bins=bins,range=xlims)[0]
        offbkg=np.histogram(var_Offbeam,bins=bins,range=xlims)[0]
        overlaybkg=np.histogram(var_Overlay,bins=bins,range=xlims)[0]
        totbkg=np.histogram(var_Offbeam,bins=bins,range=xlims)[0]+np.histogram(var_Dirt,bins=bins,range=xlims)[0]+np.histogram(var_Overlay,bins=bins,range=xlims)[0]
        bins_new=[]
        for i,bin_bkg in enumerate(totbkg):
            if(offbkg[i]>1 or overlaybkg[i]>1):
                bins_new.append(bins[i])
                
        bins_new.append(bins[-1])

        bins=bins_new
    
    fig,ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=figsize,dpi=dpi)
        
    if(discrete):
        bins = np.arange(xlims[0], xlims[1] + 1.5) - 0.5
        xlims[0]=xlims[0]-1
        xlims[1]=xlims[1]+1
        ax[0].set_xticks(bins + 0.5)
        nbins=len(bins)-1

    Offbeamnum=sum(weight_Offbeam)
    Dirtnum=sum(weight_Dirt)
    Overlaynum=sum(weight_Overlay) 
    HNL_num=sum(weight_signal)

    varis=[var_Overlay,var_Dirt,var_Offbeam]
    weights=[weight_Overlay,weight_Dirt,weight_Offbeam]
    colors=['peru',"darkorange",'deepskyblue']
    labels=[fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})",f"Beam-Off ({Offbeamnum:.1f})"]
        
    plot=plt.hist(varis,
              label=labels,
              range=xlims,bins=bins,
              histtype="stepfilled",
              stacked=True,density=density,linewidth=2,edgecolor="black",
              weights=weights, color=colors)
    

    color="red"
    
    bkg_stack=varis
    bkg_stack_w=weights
    plt.hist(var_HNL,
              label=[f"{HNL_mass} MeV HNL ({HNL_num:.1f})"],
              range=xlims,bins=bins,
              stacked=True,density=density,
              weights= weight_signal,histtype="step",color=color,lw=4)
    
    if cutline != 0.0:
        plt.axvline(x=cutline, lw=3, color='green', linestyle = 'dashed')
    
    if(logy == "True"):
        plt.yscale("log")
    else:
        plt.yscale("linear")
    
    plt.legend(loc=legloc,frameon=False)
    
    plt.xlabel(xlabel)
    plt.xlim(xlims)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
          

def HNL_scaling_calculator(samples=[], sample_norms=[]): #Prints the value which should be used to scale HNL hist to overlay hist
    
    beamoff=samples["beamoff"]
    overlay=samples["overlay"]
    dirtoverlay=samples["dirtoverlay"]
    signal=samples["signal"]
    
    weight_Offbeam=np.ones(len(beamoff["nslice"]))*sample_norms["beamoff"]
    weight_Overlay=overlay["weight"]*sample_norms["overlay"]
    weight_Dirt=dirtoverlay["weight"]*sample_norms["dirtoverlay"]
    weight_signal=np.ones(len(signal["nslice"]))*sample_norms["signal"]
    
    Offbeamnum=sum(weight_Offbeam)
    Dirtnum=sum(weight_Dirt)
    Overlaynum=sum(weight_Overlay) 
    HNL_num=sum(weight_signal)
    
    ratio_overlay_HNL = Overlaynum/HNL_num
    
    ratio_all_bkg_HNL = (Offbeamnum+Dirtnum+Overlaynum)/HNL_num
    
    print("The ratio of overlay to HNL events is " + "%.0f" % ratio_overlay_HNL + "\n")
    
    print("The ratio of all bkgs to HNL events is " + "%.0f" % ratio_all_bkg_HNL + "\n")
    

# def Presel_efficiency():

# def Plot_BDT_input():

def Plot_BDT_output(HNL_masses=[], samples=[], sample_norms=[], colours={}, ALPHA=1.0, xlims=[0,1.0],bins=20,figsize=[12,8], MergeBins=False, density=False, legloc="upper center",logy=True, savefig=False, save_str="", Run="_", logit=False, HNL_scale=1.0):
    
    if(HNL_masses==[]): raise Exception("Specify HNL sample masses")
    if(samples==[]): raise Exception("Specify samples")
    if(colours=={}): colours = {'overlay_test':Constants.sample_colours['overlay'],
                                'dirtoverlay':Constants.sample_colours['dirtoverlay'],
                                'beamoff':Constants.sample_colours['beamoff'],
                                'signal_test':Constants.sample_colours['signal']}
    
    if logy == True:
        logscale="log"
    elif logy == False:
        logscale="linear"
    
    for HNL_mass in HNL_masses:
        plt.figure(figsize=figsize,facecolor='white')
        if logit == False:
            bkg_scores=[samples['overlay_test'][f'BDT_output_{HNL_mass}MeV'],samples['dirtoverlay'][f'BDT_output_{HNL_mass}MeV'],
                   samples['beamoff'][f'BDT_output_{HNL_mass}MeV']]
        if logit == True:
            bkg_scores=[Functions.logit(samples['overlay_test'][f'BDT_output_{HNL_mass}MeV']),
                        Functions.logit(samples['dirtoverlay'][f'BDT_output_{HNL_mass}MeV']),
                        Functions.logit(samples['beamoff'][f'BDT_output_{HNL_mass}MeV'])]
        bkg_weights=[sample_norms['overlay_test'],sample_norms['dirtoverlay'],sample_norms['beamoff']]
        bkg_colors=[colours['overlay_test'],colours['dirtoverlay'],colours['beamoff']]
        labels=[fr"In-Cryo $\nu$",fr"Out-Cryo $\nu$",f"Beam-Off"]
        
        bins_list = np.histogram(bkg_scores[0],bins=bins,range=xlims)[1] #For mergebins part
              
        if(MergeBins): #remove bins with zero bkg prediction
            totbkg=np.histogram(bkg_scores[0],bins=bins,range=xlims)[0]+np.histogram(bkg_scores[1],bins=bins,range=xlims)[0]+np.histogram(bkg_scores[2],bins=bins,range=xlims)[0]
            offbkg=np.histogram(bkg_scores[2],bins=bins,range=xlims)[0]
            overlaybkg=np.histogram(bkg_scores[0],bins=bins,range=xlims)[0]
            dirtbkg=np.histogram(bkg_scores[1],bins=bins,range=xlims)[0]
            bins_new=[]
            for i,bin_bkg in enumerate(totbkg):
                if(offbkg[i]>1 or overlaybkg[i]>1):
                    bins_new.append(bins_list[i])

            bins_new.append(bins_list[-1])

            bins=bins_new
        
        plot=plt.hist(bkg_scores,
              label=labels,
              range=xlims,bins=bins,
              histtype="stepfilled",
              stacked=True,linewidth=2,edgecolor="black",
              weights=bkg_weights, color=bkg_colors, alpha=ALPHA)
        if logit == False:
            plt.hist(samples[HNL_mass][f'BDT_output_{HNL_mass}MeV'],weights=sample_norms[HNL_mass]*HNL_scale,bins=bins,range=xlims,
                     lw=4, edgecolor=colours['signal_test'], label=f'HNL {HNL_mass} MeV', histtype="step")
        if logit == True:
            plt.hist(Functions.logit(samples[HNL_mass][f'BDT_output_{HNL_mass}MeV']),
                     weights=sample_norms[HNL_mass]*HNL_scale,bins=bins,range=xlims,
                     lw=4, edgecolor=colours['signal_test'], label=f'HNL {HNL_mass} MeV', histtype="step")
        plt.legend(loc=legloc,frameon=True)
        
        plt.xlabel('BDT score', fontsize=30)
        plt.ylabel('Events', fontsize=30)
        plt.rcParams.update({'font.size': 30})
        plt.yscale(logscale)
        if savefig == True:
            plt.savefig("plots/BDT_output/" + Run + "_" + str(HNL_mass) + "_MeV_" + logscale + save_str + ".png")
        plt.show()
        
def Plot_BDT_output_systematics(HNL_masses=[], samples=[], colours={}, ALPHA=1.0, xlims=[0,1.0],figsize=[12,8], density=False, legloc="upper center",logy=True, savefig=False, save_str="", Run="_", logit=False, HNL_scale=1.0):
    """
    This should take the histograms which have already been binned and scaled and plot the total uncertainties on bkg.
    Therefore it will display what is being fed into the limit setting software.
    """
    
    if(HNL_masses==[]): raise Exception("Specify HNL sample masses")
    if(samples==[]): raise Exception("Specify samples")
    if(colours=={}): colours = {'overlay_test':Constants.sample_colours['overlay'],
                                'dirtoverlay':Constants.sample_colours['dirtoverlay'],
                                'beamoff':Constants.sample_colours['beamoff'],
                                'signal_test':Constants.sample_colours['signal']}
    
    if logy == True:
        logscale="log"
    elif logy == False:
        logscale="linear"
    
    for HNL_mass in HNL_masses:
        plt.figure(figsize=figsize,facecolor='white')
        if logit == False:
            bkg_scores=[samples['overlay_test'][f'BDT_output_{HNL_mass}MeV'],samples['dirtoverlay'][f'BDT_output_{HNL_mass}MeV'],
                   samples['beamoff'][f'BDT_output_{HNL_mass}MeV']]
        if logit == True:
            bkg_scores=[Functions.logit(samples['overlay_test'][f'BDT_output_{HNL_mass}MeV']),
                        Functions.logit(samples['dirtoverlay'][f'BDT_output_{HNL_mass}MeV']),
                        Functions.logit(samples['beamoff'][f'BDT_output_{HNL_mass}MeV'])]
        bkg_weights=[sample_norms['overlay_test'],sample_norms['dirtoverlay'],sample_norms['beamoff']]
        bkg_colors=[colours['overlay_test'],colours['dirtoverlay'],colours['beamoff']]
        labels=[fr"In-Cryo $\nu$",fr"Out-Cryo $\nu$",f"Beam-Off"]
        
        bins_list = np.histogram(bkg_scores[0],bins=bins,range=xlims)[1] #For mergebins part
        tot_uncrt = 
              
        plot=plt.hist(bkg_scores,
              label=labels,
              range=xlims,bins=bins,
              histtype="stepfilled",
              stacked=True,linewidth=2,edgecolor="black",
              weights=bkg_weights, color=bkg_colors, alpha=ALPHA)
        if logit == False:
            plt.hist(samples[HNL_mass][f'BDT_output_{HNL_mass}MeV'],weights=sample_norms[HNL_mass]*HNL_scale,bins=bins,range=xlims,
                     lw=4, edgecolor=colours['signal_test'], label=f'HNL {HNL_mass} MeV', histtype="step")
        if logit == True:
            plt.hist(Functions.logit(samples[HNL_mass][f'BDT_output_{HNL_mass}MeV']),
                     weights=sample_norms[HNL_mass]*HNL_scale,bins=bins,range=xlims,
                     lw=4, edgecolor=colours['signal_test'], label=f'HNL {HNL_mass} MeV', histtype="step")
        plt.legend(loc=legloc,frameon=True)
        
        plt.xlabel('BDT score', fontsize=30)
        plt.ylabel('Events', fontsize=30)
        plt.rcParams.update({'font.size': 30})
        plt.yscale(logscale)
        if savefig == True:
            plt.savefig("plots/BDT_output/" + Run + "_" + str(HNL_mass) + "_MeV_" + logscale + save_str + ".png")
        plt.show()
        
                    
# def Plot_systematics():
    
# def Plot_Final_Limit():