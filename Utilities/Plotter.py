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


def Plot_preselection_variable(variable, samples=[], sample_norms=[], xlabel=[],xlims=[0,0],bins=40,figsize=[10,10],dpi=100,MergeBins=False, discrete=False, HNL_mass = 0, HNLplotscale=100000,density=False,legloc="best",logy = "False", cutline = 0.0, show_ev_nums=True):
    
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
    if show_ev_nums==True:
        labels=[fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})",f"Beam-Off ({Offbeamnum:.1f})"]
        sig_label = [f"{HNL_mass} MeV HNL ({HNL_num:.1f})"]
    else:
        labels=[fr"In-Cryo $\nu$",fr"Out-Cryo $\nu$",f"Beam-Off"]
        sig_label = [f"{HNL_mass} MeV HNL"]
        
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
    
def Plot_preselection_variable_data(variable, samples=[], sample_norms=[], xlabel=[],xlims=[0,0],bins=40,figsize=[10,10],dpi=100,MergeBins=False, discrete=False, HNL_mass = 0, HNLplotscale=100000,density=False,legloc="best",logy = False, cutline = 0.0, show_ev_nums=False, CalcSys=False, xticks=[], colours_sample={}, order=[], sys_dict={}, centre_bins=False, hatch=False, ylabel="Events", Frame=True, arrow_place=[], ylimit=None, legsize=22, display=True, savefig=False, savename="test"):
    
    if(samples==[]): raise Exception("Specify samples dict") 
    if(xlabel==[]): xlabel=variable
    if(colours_sample=={}): colours_sample = {'overlay':Constants.sample_colours['overlay'],
                                              'dirtoverlay':Constants.sample_colours['dirtoverlay'],
                                              'beamoff':Constants.sample_colours['beamoff'],
                                              'signal':Constants.sample_colours['signal']}
    if(order==[]): order = ["beamoff","overlay","dirtoverlay"] #From bottom to top in stack
    if(sys_dict=={} and CalcSys==True): raise Exception("Specify systematic errors dict")
    
    beamgood=samples["beamgood"] #I should loop through samples instead, so don't always need data
    beamoff=samples["beamoff"]
    overlay=samples["overlay"]
    dirtoverlay=samples["dirtoverlay"]
    signal=samples["signal"]
    
    var_Data=beamgood[variable]
    var_Offbeam=beamoff[variable]
    var_Overlay=overlay[variable]
    var_Dirt=dirtoverlay[variable]
    var_HNL=signal[variable]
    
    variable_sample = {'overlay':var_Overlay,
                       'dirtoverlay':var_Dirt,
                       'beamoff':var_Offbeam,
                       'signal':var_HNL}
    
    weight_Offbeam=np.ones(len(var_Offbeam))*sample_norms["beamoff"]
    weight_Overlay=overlay["weight"]*sample_norms["overlay"]
    weight_Dirt=dirtoverlay["weight"]*sample_norms["dirtoverlay"]
    weight_signal=np.ones(len(var_HNL))*sample_norms["signal"]*HNLplotscale
    
    weights_sample = {'overlay':weight_Overlay,
                      'dirtoverlay':weight_Dirt,
                      'beamoff':weight_Offbeam,
                      'signal':weight_signal}
    
    if xlims[0] == 0 and xlims[1] == 0: xlims = [min(var_Overlay),max(var_Overlay)]
    
    if(isinstance(bins, int)):
        nbins=bins
        if centre_bins == True:
            bins=np.linspace(xlims[0],xlims[1],nbins+1)-0.5
        else:
            bins=np.linspace(xlims[0],xlims[1],nbins+1)
    else: nbins=len(bins)-1
    
    #all UNWEIGHTED hists
    totbkg=np.histogram(var_Offbeam,bins=bins,range=xlims)[0]+np.histogram(var_Dirt,bins=bins,range=xlims)[0]+np.histogram(var_Overlay,bins=bins,range=xlims)[0]
    offbkg=np.histogram(var_Offbeam,bins=bins,range=xlims)[0]
    overlaybkg=np.histogram(var_Overlay,bins=bins,range=xlims)[0]
    dirtbkg=np.histogram(var_Dirt,bins=bins,range=xlims)[0]
    
    #weighted hists
    offbkg_weighted=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam)[0]
    overlaybkg_weighted=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay)[0]
    dirtbkg_weighted=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt)[0]
           
    #Testing Owens way, Err = sqrt(N*S.F**2)
    mc_w=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay**2)
    off_w=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam**2)
    dirt_w=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt**2)
    
    off_err=np.sqrt(off_w[0])
    mc_err=np.sqrt(mc_w[0])
    dirt_err=np.sqrt(dirt_w[0])
    
    # stat_bkgerr=np.sqrt(offbkg_stat**2+overlaybkg_stat**2+dirtbkg_stat**2) #Adding stat errors in quadrature, my way
    stat_bkgerr=np.sqrt(off_err**2+mc_err**2+dirt_err**2) #Adding stat errors in quadrature
    tot_mcerr = stat_bkgerr
    
    if(CalcSys): #Owen's way of calculating systematics for arbitrary variable, using my new function
        frac_total = 0
        for frac_sys in sys_dict[variable]:
            frac_total += frac_sys**2 #Adding systematic sources of bkg in quadrature
        
        total_frac_sys = np.sqrt(frac_total)
        total_sys_err = frac_total*overlaybkg_weighted
        
        dirt_norm_err_fac = 1.0
        dirt_norm_err=dirtbkg_weighted*dirt_norm_err_fac
        tot_mcerr=np.sqrt( stat_bkgerr**2+total_sys_err**2+dirt_norm_err**2)
    
    if(MergeBins): #remove bins with zero bkg prediction
        bins_new=[]
        for i,bin_bkg in enumerate(totbkg):
            if(offbkg[i]>1 or overlaybkg[i]>1):
                bins_new.append(bins[i])
                
        bins_new.append(bins[-1])

        bins=bins_new
    
    # fig,ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=figsize,dpi=dpi) #Just variable plot
    fig,ax = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=figsize,dpi=dpi)
    
    plt.sca(ax[0])
        
    if(discrete):
        bins = np.arange(xlims[0], xlims[1] + 1.5) - 0.5
        xlims[0]=xlims[0]-1
        xlims[1]=xlims[1]+1
        ax[0].set_xticks(bins + 0.5)
        nbins=len(bins)-1
        
    x,y=np.histogram(var_Data,bins=bins,range=xlims,density=density)
    x1,y=np.histogram(var_Data,bins=bins,range=xlims)
    bin_center = [(y[i] + y[i+1])/2. for i in range(len(y)-1)]
    dat_val=x
    # dat_err=np.sqrt(x1)*Functions.safe_div(x,x1) #need to write one for arrays instead of single values.
    dat_err=np.sqrt(x1)*np.nan_to_num(x/x1)

    Datanum=dat_val.sum()
    Offbeamnum=sum(weight_Offbeam)
    Dirtnum=sum(weight_Dirt)
    Overlaynum=sum(weight_Overlay) 
    HNL_num=sum(weight_signal)
    
    if show_ev_nums==True:
        labels_sample = {'overlay':fr"In-Cryo $\nu$ ({Overlaynum:.1f})",
                         'dirtoverlay':fr"Out-Cryo $\nu$ ({Dirtnum:.1f})",
                         'beamoff':f"Beam-Off ({Offbeamnum:.1f})",
                         'signal':f"{HNL_mass} MeV HNL ({HNL_num:.1f})"}
        # labels=[fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})",f"Beam-Off ({Offbeamnum:.1f})"]
        sig_label = [f"{HNL_mass} MeV HNL ({HNL_num:.1f})"]
        data_label = f"NuMI Data ({Datanum:.0f})"
    else:
        labels_sample = {'overlay':fr"In-Cryo $\nu$",
                         'dirtoverlay':fr"Out-Cryo $\nu$",
                         'beamoff':f"Beam-Off",
                         'signal':f"{HNL_mass} MeV HNL"}
        # labels=[fr"In-Cryo $\nu$",fr"Out-Cryo $\nu$",f"Beam-Off"]
        sig_label = [f"{HNL_mass} MeV HNL"]
        data_label = "NuMI Data"
    
    plt.errorbar(bin_center,dat_val,yerr=dat_err,fmt='.',color='black',lw=5,capsize=5,elinewidth=3,label=data_label) #Plotting data

    varis, weights, colors, labels = [], [], [], []
    for sample in order:
        varis.append(variable_sample[sample])
        weights.append(weights_sample[sample])
        colors.append(colours_sample[sample])
        labels.append(labels_sample[sample])
            
    plot=plt.hist(varis,
              label=labels,
              range=xlims,bins=bins,
              histtype="stepfilled",
              stacked=True,density=density,linewidth=2,edgecolor="black",
              weights=weights, color=colors)
    
    upvals=np.append((plot[0][2]+tot_mcerr),(plot[0][2]+tot_mcerr)[-1])
    lowvals=np.append((plot[0][2]-tot_mcerr),(plot[0][2]-tot_mcerr)[-1])
    
    if hatch == False:
        plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)
    if hatch == True:
        plt.fill_between(y, lowvals, upvals,step="post",hatch='//',alpha=0,zorder=2)

    color=colours_sample["signal"]
    
    bkg_stack=varis
    bkg_stack_w=weights
    plt.hist(var_HNL,
              label=sig_label,
              range=xlims,bins=bins,
              stacked=True,density=density,
              weights= weight_signal,histtype="step",color=color,lw=4)
    
    if cutline != 0.0:
        plt.axvline(x=cutline, lw=3, color='green', linestyle = 'dashed')
    if arrow_place != []: #No idea how to fix this.
        plt.arrow(arrow_place[0], arrow_place[1], arrow_place[2], arrow_place[3], color='green', shape='full',fill=False, 
                  length_includes_head=True, overhang=0.0, head_width=-0.1, head_length=0.15, lw=3)
    
    if(logy == True):
        plt.yscale("log")
    else:
        plt.yscale("linear")
        
        
    plt.ylabel(ylabel)
    plt.legend(loc=legloc,frameon=Frame, prop={'size': legsize})
    
    # plt.xlabel(xlabel)
    plt.xlim(xlims)
    if ylimit != None: plt.ylim(0,ylimit)
    
    #---Sub-plot----#
    plt.sca(ax[1])
    
    fracer_data=np.nan_to_num(np.sqrt(x1)/x1)
    x_err=fracer_data*x
    fracer_mc=np.nan_to_num(tot_mcerr/plot[0][2])
    
    rat_err_data=x_err*(1/plot[0][2])
    
    rat_err_mc=fracer_mc
    rat_err=np.sqrt(rat_err_data**2)

    rat_err_mc=np.nan_to_num(rat_err_mc) #other wise the next doesnt plot pro[erly]

    upvals= np.append(1+(rat_err_mc),1+(rat_err_mc)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
    lowvals=np.append(1-(rat_err_mc),1-(rat_err_mc)[-1])


    plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)

    rat=np.nan_to_num(x/plot[0][2])
    rat[x==0]=1 #dont think this is a good way to deal with this

    rat_err=np.nan_to_num(rat*np.sqrt(fracer_mc**2+fracer_data**2))
       
    # plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data") #Had this before, but wrong I think
    plt.errorbar(bin_center,rat,yerr=fracer_data,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
    plt.ylabel("Data/MC")
    plt.axhline(1,ls='-',color='black')
    plt.axhline(1.1,ls='--',color='grey')
    plt.axhline(0.9,ls='--',color='grey')
    ylim = max(abs(np.nan_to_num(rat)))*1.1
    plt.ylim(0.7,1.3)
    # plt.ylim(0.9,1.1)
    plt.xlim(xlims)
    if xticks != []:
        plt.xticks(xticks)
    
    plt.xlabel(xlabel)
    
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    
    if savefig==True:
        plt.savefig(savename+".pdf")
        plt.savefig(savename+".png")
    if display == False:
        plt.close()
          

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
                if(overlaybkg[i]>1):
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
            plt.savefig("plots/BDT_output/BDT_testing/" + Run + "_" + str(HNL_mass) + "MeV_" + logscale + save_str + ".pdf")
            plt.savefig("plots/BDT_output/BDT_testing/" + Run + "_" + str(HNL_mass) + "MeV_" + logscale + save_str + ".png")
        plt.show()
        
def Plot_BDT_output_data(HNL_masses=[], samples=[], sample_norms=[], colours={}, ALPHA=1.0, xlims=[0,1.0],bins=20,figsize=[12,8], MergeBins=False, density=False, legloc="upper center",logy=True, savefig=False, save_str="", Run="_", logit=False, HNL_scale=1.0, dpi=100):
    
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
        
        fig,ax = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=figsize,dpi=dpi)
    
        #----primary plot----#
        plt.sca(ax[0])
        
        # plt.figure(figsize=figsize,facecolor='white')
        if logit == False:
            bkg_scores=[samples['overlay_test'][f'BDT_output_{HNL_mass}MeV'],samples['dirtoverlay'][f'BDT_output_{HNL_mass}MeV'],
                   samples['beamoff'][f'BDT_output_{HNL_mass}MeV']]
            dat_val = samples['beamgood'][f'BDT_output_{HNL_mass}MeV']
        if logit == True:
            bkg_scores=[Functions.logit(samples['overlay_test'][f'BDT_output_{HNL_mass}MeV']),
                        Functions.logit(samples['dirtoverlay'][f'BDT_output_{HNL_mass}MeV']),
                        Functions.logit(samples['beamoff'][f'BDT_output_{HNL_mass}MeV'])]
            dat_val=Functions.logit(samples['beamgood'][f'BDT_output_{HNL_mass}MeV'])
        bkg_weights=[sample_norms['overlay_test'],sample_norms['dirtoverlay'],sample_norms['beamoff']]
        bkg_colors=[colours['overlay_test'],colours['dirtoverlay'],colours['beamoff']]
        labels=[fr"In-Cryo $\nu$",fr"Out-Cryo $\nu$",f"Beam-Off"]
        data_label = "NuMI Data"
        
        x,y=np.histogram(dat_val,bins=bins,range=xlims,density=density)
        x1,y=np.histogram(dat_val,bins=bins,range=xlims)
        bin_center = [(y[i] + y[i+1])/2. for i in range(len(y)-1)]
        dat_placeholder=x
        # dat_err=np.sqrt(x1)*Functions.safe_div(x,x1) #need to write one for arrays instead of single values.
        dat_err=np.sqrt(x1)*np.nan_to_num(x/x1)
        
        plt.errorbar(bin_center,dat_placeholder,yerr=dat_err,fmt='.',color='black',lw=5,capsize=5,elinewidth=3,label=data_label) #Plotting data
        
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
        
        #----sub-plot----#
        plt.sca(ax[1])
    
        fracer_data=np.nan_to_num(np.sqrt(x1)/x1)
        x_err=fracer_data*x
        fracer_mc=np.nan_to_num(tot_mcerr/plot[0][2])

        rat_err_data=x_err*(1/plot[0][2])

        rat_err_mc=fracer_mc
        rat_err=np.sqrt(rat_err_data**2)

        rat_err_mc=np.nan_to_num(rat_err_mc) #other wise the next doesnt plot pro[erly]

        upvals= np.append(1+( rat_err_mc),1+( rat_err_mc)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
        lowvals=np.append(1-( rat_err_mc),1-( rat_err_mc)[-1])


        plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)

        rat=np.nan_to_num(x/plot[0][2])
        rat[x==0]=1 #dont think this is a good way to deal with this

        rat_err=np.nan_to_num(rat*np.sqrt(fracer_mc**2+fracer_data**2))

        plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
        plt.ylabel("Data/MC")
        plt.axhline(1,ls='-',color='black')
        plt.axhline(1.1,ls='--',color='grey')
        plt.axhline(0.9,ls='--',color='grey')
        ylim = max(abs(np.nan_to_num(rat)))*1.1
        plt.ylim(0.7,1.3)
        
        if savefig == True:
            plt.savefig("plots/BDT_output/BDT_output_" + Run + "_" + str(HNL_mass) + "_MeV_" + logscale + save_str + ".pdf")
            plt.savefig("plots/BDT_output/BDT_output_" + Run + "_" + str(HNL_mass) + "_MeV_" + logscale + save_str + ".png")
        # plt.show()
        
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
        # tot_uncrt = 
              
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