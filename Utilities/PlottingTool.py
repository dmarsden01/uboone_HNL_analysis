Gimport matplotlib.pyplot as plt
import awkward as awk
import matplotlib as mpl
import numpy as np
import math
import Utilities.SystematicsCalc as SC
import Utilities.Consts as C

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





def SysPlot(var,query,xlabel=[],sample_info=[],xlims=[0,0],nbins=40,density=False,discrete=False,whichdf="vert_df",sortvar=[],dropdupes=False,figsize=[10,10],dpi=100):
    #dirt_norm_err 1 means 100%
    if(sample_info==[]): raise Exception("Specify sample_info dict") 
    if(xlabel==[]): xlabel=var
    
    if(whichdf=="daughters"):
        query=query+" & "+"daughter==0" #double count otherwise and you keep forgetting numbskull
    
        
    if(dropdupes):
        
        if(sortvar==[]):
            overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
        else:
    
            overlay=sample_info["overlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])

    else:                                                                               

        overlay=sample_info["overlay"][whichdf].query(query)


    var_Overlay=overlay[var]

    weight_Overlay=overlay["weight"]*sample_info["overlay"]["NormScale"]

    if xlims[0] == 0 and xlims[1] == 0: xlims = [min(var_Data),max(var_Data)]

    

    sys_dict={}
    covs=["Genie","PPFX","Reint"]
    for key in covs:
    
        sys_dict[key]={}
        temp=SC.sys_err(overlay,"weights"+key,var,query,xlims,nbins,sample_info["overlay"]["NormScale"])
        sys_dict[key]["cov"]=temp[0]
        cv=temp[1]
        sys_dict[key]["n_tot"]=temp[2]
        bins=temp[3]

    
    bins_cent=(bins[:-1]+bins[1:])/2
    bins_centlong=np.tile(bins_cent,600)

    plt.step(bins_cent,np.sqrt(np.diag(sys_dict["Genie"]["cov"]))/cv,label="Genie MultiSims")
    plt.step(bins_cent,np.sqrt(np.diag(sys_dict["PPFX"]["cov"]))/cv,label="PPFX MultiSims")
    plt.step(bins_cent,np.sqrt(np.diag(sys_dict["Reint"]["cov"]))/cv,label="G4 Reint MultiSims")
    plt.step(bins_cent,np.sqrt(np.diag(sys_dict["Reint"]["cov"]+sys_dict["Genie"]["cov"]+sys_dict["PPFX"]["cov"]))/cv,label="Total",color="black")
    plt.legend(frameon=False)
    plt.ylabel("Frac Err")
    
    plt.show()
    
    
    
    for key,v in sys_dict.items():
        origmat=v["cov"].copy()
        fraccovmat=v["cov"].copy()
        cormat=v["cov"].copy()
        
        for i in range(0,len(origmat)):
            for j in range(0,len(origmat[0])):
                cormat[i,j]=origmat[i,j]/(np.sqrt(origmat[i,i])*np.sqrt(origmat[j,j]))
                fraccovmat[i,j]=origmat[i,j]/(cv[i]*cv[j])
        plt.matshow(cormat) #corrlation matrix
        plt.colorbar()
        plt.show()
        plt.matshow(fraccovmat)
        plt.colorbar()
        plt.show()

def PlotVarible(var,query,xlabel=[],sample_info=[],xlims=[0,0],bins=40,MergeBins=False,StackSig=False,density=False,unblind=False,HNLplotscale="",HNLnumscale=1,datsigrar=3,HNLmasses=[],discrete=False,ord_flip=False,whichdf="vert_df",sortvar=[],dropdupes=False,ratio=False,CalcSys=False,PlotWellReco=True,figsize=[10,10],dpi=100,dirt_norm_err_fac=1,legnummode="inlim",legloc="best"):
    #dirt_norm_err 1 means 100%
    if(sample_info==[]): raise Exception("Specify sample_info dict") 
    if(xlabel==[]): xlabel=var
#     plt.figure(figsize=[12,12])

    q_control_region= "flash_time<6.900"

    if(not unblind):

        query=query+" & "+q_control_region
    
    if(whichdf=="daughters"):
        query=query+" & "+"daughter==0" #double count otherwise and you keep forgetting numbskull
        PlotWellReco=False
    
    
   
        
    if(dropdupes):
        
        if(sortvar==[]):
            beamgood=sample_info["beamgood"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            beamoff=sample_info["beamoff"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])

#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]
        else:
    
            beamgood=sample_info["beamgood"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            beamoff=sample_info["beamoff"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            overlay=sample_info["overlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])

#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]


    else:                                                                               
        beamgood=sample_info["beamgood"][whichdf].query(query)
        beamoff=sample_info["beamoff"][whichdf].query(query)
        overlay=sample_info["overlay"][whichdf].query(query)
        dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query)


#         weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#         weight_Overlay=sample_info["overlay"][whichdf].query(query)["weight"]*sample_info["overlay"]["NormScale"]
#         weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query)["weight"]*sample_info["dirtoverlay"]["NormScale"]

        
    var_Data=beamgood[var]
    var_Offbeam=beamoff[var]
    var_Overlay=overlay[var]
    var_Dirt=dirtoverlay[var]
    
    weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
    weight_Overlay=overlay["weight"]*sample_info["overlay"]["NormScale"]
    weight_Dirt=dirtoverlay["weight"]*sample_info["dirtoverlay"]["NormScale"]

    if xlims[0] == 0 and xlims[1] == 0: xlims = [min(var_Data),max(var_Data)]

#     if xlims[0]<min(var_Data): xlims[0]=min(var_Data)
#     if xlims[1]>max(var_Data): xlims[1]=max(var_Data)
        
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
                
#             else:
#                 print("empty bin,",bins[i],bins[i+1])
        bins_new.append(bins[-1])

        bins=bins_new
        

#     if xlims[0] == 0 and xlims[1] == 0: xlims = [min(min(va),min(var_MC)),max(max(var_Data),max(var_MC))]

    
    fig,ax = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=figsize,dpi=dpi)
    
    
    
    plt.sca(ax[0])
    
    if(dropdupes): plt.ylabel("Events")
    else: plt.ylabel("Vertices")
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
    dat_err=np.sqrt(x1)*np.nan_to_num(x/x1)
    
    if(legnummode=="inlim"): #this decides if the legend shows number just for the bins shown or the whole distrubtuin
        Offbeamnum=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam)[0].sum()
        Dirtnum=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt)[0].sum()
        Overlaynum=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay)[0].sum()
        Datanum=dat_val.sum()
    else:
        Offbeamnum=sum(weight_Offbeam)
        Dirtnum=sum(weight_Dirt)
        Overlaynum=sum(weight_Overlay)
        Datanum=len(var_Data)
        
    plt.plot([], [], ' ', label=f"Data/MC = {Datanum/(Dirtnum+Overlaynum+Offbeamnum):.2f}")
    
    plt.errorbar(bin_center,dat_val,yerr=dat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label=f"NuMI Data ({Datanum:.0f})")
    if(ord_flip):
        varis=[var_Offbeam,var_Overlay,var_Dirt]
        weights=[weight_Offbeam,weight_Overlay,weight_Dirt]
        colors=['sandybrown','seagreen',"darkgreen"]
        labels=[f"Beam-Off ({Offbeamnum:.1f})",fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})"] 
    else:
        varis=[var_Overlay,var_Dirt,var_Offbeam]
        weights=[weight_Overlay,weight_Dirt,weight_Offbeam]
        colors=['seagreen',"darkgreen",'sandybrown']
        labels=[fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})",f"Beam-Off ({Offbeamnum:.1f})"]


    plot=plt.hist(varis,
                  range=xlims,bins=bins,
                  histtype="stepfilled",
                  stacked=True,density=density,linewidth=2,edgecolor="darkblue",
                  weights=weights, color=colors)

    plot=plt.hist(varis,
              label=labels,
              range=xlims,bins=bins,
              histtype="stepfilled",
              stacked=True,density=density,linewidth=2,edgecolor="darkblue",
              weights=weights, color=colors)
#     mc=np.histogram(var_Overlay,bins=bins,range=xlims)
#     off=np.histogram(var_Offbeam,bins=bins,range=xlims)
    dirt=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt)

    mc_w=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay**2)
    off_w=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam**2)
    dirt_w=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt**2)
    
    off_err=np.sqrt(off_w[0])
    mc_err=np.sqrt(mc_w[0])
    dirt_err=np.sqrt(dirt_w[0])

    
    stat_mcerr=np.sqrt(off_err**2+mc_err**2+dirt_err**2)
    tot_mcerr=np.sqrt(off_err**2+mc_err**2+dirt_err**2)
#     print(tot_mcerr/plot[0][2])
    if(CalcSys):
        cov_gen,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsGenie",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DoneGenie")
        cov_PPFX,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsPPFX",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DonePPFX")
        cov_Reint,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsReint",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DoneReint")
        cov_mc_stat   = np.zeros([len(stat_mcerr), len(stat_mcerr)])
        cov_mc_stat[np.diag_indices_from(cov_mc_stat)]=stat_mcerr**2
        
        
#         print("gen",np.diag(cov_gen))
#         print("stat",np.diag(cov_mc_stat))
#         print("PPFX",np.diag(cov_PPFX))
        
        dirt_norm_err=dirt[0]*dirt_norm_err_fac
        tot_mcerr=np.sqrt( np.diag((cov_gen+ cov_mc_stat + cov_PPFX+cov_Reint))+dirt_norm_err**2) 
        
#     print(tot_mcerr/plot[0][2])
    upvals= np.append((plot[0][2]+tot_mcerr),(plot[0][2]+tot_mcerr)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
    lowvals=np.append((plot[0][2]-tot_mcerr),(plot[0][2]-tot_mcerr)[-1])

    plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)

    for mass in HNLmasses:
#         theta_u2=sample_info[mass]["theta_u2"]
        
        if(dropdupes):
            
            if(sortvar==[]):
                var_HNL=sample_info[mass][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])[var]
                if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco).drop_duplicates(subset=["run","evt","sub"])[var]
            else:
                
                var_HNL=sample_info[mass][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])[var]
                if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])[var]
#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]
        else:
    
            
            var_HNL=sample_info[mass][whichdf].query(query)[var]
            if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco)[var]
        
        
        if(legnummode=="inlim"):
            HNL_num=np.histogram(var_HNL,bins=bins,range=xlims)[0].sum()*sample_info[mass]["NormScale"]
            if (PlotWellReco): HNL_wrnum=np.histogram(var_HNL_wellreco,bins=bins,range=xlims)[0].sum()*sample_info[mass]["NormScale"]
        else:
            HNL_num=len(var_HNL)*sample_info[mass]["NormScale"]
            if (PlotWellReco): HNL_wrnum=len(var_HNL_wellreco)*sample_info[mass]["NormScale"]
        
        
        if(HNLplotscale==""): HNLplotscale=len(var_Data)/(datsigrar*(len(var_HNL)*sample_info[mass]["NormScale"]))
#         if(HNLplotscale==""): HNLplotscale=len(var_Data)/(datsigrar*(len(var_HNL)*sample_info[mass]["NormScale"]))
        
        bkg_stack=[]
        bkg_stack_w=[]
        color="red"
        color_wr="purple"
        if(StackSig):
            bkg_stack=varis
            bkg_stack_w=weights
            color=["#FF000000","#FF000000","#FF000000","red"]
            color_wr=["#FF000000","#FF000000","#FF000000","purple"]
        plt.hist(bkg_stack+[var_HNL],
#               label=[f"HNL ({mass} MeV) \n $|U_{{\mu4}}|^2="+sci_notation(sample_info["300"]["theta_u2"]) +f" (x{HNLplotscale})"],
              label=[f"HNL ({HNL_num*HNLnumscale:.1f})"],
              range=xlims,bins=bins,
              stacked=True,density=density,
              weights= bkg_stack_w+[np.ones(len(var_HNL))*sample_info[mass]["NormScale"]*HNLplotscale],histtype="step",color=color,lw=4)

        if(PlotWellReco):
            plt.hist(bkg_stack+[var_HNL_wellreco],
                  label=[f"HNL (Complete) ({HNL_wrnum*HNLnumscale:.1f})"],
                  range=xlims,bins=bins,
                  stacked=True,density=density,
                  weights=bkg_stack_w+[np.ones(len(var_HNL_wellreco))*sample_info[mass]["NormScale"]*HNLplotscale],histtype="step",color=color_wr,lw=4)
       
    plt.xlim(xlims)
    plt.legend(loc=legloc,frameon=False)
    plt.sca(ax[1])
    
    
    
    fracer_data=np.sqrt(x1)/x1
    x_err=fracer_data*x
    fracer_mc=tot_mcerr/plot[0][2]


    if(ratio):
        rat=x/plot[0][2]
        rat[x==0]=1 #dont think this is a good way to deal with this

        rat_err=rat*np.sqrt(fracer_mc**2+fracer_data**2)


        plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
        plt.ylabel("Data/MC")
        plt.axhline(1,ls='-',color='black')
        plt.axhline(1.1,ls='--',color='grey')
        plt.axhline(0.9,ls='--',color='grey')
        ylim = max(abs(np.nan_to_num(rat)))*1.1
        plt.ylim(0.7,1.3)

    
    else:
        
        rat=(x-plot[0][2])/plot[0][2]
        rat[x==0]=0
        plt.ylabel("(Data-Pred)/Pred",fontsize=18)
#         print(x-plot[0][2])
#         print(plot[0][2])
#         print("rat",rat)
        rat_err_data=x_err*(1/plot[0][2])
    
    
        #rat_err_mc=(x/plot[0][2]**2)*tot_mcerr # old propagated
        rat_err_mc=fracer_mc
#         rat_err=np.sqrt(rat_err_data**2+rat_err_mc**2)
        rat_err=np.sqrt(rat_err_data**2)
        
#         upvals= np.append(+(tot_mcerr/x),+(tot_mcerr/x)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
#         lowvals=np.append(-(tot_mcerr/x),-(tot_mcerr/x)[-1])
    
        
        rat_err_mc=np.nan_to_num(rat_err_mc) #other wise the next doesnt plot pro[erly]
        
        
        upvals= np.append(+( rat_err_mc),+( rat_err_mc)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
        lowvals=np.append(-( rat_err_mc),-( rat_err_mc)[-1])
    
    
        plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)
        rat[ rat==0 ] = np.nan
        rat[ rat==math.inf ] = np.nan
        plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
        ylim = max(abs(np.nan_to_num(rat)))*1.2
#         print(ylim)
        plt.ylim(-ylim,ylim)
        plt.axhline(0,ls='-',color='black')
        
    
    plt.xlabel(xlabel)
    plt.xlim(xlims)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.sca(ax[0])
    
    return bins,dat_val,dat_err,plot[0][2],tot_mcerr


def PlotVarible_bdt(var,query,legfontsize=20,ptype="HNL",xlabel=[],sample_info=[],xlims=[0,0],bins=40,MergeBins=False,StackSig=False,density=False,unblind=False,HNLplotscale="",HNLnumscale=1,extrascale=5,datsigrar=3,HNLmasses=[],discrete=False,ord_flip=False,whichdf="vert_df",sortvar=[],dropdupes=False,ratio=False,CalcSys=False,PlotWellReco=True,figsize=[10,10],dpi=100,dirt_norm_err_fac=1,legnummode="inlim",legloc="best"):
    #dirt_norm_err 1 means 100%
    if(sample_info==[]): raise Exception("Specify sample_info dict") 
    if(xlabel==[]): xlabel=var
#     plt.figure(figsize=[12,12])

    q_control_region= "flash_time<6.900"

    if(not unblind):

        query=query+" & "+q_control_region
    
    if(whichdf=="daughters"):
        query=query+" & "+"daughter==0" #double count otherwise and you keep forgetting numbskull
        PlotWellReco=False
    
    
   
        
    if(dropdupes):
        
        if(sortvar==[]):
            beamgood=sample_info["beamgood"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            beamoff=sample_info["beamoff"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])

#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]
        else:
    
            beamgood=sample_info["beamgood"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            beamoff=sample_info["beamoff"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            overlay=sample_info["overlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])

#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]


    else:                                                                               
        beamgood=sample_info["beamgood"][whichdf].query(query)
        beamoff=sample_info["beamoff"][whichdf].query(query)
        overlay=sample_info["overlay"][whichdf].query(query)
        dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query)


#         weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#         weight_Overlay=sample_info["overlay"][whichdf].query(query)["weight"]*sample_info["overlay"]["NormScale"]
#         weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query)["weight"]*sample_info["dirtoverlay"]["NormScale"]

        
    var_Data=beamgood[var]
    var_Offbeam=beamoff[var]
    var_Overlay=overlay[var]
    var_Dirt=dirtoverlay[var]
    
    weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
    weight_Overlay=overlay["weight"]*sample_info["overlay"]["NormScale"]
    weight_Dirt=dirtoverlay["weight"]*sample_info["dirtoverlay"]["NormScale"]

    if xlims[0] == 0 and xlims[1] == 0: xlims = [min(var_Data),max(var_Data)]

#     if xlims[0]<min(var_Data): xlims[0]=min(var_Data)
#     if xlims[1]>max(var_Data): xlims[1]=max(var_Data)
        
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
                
#             else:
#                 print("empty bin,",bins[i],bins[i+1])
        bins_new.append(bins[-1])

        bins=bins_new
        

#     if xlims[0] == 0 and xlims[1] == 0: xlims = [min(min(va),min(var_MC)),max(max(var_Data),max(var_MC))]

    
    fig,ax = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=figsize,dpi=dpi)
    
    
    
    plt.sca(ax[0])
    
    if(dropdupes): plt.ylabel("Events")
    else: plt.ylabel("Vertices")
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
    dat_err=np.sqrt(x1)*np.nan_to_num(x/x1)
    
    if(legnummode=="inlim"): #this decides if the legend shows number just for the bins shown or the whole distrubtuin
        Offbeamnum=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam)[0].sum()
        Dirtnum=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt)[0].sum()
        Overlaynum=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay)[0].sum()
        Datanum=dat_val.sum()
    else:
        Offbeamnum=sum(weight_Offbeam)
        Dirtnum=sum(weight_Dirt)
        Overlaynum=sum(weight_Overlay)
        Datanum=len(var_Data)
        
    plt.plot([], [], ' ', label=f"Data/MC = {Datanum/(Dirtnum+Overlaynum+Offbeamnum):.2f}")
    
    plt.errorbar(bin_center,dat_val,yerr=dat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label=f"NuMI Data ({Datanum:.0f})")
    if(ord_flip):
        varis=[var_Offbeam,var_Overlay,var_Dirt]
        weights=[weight_Offbeam,weight_Overlay,weight_Dirt]
        colors=['sandybrown','seagreen',"darkgreen"]
        labels=[f"Beam-Off ({Offbeamnum:.1f})",fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})"] 
    else:
        varis=[var_Overlay,var_Dirt,var_Offbeam]
        weights=[weight_Overlay,weight_Dirt,weight_Offbeam]
        colors=['seagreen',"darkgreen",'sandybrown']
        labels=[fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})",f"Beam-Off ({Offbeamnum:.1f})"]


    plot=plt.hist(varis,
                  range=xlims,bins=bins,
                  histtype="stepfilled",
                  stacked=True,density=density,linewidth=2,edgecolor="darkblue",
                  weights=weights, color=colors)

    plot=plt.hist(varis,
              label=labels,
              range=xlims,bins=bins,
              histtype="stepfilled",
              stacked=True,density=density,linewidth=2,edgecolor="darkblue",
              weights=weights, color=colors)
#     mc=np.histogram(var_Overlay,bins=bins,range=xlims)
#     off=np.histogram(var_Offbeam,bins=bins,range=xlims)
    dirt=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt)

    mc_w=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay**2)
    off_w=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam**2)
    dirt_w=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt**2)
    
    off_err=np.sqrt(off_w[0])
    mc_err=np.sqrt(mc_w[0])
    dirt_err=np.sqrt(dirt_w[0])

    
    stat_mcerr=np.sqrt(off_err**2+mc_err**2+dirt_err**2)
    tot_mcerr=np.sqrt(off_err**2+mc_err**2+dirt_err**2)
#     print(tot_mcerr/plot[0][2])
    if(CalcSys):
        cov_gen,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsGenie",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DoneGenie")
        cov_PPFX,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsPPFX",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DonePPFX")
        cov_Reint,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsReint",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DoneReint")
        cov_mc_stat   = np.zeros([len(stat_mcerr), len(stat_mcerr)])
        cov_mc_stat[np.diag_indices_from(cov_mc_stat)]=stat_mcerr**2
        
        
#         print("gen",np.diag(cov_gen))
#         print("stat",np.diag(cov_mc_stat))
#         print("PPFX",np.diag(cov_PPFX))
        
        dirt_norm_err=dirt[0]*dirt_norm_err_fac
        tot_mcerr=np.sqrt( np.diag((cov_gen+ cov_mc_stat + cov_PPFX+cov_Reint))+dirt_norm_err**2) 
        
#     print(tot_mcerr/plot[0][2])
    upvals= np.append((plot[0][2]+tot_mcerr),(plot[0][2]+tot_mcerr)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
    lowvals=np.append((plot[0][2]-tot_mcerr),(plot[0][2]-tot_mcerr)[-1])

    plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)

    for mass in HNLmasses:
#         theta_u2=sample_info[mass]["theta_u2"]
        
        if(dropdupes):
            
            if(sortvar==[]):
                var_HNL=sample_info[mass][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])[var]
                if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco).drop_duplicates(subset=["run","evt","sub"])[var]
            else:
                
                var_HNL=sample_info[mass][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])[var]
                if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])[var]
#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]
        else:
    
            
            var_HNL=sample_info[mass][whichdf].query(query)[var]
            if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco)[var]
        
        
        if(legnummode=="inlim"):
            HNL_num=np.histogram(var_HNL,bins=bins,range=xlims)[0].sum()*sample_info[mass]["NormScale"]
            if (PlotWellReco): HNL_wrnum=np.histogram(var_HNL_wellreco,bins=bins,range=xlims)[0].sum()*sample_info[mass]["NormScale"]
        else:
            HNL_num=len(var_HNL)*sample_info[mass]["NormScale"]
            if (PlotWellReco): HNL_wrnum=len(var_HNL_wellreco)*sample_info[mass]["NormScale"]
        
        
        if(HNLplotscale==""): HNLplotscale=len(var_Data)/(datsigrar*(len(var_HNL)*sample_info[mass]["NormScale"]))
#         if(HNLplotscale==""): HNLplotscale=len(var_Data)/(datsigrar*(len(var_HNL)*sample_info[mass]["NormScale"]))
        
        bkg_stack=[]
        bkg_stack_w=[]
        if(ptype=="HNL"): color="red"
        else: color= "tab:blue"
        color_wr="purple"
        if(StackSig):
            bkg_stack=varis
            bkg_stack_w=weights
            color=["#FF000000","#FF000000","#FF000000",color]
            color_wr=["#FF000000","#FF000000","#FF000000","purple"]
        plt.hist(bkg_stack+[var_HNL],
#               label=[f"HNL ({mass} MeV) \n $|U_{{\mu4}}|^2="+sci_notation(sample_info["300"]["theta_u2"]) +f" (x{HNLplotscale})"],
              label=["","","",fr"{ptype} ({HNL_num*HNLnumscale:.1f}) $\times{extrascale}$"],
              range=xlims,bins=bins,
              stacked=True,density=density,
              weights= bkg_stack_w+[np.ones(len(var_HNL))*sample_info[mass]["NormScale"]*HNLplotscale],color=color,histtype="step",lw=4)

        if(PlotWellReco):
            plt.hist(bkg_stack+[var_HNL_wellreco],
                  label=[f"HNL (Complete) ({HNL_wrnum*HNLnumscale:.1f})"],
                  range=xlims,bins=bins,
                  stacked=True,density=density,
                  weights=bkg_stack_w+[np.ones(len(var_HNL_wellreco))*sample_info[mass]["NormScale"]*HNLplotscale],histtype="step",color=color_wr,lw=4)
       
    plt.xlim(xlims)
    plt.legend(loc=legloc,frameon=False,fontsize=legfontsize)
    plt.sca(ax[1])
    
    
    
    fracer_data=np.sqrt(x1)/x1
    x_err=fracer_data*x
    fracer_mc=tot_mcerr/plot[0][2]


    if(ratio):
        rat=x/plot[0][2]
        rat[x==0]=1 #dont think this is a good way to deal with this

        rat_err=rat*np.sqrt(fracer_mc**2+fracer_data**2)


        plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
        plt.ylabel("Data/MC")
        plt.axhline(1,ls='-',color='black')
        plt.axhline(1.1,ls='--',color='grey')
        plt.axhline(0.9,ls='--',color='grey')
        ylim = max(abs(np.nan_to_num(rat)))*1.1
        plt.ylim(0.7,1.3)

    
    else:
        
        rat=(x-plot[0][2])/plot[0][2]
        rat[x==0]=0
        plt.ylabel("(Data-Pred)/Pred",fontsize=18)
#         print(x-plot[0][2])
#         print(plot[0][2])
#         print("rat",rat)
        rat_err_data=x_err*(1/plot[0][2])
    
    
        #rat_err_mc=(x/plot[0][2]**2)*tot_mcerr # old propagated
        rat_err_mc=fracer_mc
#         rat_err=np.sqrt(rat_err_data**2+rat_err_mc**2)
        rat_err=np.sqrt(rat_err_data**2)
        
#         upvals= np.append(+(tot_mcerr/x),+(tot_mcerr/x)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
#         lowvals=np.append(-(tot_mcerr/x),-(tot_mcerr/x)[-1])
    
        
        rat_err_mc=np.nan_to_num(rat_err_mc) #other wise the next doesnt plot pro[erly]
        
        
        upvals= np.append(+( rat_err_mc),+( rat_err_mc)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
        lowvals=np.append(-( rat_err_mc),-( rat_err_mc)[-1])
    
    
        plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)
        rat[ rat==0 ] = np.nan
        rat[ rat==math.inf ] = np.nan
        plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
        ylim = max(abs(np.nan_to_num(rat)))*1.2
#         print(ylim)
        plt.ylim(-ylim,ylim)
        plt.axhline(0,ls='-',color='black')
        
    
    plt.xlabel(xlabel)
    plt.xlim(xlims)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.sca(ax[0])
    
    return bins,dat_val,dat_err,plot[0][2],tot_mcerr


def PlotVarible_new(var,query,xlabel=[],sample_info=[],xlims=[0,0],bins=40,MergeBins=False,StackSig=False,density=False,unblind=False,HNLplotscale="",datsigrar=3,HNLmasses=[],discrete=False,ord_flip=False,whichdf="vert_df",sortvar=[],dropdupes=False,ratio=False,CalcSys=False,PlotWellReco=True,figsize=[10,10],dpi=100,dirt_norm_err_fac=1,legnummode="inlim",legloc="best"):
    #dirt_norm_err 1 means 100%
    if(sample_info==[]): raise Exception("Specify sample_info dict") 
    if(xlabel==[]): xlabel=var
#     plt.figure(figsize=[12,12])

    q_control_region= "flash_time<6.900"

    if(not unblind):

        query=query+" & "+q_control_region
    
    if(whichdf=="daughters"):
        query=query+" & "+"daughter==0" #double count otherwise and you keep forgetting numbskull
        PlotWellReco=False
    
    
   
        
    if(dropdupes):
        
        if(sortvar==[]):
            beamgood=sample_info["beamgood"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            beamoff=sample_info["beamoff"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])

#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]
        else:
    
            beamgood=sample_info["beamgood"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            beamoff=sample_info["beamoff"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            overlay=sample_info["overlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])

#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]


    else:                                                                               
        beamgood=sample_info["beamgood"][whichdf].query(query)
        beamoff=sample_info["beamoff"][whichdf].query(query)
        overlay=sample_info["overlay"][whichdf].query(query)
        dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query)


#         weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#         weight_Overlay=sample_info["overlay"][whichdf].query(query)["weight"]*sample_info["overlay"]["NormScale"]
#         weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query)["weight"]*sample_info["dirtoverlay"]["NormScale"]

        
    var_Data=beamgood[var]
    var_Offbeam=beamoff[var]
    var_Overlay=overlay[var]
    var_Dirt=dirtoverlay[var]
    
    weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
    weight_Overlay=overlay["weight"]*sample_info["overlay"]["NormScale"]
    weight_Dirt=dirtoverlay["weight"]*sample_info["dirtoverlay"]["NormScale"]

    if xlims[0] == 0 and xlims[1] == 0: xlims = [min(var_Data),max(var_Data)]

#     if xlims[0]<min(var_Data): xlims[0]=min(var_Data)
#     if xlims[1]>max(var_Data): xlims[1]=max(var_Data)
        
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
            print(offbkg[i],overlaybkg[i])
            if(offbkg[i]>0 or overlaybkg[i]>0):
                bins_new.append(bins[i])
                
#             else:
#                 print("empty bin,",bins[i],bins[i+1])
        bins_new.append(bins[-1])

        bins=bins_new
        

#     if xlims[0] == 0 and xlims[1] == 0: xlims = [min(min(va),min(var_MC)),max(max(var_Data),max(var_MC))]

    
    fig,ax = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=figsize,dpi=dpi)
    
    
    
    plt.sca(ax[0])
    
    if(dropdupes): plt.ylabel("Events")
    else: plt.ylabel("Vertices")
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
    dat_err=np.sqrt(x1)*np.nan_to_num(x/x1)
    
    if(legnummode=="inlim"): #this decides if the legend shows number just for the bins shown or the whole distrubtuin
        Offbeamnum=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam)[0].sum()
        Dirtnum=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt)[0].sum()
        Overlaynum=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay)[0].sum()
        Datanum=dat_val.sum()
    else:
        Offbeamnum=sum(weight_Offbeam)
        Dirtnum=sum(weight_Dirt)
        Overlaynum=sum(weight_Overlay)
        Datanum=len(var_Data)
        
    plt.plot([], [], ' ', label=f"Data/MC = {Datanum/(Dirtnum+Overlaynum+Offbeamnum):.2f}")
    
    plt.errorbar(bin_center,dat_val,yerr=dat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label=f"NuMI Data ({Datanum:.0f})")
    if(ord_flip):
        varis=[var_Offbeam,var_Overlay,var_Dirt]
        weights=[weight_Offbeam,weight_Overlay,weight_Dirt]
        colors=['sandybrown','seagreen',"darkgreen"]
        labels=[f"Beam-Off ({Offbeamnum:.1f})",fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})"] 
    else:
        varis=[var_Overlay,var_Dirt,var_Offbeam]
        weights=[weight_Overlay,weight_Dirt,weight_Offbeam]
        colors=['seagreen',"darkgreen",'sandybrown']
        labels=[fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})",f"Beam-Off ({Offbeamnum:.1f})"]


    plot=plt.hist(varis,
                  range=xlims,bins=bins,
                  histtype="stepfilled",
                  stacked=True,density=density,linewidth=2,edgecolor="darkblue",
                  weights=weights, color=colors)

    plot=plt.hist(varis,
              label=labels,
              range=xlims,bins=bins,
              histtype="stepfilled",
              stacked=True,density=density,linewidth=2,edgecolor="darkblue",
              weights=weights, color=colors)
#     mc=np.histogram(var_Overlay,bins=bins,range=xlims)
#     off=np.histogram(var_Offbeam,bins=bins,range=xlims)
    dirt=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt)

    mc_w=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay**2)
    off_w=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam**2)
    dirt_w=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt**2)
    
    off_err=np.sqrt(off_w[0])
    mc_err=np.sqrt(mc_w[0])
    dirt_err=np.sqrt(dirt_w[0])

    
    stat_mcerr=np.sqrt(off_err**2+mc_err**2+dirt_err**2)
    tot_mcerr=np.sqrt(off_err**2+mc_err**2+dirt_err**2)
#     print(tot_mcerr/plot[0][2])
    if(CalcSys):
        cov_gen,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsGenie",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DoneGenie")
        cov_PPFX,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsPPFX",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DonePPFX")
        cov_Reint,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsReint",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DoneReint")
        cov_mc_stat   = np.zeros([len(stat_mcerr), len(stat_mcerr)])
        cov_mc_stat[np.diag_indices_from(cov_mc_stat)]=stat_mcerr**2
        
        
#         print("gen",np.diag(cov_gen))
#         print("stat",np.diag(cov_mc_stat))
#         print("PPFX",np.diag(cov_PPFX))
        
        dirt_norm_err=dirt[0]*dirt_norm_err_fac
        tot_mcerr=np.sqrt( np.diag((cov_gen+ cov_mc_stat + cov_PPFX+cov_Reint))+dirt_norm_err**2) 
        
#     print(tot_mcerr/plot[0][2])
    upvals= np.append((plot[0][2]+tot_mcerr),(plot[0][2]+tot_mcerr)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
    lowvals=np.append((plot[0][2]-tot_mcerr),(plot[0][2]-tot_mcerr)[-1])

    plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)

    for mass in HNLmasses:
#         theta_u2=sample_info[mass]["theta_u2"]
        
        if(dropdupes):
            
            if(sortvar==[]):
                var_HNL=sample_info[mass][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])[var]
                if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco).drop_duplicates(subset=["run","evt","sub"])[var]
            else:
                
                var_HNL=sample_info[mass][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])[var]
                if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])[var]
#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]
        else:
    
            
            var_HNL=sample_info[mass][whichdf].query(query)[var]
            if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco)[var]
        
        
        if(legnummode=="inlim"):
            HNL_num=np.histogram(var_HNL,bins=bins,range=xlims)[0].sum()*sample_info[mass]["NormScale"]
            if (PlotWellReco): HNL_wrnum=np.histogram(var_HNL_wellreco,bins=bins,range=xlims)[0].sum()*sample_info[mass]["NormScale"]
        else:
            HNL_num=len(var_HNL)*sample_info[mass]["NormScale"]
            if (PlotWellReco): HNL_wrnum=len(var_HNL_wellreco)*sample_info[mass]["NormScale"]
        
        
        if(HNLplotscale==""): HNLplotscale=len(var_Data)/(datsigrar*(len(var_HNL)*sample_info[mass]["NormScale"]))
#         if(HNLplotscale==""): HNLplotscale=len(var_Data)/(datsigrar*(len(var_HNL)*sample_info[mass]["NormScale"]))
        
        bkg_stack=[]
        bkg_stack_w=[]
        color="red"
        color_wr="purple"
        if(StackSig):
            bkg_stack=varis
            bkg_stack_w=weights
            color=["#FF000000","#FF000000","#FF000000","red"]
            color_wr=["#FF000000","#FF000000","#FF000000","purple"]
        plt.hist(bkg_stack+[var_HNL],
#               label=[f"HNL ({mass} MeV) \n $|U_{{\mu4}}|^2="+sci_notation(sample_info["300"]["theta_u2"]) +f" (x{HNLplotscale})"],
              label=[f"HNL ({HNL_num:.1f})"],
              range=xlims,bins=bins,
              stacked=True,density=density,
              weights= bkg_stack_w+[np.ones(len(var_HNL))*sample_info[mass]["NormScale"]*HNLplotscale],histtype="step",color=color,lw=4)

        if(PlotWellReco):
            plt.hist(bkg_stack+[var_HNL_wellreco],
                  label=[f"HNL (Complete) ({HNL_wrnum:.1f})"],
                  range=xlims,bins=bins,
                  stacked=True,density=density,
                  weights=bkg_stack_w+[np.ones(len(var_HNL_wellreco))*sample_info[mass]["NormScale"]*HNLplotscale],histtype="step",color=color_wr,lw=4)
       
    plt.xlim(xlims)
    plt.legend(loc=legloc,frameon=False, fontsize=legfontsize)
    plt.sca(ax[1])
    
    
    
    fracer_data=np.sqrt(x1)/x1
    x_err=fracer_data*x
    fracer_mc=tot_mcerr/plot[0][2]


    if(ratio):
        rat=x/plot[0][2]
        rat[x==0]=1 #dont think this is a good way to deal with this

        rat_err=rat*np.sqrt(fracer_mc**2+fracer_data**2)


        plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
        plt.ylabel("Data/MC")
        plt.axhline(1,ls='-',color='black')
        plt.axhline(1.1,ls='--',color='grey')
        plt.axhline(0.9,ls='--',color='grey')
        ylim = max(abs(np.nan_to_num(rat)))*1.1
        plt.ylim(0.7,1.3)

    
    else:
        
        rat=(x-plot[0][2])/plot[0][2]
        rat[x==0]=0
        plt.ylabel("(Data-Pred)/Pred",fontsize=18)
#         print(x-plot[0][2])
#         print(plot[0][2])
#         print("rat",rat)
        rat_err_data=x_err*(1/plot[0][2])
    
    
        #rat_err_mc=(x/plot[0][2]**2)*tot_mcerr # old propagated
        rat_err_mc=fracer_mc
#         rat_err=np.sqrt(rat_err_data**2+rat_err_mc**2)
        rat_err=np.sqrt(rat_err_data**2)
        
#         upvals= np.append(+(tot_mcerr/x),+(tot_mcerr/x)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
#         lowvals=np.append(-(tot_mcerr/x),-(tot_mcerr/x)[-1])
    
        
        rat_err_mc=np.nan_to_num(rat_err_mc) #other wise the next doesnt plot pro[erly]
        
        
        upvals= np.append(+( rat_err_mc),+( rat_err_mc)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
        lowvals=np.append(-( rat_err_mc),-( rat_err_mc)[-1])
    
    
        plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)
        rat[ rat==0 ] = np.nan
        rat[ rat==math.inf ] = np.nan
        plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
        ylim = max(abs(np.nan_to_num(rat)))*1.2
#         print(ylim)
        plt.ylim(-ylim,ylim)
        plt.axhline(0,ls='-',color='black')
        
    
    plt.xlabel(xlabel)
    plt.xlim(xlims)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.sca(ax[0])
    
    return bins,dat_val,dat_err,plot[0][2],tot_mcerr
    
def PlotVarible_nodat(var,query,xlabel=[],sample_info=[],xlims=[0,0],bins=40,MergeBins=False,StackSig=False,density=False,unblind=False,HNLplotscale="",datsigrar=3,HNLmasses=[],discrete=False,ord_flip=False,whichdf="vert_df",sortvar=[],dropdupes=False,ratio=False,CalcSys=False,PlotWellReco=True,figsize=[10,10],dpi=100,dirt_norm_err_fac=1,legnummode="inlim",legloc="best"):
    #dirt_norm_err 1 means 100%
    if(sample_info==[]): raise Exception("Specify sample_info dict") 
    if(xlabel==[]): xlabel=var
#     plt.figure(figsize=[12,12])

    q_control_region= "flash_time<6.900"

    if(not unblind):

        query=query+" & "+q_control_region
    
    if(whichdf=="daughters"):
        query=query+" & "+"daughter==0" #double count otherwise and you keep forgetting numbskull
        PlotWellReco=False
    
    
   
        
    if(dropdupes):
        
        if(sortvar==[]):
            beamgood=sample_info["beamgood"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            beamoff=sample_info["beamoff"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])
            dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])

#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]
        else:
    
            beamgood=sample_info["beamgood"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            beamoff=sample_info["beamoff"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            overlay=sample_info["overlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])
            dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])

#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]


    else:                                                                               
        beamgood=sample_info["beamgood"][whichdf].query(query)
        beamoff=sample_info["beamoff"][whichdf].query(query)
        overlay=sample_info["overlay"][whichdf].query(query)
        dirtoverlay=sample_info["dirtoverlay"][whichdf].query(query)


#         weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#         weight_Overlay=sample_info["overlay"][whichdf].query(query)["weight"]*sample_info["overlay"]["NormScale"]
#         weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query)["weight"]*sample_info["dirtoverlay"]["NormScale"]

        
    var_Data=beamgood[var]
    var_Offbeam=beamoff[var]
    var_Overlay=overlay[var]
    var_Dirt=dirtoverlay[var]
    
    weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
    weight_Overlay=overlay["weight"]*sample_info["overlay"]["NormScale"]
    weight_Dirt=dirtoverlay["weight"]*sample_info["dirtoverlay"]["NormScale"]

    if xlims[0] == 0 and xlims[1] == 0: xlims = [min(var_Data),max(var_Data)]

#     if xlims[0]<min(var_Data): xlims[0]=min(var_Data)
#     if xlims[1]>max(var_Data): xlims[1]=max(var_Data)
        
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
                
#             else:
#                 print("empty bin,",bins[i],bins[i+1])
        bins_new.append(bins[-1])

        bins=bins_new
        

#     if xlims[0] == 0 and xlims[1] == 0: xlims = [min(min(va),min(var_MC)),max(max(var_Data),max(var_MC))]

    
    fig,ax = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=figsize,dpi=dpi)
    
    
    
    plt.sca(ax[0])
    
    if(dropdupes): plt.ylabel("Events")
    else: plt.ylabel("Vertices")
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
    dat_err=np.sqrt(x1)*np.nan_to_num(x/x1)
    
    if(legnummode=="inlim"): #this decides if the legend shows number just for the bins shown or the whole distrubtuin
        Offbeamnum=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam)[0].sum()
        Dirtnum=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt)[0].sum()
        Overlaynum=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay)[0].sum()
        Datanum=dat_val.sum()
    else:
        Offbeamnum=sum(weight_Offbeam)
        Dirtnum=sum(weight_Dirt)
        Overlaynum=sum(weight_Overlay)
        Datanum=len(var_Data)
        
#     plt.plot([], [], ' ', label=f"Data/MC = {Datanum/(Dirtnum+Overlaynum+Offbeamnum):.2f}")
    
#     plt.errorbar(bin_center,dat_val,yerr=dat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label=f"NuMI Data ({Datanum:.0f})")
    if(ord_flip):
        varis=[var_Offbeam,var_Overlay,var_Dirt]
        weights=[weight_Offbeam,weight_Overlay,weight_Dirt]
        colors=['sandybrown','seagreen',"darkgreen"]
        labels=[f"Beam-Off ({Offbeamnum:.1f})",fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})"] 
    else:
        varis=[var_Overlay,var_Dirt,var_Offbeam]
        weights=[weight_Overlay,weight_Dirt,weight_Offbeam]
        colors=['seagreen',"darkgreen",'sandybrown']
        labels=[fr"In-Cryo $\nu$ ({Overlaynum:.1f})",fr"Out-Cryo $\nu$ ({Dirtnum:.1f})",f"Beam-Off ({Offbeamnum:.1f})"]


    plot=plt.hist(varis,
                  range=xlims,bins=bins,
                  histtype="stepfilled",
                  stacked=True,density=density,linewidth=2,edgecolor="darkblue",
                  weights=weights, color=colors)

    plot=plt.hist(varis,
              label=labels,
              range=xlims,bins=bins,
              histtype="stepfilled",
              stacked=True,density=density,linewidth=2,edgecolor="darkblue",
              weights=weights, color=colors)
#     mc=np.histogram(var_Overlay,bins=bins,range=xlims)
#     off=np.histogram(var_Offbeam,bins=bins,range=xlims)
    dirt=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt)

    mc_w=np.histogram(var_Overlay,bins=bins,range=xlims,weights=weight_Overlay**2)
    off_w=np.histogram(var_Offbeam,bins=bins,range=xlims,weights=weight_Offbeam**2)
    dirt_w=np.histogram(var_Dirt,bins=bins,range=xlims,weights=weight_Dirt**2)
    
    off_err=np.sqrt(off_w[0])
    mc_err=np.sqrt(mc_w[0])
    dirt_err=np.sqrt(dirt_w[0])

    
    stat_mcerr=np.sqrt(off_err**2+mc_err**2+dirt_err**2)
    tot_mcerr=np.sqrt(off_err**2+mc_err**2+dirt_err**2)
#     print(tot_mcerr/plot[0][2])
    if(CalcSys):
        cov_gen,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsGenie",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DoneGenie")
        cov_PPFX,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsPPFX",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DonePPFX")
        cov_Reint,cv,n_tot,bins_sys=SC.sys_err(overlay,"weightsReint",var,query,xlims,bins,sample_info["overlay"]["NormScale"])
        print("DoneReint")
        cov_mc_stat   = np.zeros([len(stat_mcerr), len(stat_mcerr)])
        cov_mc_stat[np.diag_indices_from(cov_mc_stat)]=stat_mcerr**2
        
        
#         print("gen",np.diag(cov_gen))
#         print("stat",np.diag(cov_mc_stat))
#         print("PPFX",np.diag(cov_PPFX))
        
        dirt_norm_err=dirt[0]*dirt_norm_err_fac
        tot_mcerr=np.sqrt( np.diag((cov_gen+ cov_mc_stat + cov_PPFX+cov_Reint))+dirt_norm_err**2) 
        
#     print(tot_mcerr/plot[0][2])
    upvals= np.append((plot[0][2]+tot_mcerr),(plot[0][2]+tot_mcerr)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
    lowvals=np.append((plot[0][2]-tot_mcerr),(plot[0][2]-tot_mcerr)[-1])

    plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)

    for mass in HNLmasses:
#         theta_u2=sample_info[mass]["theta_u2"]
        
        if(dropdupes):
            
            if(sortvar==[]):
                var_HNL=sample_info[mass][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])[var]
                if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco).drop_duplicates(subset=["run","evt","sub"])[var]
            else:
                
                var_HNL=sample_info[mass][whichdf].query(query).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])[var]
                if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco).sort_values(sortvar,ascending=False).drop_duplicates(subset=["run","evt","sub"])[var]
#             weight_Offbeam=np.ones(len(var_Offbeam))*sample_info["beamoff"]["NormScale"]
#             weight_Overlay=sample_info["overlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["overlay"]["NormScale"]
#             weight_Dirt=sample_info["dirtoverlay"][whichdf].query(query).drop_duplicates(subset=["run","evt","sub"])["weight"]*sample_info["dirtoverlay"]["NormScale"]
        else:
    
            
            var_HNL=sample_info[mass][whichdf].query(query)[var]
            if(PlotWellReco): var_HNL_wellreco=sample_info[mass][whichdf].query(query+" & "+C.q_wellreco)[var]
        
        
        if(legnummode=="inlim"):
            HNL_num=np.histogram(var_HNL,bins=bins,range=xlims)[0].sum()*sample_info[mass]["NormScale"]
            if (PlotWellReco): HNL_wrnum=np.histogram(var_HNL_wellreco,bins=bins,range=xlims)[0].sum()*sample_info[mass]["NormScale"]
        else:
            HNL_num=len(var_HNL)*sample_info[mass]["NormScale"]
            if (PlotWellReco): HNL_wrnum=len(var_HNL_wellreco)*sample_info[mass]["NormScale"]
        
        
        if(HNLplotscale==""): HNLplotscale=len(var_Data)/(datsigrar*(len(var_HNL)*sample_info[mass]["NormScale"]))
#         if(HNLplotscale==""): HNLplotscale=len(var_Data)/(datsigrar*(len(var_HNL)*sample_info[mass]["NormScale"]))
        
        bkg_stack=[]
        bkg_stack_w=[]
        color="red"
        color_wr="purple"
        if(StackSig):
            bkg_stack=varis
            bkg_stack_w=weights
            color=["#FF000000","#FF000000","#FF000000","red"]
            color_wr=["#FF000000","#FF000000","#FF000000","purple"]
        plt.hist(bkg_stack+[var_HNL],
#               label=[f"HNL ({mass} MeV) \n $|U_{{\mu4}}|^2="+sci_notation(sample_info["300"]["theta_u2"]) +f" (x{HNLplotscale})"],
              label=[f"HNL ({HNL_num:.1f})"],
              range=xlims,bins=bins,
              stacked=True,density=density,
              weights= bkg_stack_w+[np.ones(len(var_HNL))*sample_info[mass]["NormScale"]*HNLplotscale],histtype="step",color=color,lw=4)

        if(PlotWellReco):
            plt.hist(bkg_stack+[var_HNL_wellreco],
                  label=[f"HNL (Complete) ({HNL_wrnum:.1f})"],
                  range=xlims,bins=bins,
                  stacked=True,density=density,
                  weights=bkg_stack_w+[np.ones(len(var_HNL_wellreco))*sample_info[mass]["NormScale"]*HNLplotscale],histtype="step",color=color_wr,lw=4)
       
    plt.xlim(xlims)
    plt.legend(loc=legloc,frameon=False)
    plt.sca(ax[1])
    
    
    
    fracer_data=np.sqrt(x1)/x1
    x_err=fracer_data*x
    fracer_mc=tot_mcerr/plot[0][2]


    if(ratio):
        rat=x/plot[0][2]
        rat[x==0]=1 #dont think this is a good way to deal with this

        rat_err=rat*np.sqrt(fracer_mc**2+fracer_data**2)


#         plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
        plt.ylabel("Data/MC")
        plt.axhline(1,ls='-',color='black')
        plt.axhline(1.1,ls='--',color='grey')
        plt.axhline(0.9,ls='--',color='grey')
        ylim = max(abs(np.nan_to_num(rat)))*1.1
        plt.ylim(0.7,1.3)

    
    else:
        
        rat=(x-plot[0][2])/plot[0][2]
        rat[x==0]=0
        plt.ylabel("Fractional Error",fontsize=18)
#         print(x-plot[0][2])
#         print(plot[0][2])
#         print("rat",rat)
        rat_err_data=x_err*(1/plot[0][2])
    
    
        #rat_err_mc=(x/plot[0][2]**2)*tot_mcerr # old propagated
        rat_err_mc=fracer_mc
#         rat_err=np.sqrt(rat_err_data**2+rat_err_mc**2)
        rat_err=np.sqrt(rat_err_data**2)
        
#         upvals= np.append(+(tot_mcerr/x),+(tot_mcerr/x)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
#         lowvals=np.append(-(tot_mcerr/x),-(tot_mcerr/x)[-1])
    
        
        rat_err_mc=np.nan_to_num(rat_err_mc) #other wise the next doesnt plot pro[erly]
        
        
        upvals= np.append(+( rat_err_mc),+( rat_err_mc)[-1]) #hate this but need to repeat last value to get bar on last bin to work, saw it here https://matplotlib.org/stable/gallery/lines_bars_and_markers/filled_step.html
        lowvals=np.append(-( rat_err_mc),-( rat_err_mc)[-1])
    
    
        plt.fill_between(y, lowvals, upvals,step="post",color="grey",alpha=0.3,zorder=2)
        rat[ rat==0 ] = np.nan
        rat[ rat==math.inf ] = np.nan
#         plt.errorbar(bin_center,rat,yerr=rat_err,fmt='.',color='black',lw=3,capsize=3,elinewidth=1,label="data")
        ylim = max(abs(np.nan_to_num( upvals)))*1.2
#         print(ylim)
        plt.ylim(-ylim,ylim)
        plt.axhline(0,ls='-',color='black')
        
    
    plt.xlabel(xlabel)
    plt.xlim(xlims)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.sca(ax[0])
    
    return bins,dat_val,dat_err,plot[0][2],tot_mcerr
    
       
    

