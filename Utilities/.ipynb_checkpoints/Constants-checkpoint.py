
Preselection_dict = {"nslice":"nslice==1",
"flash_time":"flash_time > 6.55 and flash_time < 16.5",
"nu_flashmatch_score":"nu_flashmatch_score < 15",#may need to take out.
"NeutrinoEnergy2":"NeutrinoEnergy2 < 500",
"contained_fraction":"contained_fraction > 0.9",
"trk_score_v":"trk_score_v < 0.97",
"n_pfps":"n_pfps < 6"
}

Preselection_dict_for_plot = {"nslice":"nslice==1",
"flash_time":"flash_time > 6.55 and flash_time < 16.5",
"nu_flashmatch_score":"nu_flashmatch_score < 15",#may need to take out.
"NeutrinoEnergy2":"NeutrinoEnergy2 < 500",
"contained_fraction":"contained_fraction > 0.9",
"trk_score":"trk_score < 0.97",
"n_pfps":"n_pfps < 6"
}

Preselection_dict_crtveto = {"crtveto":"crtveto==0"}

HNL_mass_samples = [20, 50, 100, 150, 180, 200]

theta_mu_4 = 1e-4 #This is the same for all samples

pi0_scaling_factor = 0.759

run1_POT_scaling_dict = {20: 2.1387912933253144e-12, 
                         50: 6.590170245207509e-10, 
                         100: 6.913852703765691e-08, 
                         150: 3.301953636715338e-06, 
                         180: 3.479531678617752e-05, 
                         200: 9.077313301681812e-05}

SF_overlay_run1 = 0.08559729531465318
SF_dirt_run1 = 0.08961042442406035
SF_EXT_run1 = 0.5612087579382191
SF_signal_run1_100MeV = 6.913852627300102e-08

run3_POT_scaling_dict = {20: 1.7657423469953374e-12, 
                         50: 5.448786613287057e-10, 
                         100: 5.5322831399190343e-08, 
                         150: 2.623088564712337e-06, 
                         180: 2.624339736997261e-05, 
                         200: 7.283055184784342e-05}

SF_overlay_run3 = 0.2513368817255014
SF_dirt_run3 = 0.16953052634982632
SF_EXT_run3 = 0.3089104916683624

sample_colours = {'overlay':'peru',
                  'dirtoverlay':'darkorange',
                  'beamoff':'deepskyblue',
                  'signal':'red'}

##TPC
TPCxlo=0
TPCxhi=256
TPCylo=-115
TPCyhi=117
TPCzlo=0
TPCzhi=1037



bdt_cols_all=[
    'n_tracks',
 'n_showers',
    'topological_score',
    'NeutrinoEnergy2',
    'trk1_len_v','trk2_len_v',
    "trk1_dect_phi",
    "trk1_dect_theta",
    'trk1_llr_pid_score_v',  'trk2_llr_pid_score_v',  
    'trk1_score_v', 'trk2_score_v',
          'cosangle',
       'max_x','min_x',  
    'max_y','max_z', 
         'min_y', 'min_z']

bdt_cols_hps=bdt_cols_all+[
    
    "Inv_mumu",
    "Ang_with_scalar_mumu"
]

bdt_cols_hnl=bdt_cols_all+[
    "Inv_mupi",
    "Ang_with_scalar_mupi"
]

#these are the additional scalings (to mixing angle squared. ie theta^2 ie mumix^2) produced by Collie which must be applied to find the final exlclued normlastion.

scaling_dict={
    "212":30.83,
    "215":12.10,
    "230":2.83,
    "245":1.36,
    "260":1.18,
    "269":0.85,
    "275":0.82,
    "279":0.99,
    "246":12.91,
    "250":8.57,
    "277":3.05,
    "304":1.46,
    "331":0.85,
    "358":0.54,
    "385":0.92,
    "263p5":3.86,
    "290p5":1.91,
    "317p5":1.18,
    "344p5":0.67,
    "371p5":0.81}



