sample_type = {"overlay":"MC",
              "dirtoverlay":"MC",
              "beamoff":"data",
              "signal":"MC_signal"}

root_dir = 'nuselection'
main_tree = 'NeutrinoSelectionFilter'

Preselection_dict = {"nslice":"nslice==1",
"flash_time":"flash_time > 6.55 and flash_time < 16.5",
"nu_flashmatch_score":"nu_flashmatch_score < 15",#may need to take out.
"NeutrinoEnergy2":"NeutrinoEnergy2 < 500",
"contained_fraction":"contained_fraction > 0.9",
"trk_score_v":"trk_score_v < 0.97",
"n_pfps":"n_pfps < 6"
}#want to add topological score, maybe a cut to keep below 0.96

Preselection_dict_for_plot = {"nslice":"nslice==1",
"flash_time":"flash_time > 6.55 and flash_time < 16.5",
"nu_flashmatch_score":"nu_flashmatch_score < 15",#may need to take out.
"NeutrinoEnergy2":"NeutrinoEnergy2 < 500",
"contained_fraction":"contained_fraction > 0.9",
"trk_score":"trk_score < 0.97",
"n_pfps":"n_pfps < 6"
}

Preselection_dict_crtveto = {"crtveto":"crtveto==0"}

HNL_mass_samples = [20, 50, 100, 150, 180, 200] #This are the current mass points for decays to e+e-
HNL_mass_pi0_samples = [150]

theta_mu_4 = 1e-4 #This is the same for all samples and was set in the production .fcls

pi0_scaling_factor = 0.759

run1_POT_scaling_dict = {20: 2.1387912933253144e-12, 
                         50: 6.590170245207509e-10, 
                         100: 6.913852703765691e-08, 
                         150: 3.301953636715338e-06, 
                         180: 3.479531678617752e-05, 
                         200: 9.077313301681812e-05,
                         "150_pi0": 4.159838350914407e-06}

SF_overlay_run1 = 0.08559729531465318
SF_dirt_run1 = 0.08961042442406035
SF_EXT_run1 = 0.5612087579382191
SF_signal_run1_100MeV = 6.913852627300102e-08

run1_overlay_detvar_POT = {"WireModX":3.67e20, #Should check these with POT script, just taken from spreadsheet atm
                           "WireModYZ":3.7e20,
                           "WireModThetaXZ":3.68e20,
                           "WireModThetaYZ":3.69e20,
                           "WireModdEdX":3.67e20,
                           "LYDown":3.57e20,
                           "LYRayleigh":3.69e20,
                           "LYAttenuation":3.64e20,
                           "SCE":3.71e20,
                           "Recomb2":3.72e20,
                           "CV":3.69e20}
                           
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

Multisim_univs = {"weightsPPFX":600,
                  "weightsGenie":600,
                  "weightsReint":1000}

Detector_variations = ["WireModX", "WireModYZ", "WireModThetaXZ", "WireModThetaYZ", "WireModdEdX",
                       "LYDown", "LYRayleigh", "LYAttenuation", "SCE", "Recomb2", "CV"]

##TPC
TPCxlo=0
TPCxhi=256
TPCylo=-115
TPCyhi=117
TPCzlo=0
TPCzhi=1037

