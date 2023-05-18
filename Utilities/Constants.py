sample_type = {"overlay":"MC",
              "dirtoverlay":"MC",
              "beamoff":"data",
              "signal":"MC_signal",
              "beamgood":"data"}

root_dir = 'nuselection'
main_tree = 'NeutrinoSelectionFilter'

average_HNL_direction = [0.32, 0.74, -0.59] #Calculated from the weighted files "truth" HNL information

#------------------------------------------#
#------Pre-selection related constants-----#
#------------------------------------------#

#From Owen
max_x_cut=253
min_x_cut=9
max_y_cut=112
min_y_cut=-112
max_z_cut=1020
min_z_cut=14

##TPC
TPCxlo=0
TPCxhi=256
TPCylo=-115
TPCyhi=117
TPCzlo=0
TPCzhi=1037

max_extent_cut=f"min_y>{min_y_cut} and max_y<{max_y_cut} and min_z>{min_z_cut} and max_z<{max_z_cut}" + " and "+f"min_x>{min_x_cut} and max_x<{max_x_cut}"

Preselection_dict = {"nslice":"nslice==1",
"flash_time":"flash_time > 6.55 and flash_time < 16.5",
"nu_flashmatch_score":"nu_flashmatch_score < 15",#may need to take out.
"NeutrinoEnergy2":"NeutrinoEnergy2 < 500",
"Fiducial_cut":max_extent_cut,
"contained_fraction":"contained_fraction > 0.9",
"trk_score_v":"trk_score_v < 0.97",
"n_pfps":"n_pfps < 6"
} #Should I add topological score cut?


Preselection_dict_pi0 = {"nslice":"nslice==1",
"flash_time":"flash_time > 6.55 and flash_time < 16.5",
"nu_flashmatch_score":"nu_flashmatch_score < 15",
"NeutrinoEnergy2":"NeutrinoEnergy2 < 400",
"contained_fraction":"contained_fraction > 0.9",
"trk_score_v":"trk_score_v < 0.97",
"n_pfps":"n_pfps < 6",
"topological_score":"topological_score < 0.95",
"shr_tkfit_phi_v":"shr_tkfit_phi_v < -1.0 or shr_tkfit_phi_v > 0.5"
}

Preselection_dict_for_plot = {"nslice":"nslice==1",
"flash_time":"flash_time > 6.55 and flash_time < 16.5",
"nu_flashmatch_score":"nu_flashmatch_score < 15",
"NeutrinoEnergy2":"NeutrinoEnergy2 < 500",
"contained_fraction":"contained_fraction > 0.9",
"trk_score":"trk_score < 0.97",
"n_pfps":"n_pfps < 6"
}

Preselection_dict_crtveto = {"crtveto":"crtveto==0"}


#------------------------------------#
#------Scaling related constants-----#
#------------------------------------#

HNL_mass_samples = [2, 10, 20, 50, 100, 150] #These are the current mass points for decays to e+e-, cautious of 2MeV sample sclaing
HNL_ee_samples_names = ["2_ee", "10_ee", "20_ee", "50_ee", "100_ee", "150_ee"]
HNL_mass_pi0_samples = [150, 180, 200, 220, 240, 245]
HNL_mass_pi0_samples_names = ["150_pi0", "180_pi0", "200_pi0", "220_pi0", "240_pi0", "245_pi0"]
HNL_ee_dirac_mass_samples = [10, 100, 150]
HNL_pi0_dirac_mass_samples = [150, 200, 245]
HNL_ee_dirac_names = ["10_ee_dirac", "100_ee_dirac", "150_ee_dirac"]
HNL_pi0_dirac_names = ["150_pi0_dirac", "200_pi0_dirac", "245_pi0_dirac"]
Run1_ee_DetVar_samples = [50, 100, 150]
Run3_ee_DetVar_samples = [2, 10, 20, 50, 100]
Run1_pi0_DetVar_samples = []
Run3_pi0_DetVar_samples = [150, 180, 220, 240,245] #200 is broken for some reason

theta_mu_4 = 1e-4 #This is the same for *most* of the samples and was set in the production .fcls, now made a new dict because some are higher than this
theta_mu_4_dict = {2:1e-1,
                   10:1e-2,
                   20:1e-4,
                   50:1e-4,
                   100:1e-4,
                   150:1e-4,
                   180:1e-4,
                   200:1e-4,
                   220:1e-4,
                   240:1e-4,
                   245:1e-4}

pi0_scaling_factor = 0.759 #For standard NuMI overlay files

Run1_POT = 2e20 #Taken from the NuMI samples page
Run3_POT = 5.0e20 #Taken from the NuMI samples page
OnBeam_EA9CNT_wcut_run1 = 5268051.0 #"Triggers" taken from the NuMI samples page
OnBeam_EA9CNT_wcut_run3 = 10363728.0 #"Triggers" taken from the NuMI samples page
BeamOff_scaling_for_nus = 0.98 #This Factor described by Owen: Look in script 1_POT...
DIRT_run1_scaling = 0.75 #NOT SURE where this comes from, apparently it is standard procedure for NuMI DIRT
DIRT_run3_scaling = 0.35 #NOT SURE where this comes from, apparently it is standard procedure for NuMI DIRT
NuMI_KDAR_scaling_run1 = 8.0 #This comes from the discrepancy between numu flux from KDAR dump between Geant4 and MiniBooNE measurement. Taken from Owen's thesis
NuMI_KDAR_scaling_run3 = 8.6

run1_event_numbers = {2: 36059,
                     10: 35032,
                     20: 36815,
                     50: 35518, 
                     100: 36881, 
                     150: 35254, #Fixed for new sample
                     "2_ee": 36059,
                     "10_ee": 35032,
                     "20_ee": 36815,
                     "50_ee": 35518, 
                     "100_ee": 36881, 
                     "150_ee": 35254,
                     180: 36670, 
                     200: 35325,
                     220: 36849,
                     240: 35400,
                     245: 36877,
                     "150_pi0": 30261,
                     "180_pi0": 36452,
                     "200_pi0": 35238,
                     "220_pi0": 36868,
                     "240_pi0": 35264,
                     "245_pi0": 36698,
                     '10_ee_dirac':34373,
                     '100_ee_dirac':34114,
                     "150_ee_dirac":31355,
                     "150_pi0_dirac":35110,
                     "200_pi0_dirac":33397,
                     "245_pi0_dirac":34878,
                     "overlay": 914729,
                     "dirtoverlay": 569506,
                     "beamoff": 904362,
                     "beamgood": 610496}

run3_event_numbers = {2: 45159,
                     10: 44463,
                     20: 46022,
                     50: 44579, 
                     100: 45304, 
                     150: 44031, #Fixed for new sample
                     "2_ee": 45159,
                     "10_ee": 44463,
                     "20_ee": 46022,
                     "50_ee": 44579, 
                     "100_ee": 45304, 
                     "150_ee": 44031,
                     180: 46038, 
                     200: 45003,
                     220: 46093,
                     240: 44894,
                     245: 46114,
                     "150_pi0": 45827,
                     "180_pi0": 44811,
                     "200_pi0": 43943,
                     "220_pi0": 44503,
                     "240_pi0": 45704,
                     "245_pi0": 44073,
                     "overlay": 748702,
                     "dirtoverlay": 389264,
                     "beamoff": 3211097,
                     "beamgood": 1104349}  

run1_sum_weights = {"overlay": 867669.0,"dirtoverlay": 535760.1}
run3_sum_weights = {"overlay": 704253.8,"dirtoverlay": 363852.9}

run1_POT_scaling_dict = {2: 7.116758217611076e-08, #All ee samples (up to 150MeV) now corrected for sampling fix
                         10: 4.421195517726753e-07,
                         20: 2.6922986790110916e-13,
                         50: 8.776574871312947e-11, 
                         100: 8.380344432250804e-09, 
                         150: 1.594402204319872e-07, #Fixed to be corret scaling with sampling fix 
                         "2_ee": 7.116758217611076e-08, #All ee samples (up to 150MeV) now corrected for sampling fix
                         "10_ee": 4.421195517726753e-07,
                         "20_ee": 2.6922986790110916e-13,
                         "50_ee": 8.776574871312947e-11, 
                         "100_ee": 8.380344432250804e-09, 
                         "150_ee": 1.594402204319872e-07,
                         180: 3.373080077973945e-05, #Wrong, but not using 180-245MeV ee samples anymore
                         200: 9.293586743019126e-05,
                         220: 0.00020120909708414096,
                         240: 0.0004102200024715987,
                         245: 0.0004631133119186411,
                         "150_pi0": 2.348836022876738e-06,
                         "180_pi0": 2.773198013287162e-05,
                         "200_pi0": 8.011550298626092e-05,
                         "220_pi0": 0.0001701935111423321,
                         "240_pi0": 0.00034948084037732583,
                         "245_pi0": 0.00039327366427535133,
                         "10_ee_dirac":2.259966244579603e-07,
                         "100_ee_dirac":4.51710913729398e-09,
                         "150_ee_dirac":8.919576191863471e-08,
                         "150_pi0_dirac":1.0222965933561835e-06,
                         "200_pi0_dirac":4.233726447611382e-05,
                         "245_pi0_dirac":0.00020749021228398003,
                         "overlay": 0.08559729531465318,
                         "dirtoverlay": 0.08961042442406035,
                         "beamoff": 0.5612087579382191,
                         "beamgood": 1.0}

run3_POT_scaling_dict = {2: 1.4264872844261716e-07, #All ee samples (up to 150MeV) now corrected for sampling fix
                         10: 8.612201693171446e-07,
                         20: 5.457161021952394e-13, 
                         50: 1.784272446212085e-10, 
                         100: 1.709941880563678e-08, 
                         150: 3.179231163470327e-07, #Fixed to be corret scaling with sampling fix
                         "2_ee": 1.4264872844261716e-07, #All ee samples (up to 150MeV) now corrected for sampling fix
                         "10_ee": 8.612201693171446e-07,
                         "20_ee": 5.457161021952394e-13, 
                         "50_ee": 1.784272446212085e-10, 
                         "100_ee": 1.709941880563678e-08, 
                         "150_ee": 3.179231163470327e-07,
                         180: 6.790004454782067e-05, #Wrong, but not using 180-245MeV ee samples anymore 
                         200: 0.0001841720407411439,
                         220: 0.00040388981532410006,
                         240: 0.0008189646415804218,
                         245: 0.000833548258426593,
                         "150_pi0": 3.894780273065507e-06,
                         "180_pi0": 5.702796319652067e-05,
                         "200_pi0": 0.00016293389328773464,
                         "220_pi0": 0.00035905258177289234,
                         "240_pi0": 0.0006812634124592591,
                         "245_pi0": 0.0009711075263359489,
                         "overlay": 0.2513368817255014,
                         "dirtoverlay": 0.16953052634982632,
                         "beamoff": 0.3089104916683624,
                         "beamgood": 1.0}

run1_overlay_detvar_POT = {"WireModX":3.67e20, #Just taken from spreadsheet, but matches the POT saved in the Trees
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

run3_overlay_detvar_POT = {"WireModX":3.24e20, #Taken from POT counter
                           "WireModYZ":3.36e20,
                           "WireModThetaXZ":3.2e20,
                           "WireModThetaYZ":3.36e20,
                           "WireModdEdX":3.17e20,
                           "LYDown":2.81e20,
                           "LYRayleigh":2.81e20,
                           "LYAttenuation":3.31e20,
                           "SCE":3.33e20,
                           "Recomb2":3.3e20,
                           "CV":3.72e20}

Multisim_univs = {"weightsPPFX":600,
                  "weightsGenie":600,
                  "weightsReint":1000}

Detector_variations = ["WireModX", "WireModYZ", "WireModThetaXZ", "WireModThetaYZ", "WireModdEdX",
                       "LYDown", "LYRayleigh", "LYAttenuation", "SCE", "Recomb2", "CV"]

run1_pi0_file_loc_dict = {"150_pi0": '../NuMI_signal/KDAR_dump/sfnues/pi0/sfnues_KDAR_dump_150_pi0_Umu4_majorana_FHC.root',
                          "180_pi0": '../NuMI_signal/KDAR_dump/sfnues/pi0/sfnues_KDAR_dump_180_pi0_Umu4_majorana_FHC.root',
                          "200_pi0": '../NuMI_signal/KDAR_dump/sfnues/pi0/sfnues_KDAR_dump_200_pi0_Umu4_majorana_FHC.root',
                          "220_pi0": '../NuMI_signal/KDAR_dump/sfnues/pi0/sfnues_KDAR_dump_220_pi0_Umu4_majorana_FHC.root',
                          "240_pi0": '../NuMI_signal/KDAR_dump/sfnues/pi0/sfnues_KDAR_dump_240_pi0_Umu4_majorana_FHC.root',
                          "245_pi0": '../NuMI_signal/KDAR_dump/sfnues/pi0/sfnues_KDAR_dump_245_pi0_Umu4_majorana_FHC.root'}

#--------------------------------------#
#------Plotting related constants------#
#--------------------------------------#

sample_colours = {'overlay':'mediumblue',
                  'dirtoverlay':'cornflowerblue',
                  'beamoff':'limegreen',
                  'signal':'red'}

presel_var_names = {'nslice':"Neutrino slice", 'flash_time':"Flash time", 'nu_flashmatch_score':"Flashmatch score", 
                    'NeutrinoEnergy2':"Energy in slice", 'Fiducial_cut':"Fiducial volume", 
                    'contained_fraction':"Containted fraction", 'trk_score_v':"Track score",'trk_score':"Track score",
                    'n_pfps':"Object multiplicity", 'crtveto':"CRT veto"}

variable_names_dict = {'nslice':"Neutrino slice", #Should contain all names for possible variables to be plotted
                       'flash_time': r"Flash time [$\mu$s]",
                       'nu_flashmatch_score': "Flashmatch score",
                       'NeutrinoEnergy2':"Energy in slice [MeV]",
                       'contained_fraction':"Containted fraction",
                       'trk_score':"Track score",
                       'trk_score_v':"Track score",
                       'n_pfps':"Object multiplicity",
                       'crtveto':"CRT veto",
                       'shrclusdir2':"Shower cluster direction [degrees]",
                       'n_tracks':"Reconstructed tracks",
                       'n_showers':"Reconstructed showers",
                       'trk_energy':"Highest reconstructed track energy [GeV]",
                       'shr_theta_v':"Shower theta [radians]",
                       'contained_sps_ratio':"Contained space points fraction",
                       'shr_px_v':"Shower x-momentum fraction",
                       'trk_end_x_v':"Track end x [cm]",
                       'pfnplanehits_V':"V plane hits",
                       'pfnplanehits_U':"U plane hits",
                       'pfnplanehits_Y':"Y plane hits",
                       'shr_phi_v':r"Shower $\phi$ [radians]",
                       'shr_pz_v':"Shower z-momentum fraction",
                       'trk_theta_v':r"Track $\theta$ [radians]",
                       'trk_phi_v':r"Track $\phi$ [radians]",
                       'trk_energy_hits_tot':"Total track energy from hits [GeV]",
                       'trk_energy_tot':"Total track energy [GeV]",
                       'trk_dir_z_v':"Track z-momentum fraction",
                       'SliceCaloEnergy2':"Calorimatric energy in slice [MeV]",
                       'shr_energy_tot':"Total shower energy [GeV]",
                       'trk_chipr_best':'Proton chi score',
                       'trk_calo_energy_u_v':'Track energies U plane [MeV]',
                       "Fiducial_cut": "Fiducial volume",
                       "min_x": "Minimum extent in x [cm]","max_x": "Maximum extent in x [cm]",
                       "min_y": "Minimum extent in y [cm]","max_y": "Maximum extent in y [cm]",
                       "min_z": "Minimum extent in z [cm]","max_z": "Maximum extent in z [cm]"}

limit_locs = {"SIN":'limit_files/SIN_BR_ratio.csv',
              "PIENU":'limit_files/PIENU_2019.csv',
              "PS191":'limit_files/PS191_1988.csv',
              "KEK":'limit_files/KEK_1984_E89_plus_E104.csv',
              "KEK_1982":'limit_files/KEK_1982.csv',
              "E949":'limit_files/E949_2015.csv'}

limit_colours = {"SIN":"C17",
                 "PIENU":"C1",
                 "PS191":"C2",
                 "KEK":"C3", #1984 "combined" limit
                 "KEK_1982":"C3",
                 "E949":"C4"}




#----------------------------------------------------------#
#------OLD constants, for depricated code and samples------#
#----------------------------------------------------------#

sampling_fix_scalings = {2: 0.0313, #Sceptical of 2MeV POT counting because I used such a large mixing angle. Not sure if correct
                         10: 0.0987, #This dict is only here for reference, it is now folded into the new "POT_scaling" dicts 
                         20: 0.1042,
                         50: 0.1060, 
                         100: 0.1043,
                         150: 0.0452,
                         "150_pi0": 0.5646,
                         "180_pi0": 0.8187,
                         "200_pi0": 0.8541,
                         "220_pi0": 0.8544,
                         "240_pi0": 0.8484,
                         "245_pi0": 0.8583}

OLD_HNL_mass_samples = [2, 10, 20, 50, 100, 150, 180, 200, 220, 240, 245]

OLD_Preselection_dict = {"nslice":"nslice==1",
"flash_time":"flash_time > 6.55 and flash_time < 16.5",
"nu_flashmatch_score":"nu_flashmatch_score < 15",#may need to take out.
"NeutrinoEnergy2":"NeutrinoEnergy2 < 500",
"contained_fraction":"contained_fraction > 0.9",
"trk_score_v":"trk_score_v < 0.97",
"n_pfps":"n_pfps < 6"
}#want to add topological score, maybe a cut to keep below 0.96.

SF_overlay_run1 = 0.08559729531465318
SF_dirt_run1 = 0.08961042442406035
SF_EXT_run1 = 0.5612087579382191

SF_overlay_run3 = 0.2513368817255014
SF_dirt_run3 = 0.16953052634982632
SF_EXT_run3 = 0.3089104916683624

event_number_overlay_run1 = 914729
event_number_dirt_run1 = 569506
event_number_EXT_run1 = 904362

event_number_overlay_run3 = 748702
event_number_dirt_run3 = 389264
event_number_EXT_run3 = 3211097

OLD_run1_POT_scaling_dict = {2: 2.273110163825259e-06, #Before the sampling fix, which significantly reduces event rate
                         10: 4.478610916645953e-06,
                         20: 2.583047654702099e-12,
                         50: 8.280121178469969e-10, 
                         100: 8.031518427782856e-08, 
                         150: 3.5533992790484712e-06, 
                         180: 3.373080077973945e-05, 
                         200: 9.293586743019126e-05,
                         220: 0.00020120909708414096,
                         240: 0.0004102200024715987,
                         245: 0.0004631133119186411,
                         "150_pi0": 4.159838350914407e-06,
                         "180_pi0": 3.38744768337177e-05,
                         "200_pi0": 9.379957491349839e-05,
                         "220_pi0": 0.00019919442304235102,
                         "240_pi0": 0.00041190884760711096,
                         "245_pi0": 0.000458175050366498,
                         "overlay": 0.08559729531465318,
                         "dirtoverlay": 0.08961042442406035,
                         "beamoff": 0.5612087579382191,
                         "beamgood": 1.0}

OLD_run3_POT_scaling_dict = {2: 4.55623564219535e-06,
                         10: 8.724043160892894e-06,
                         20: 5.235714406049288e-12, 
                         50: 1.6833437060203603e-09, 
                         100: 1.638766739864928e-07, 
                         150: 7.045112441597867e-06, 
                         180: 6.790004454782067e-05, 
                         200: 0.0001841720407411439,
                         220: 0.00040388981532410006,
                         240: 0.0008189646415804218,
                         245: 0.0009335050190166368,
                         "150_pi0": 6.897738365081695e-06,
                         "180_pi0": 6.965937552669846e-05,
                         "200_pi0": 0.000190763701901887,
                         "220_pi0": 0.0004202350100663058,
                         "240_pi0": 0.0008029579728605566,
                         "245_pi0": 0.0009711075263359489,
                         "overlay": 0.2513368817255014,
                         "dirtoverlay": 0.16953052634982632,
                         "beamoff": 0.3089104916683624,
                         "beamgood": 1.0}

OLDEST_run1_POT_scaling_dict = {10:4.478610916645953e-06, #Before new integration result incorporated.
                         20: 2.1387912933253144e-12,
                         50: 6.590170245207509e-10, 
                         100: 6.913852703765691e-08, 
                         150: 3.301953636715338e-06, 
                         180: 3.479531678617752e-05, 
                         200: 9.077313301681812e-05,
                         "150_pi0": 4.159838350914407e-06}

OLDEST_run3_POT_scaling_dict = {20: 1.7657423469953374e-12, 
                         50: 5.448786613287057e-10, 
                         100: 5.5322831399190343e-08, 
                         150: 2.623088564712337e-06, 
                         180: 2.624339736997261e-05, 
                         200: 7.283055184784342e-05}

OLDER_generator_mass_points = [20, 50, 100, 150, 180, 200] #These are the samples created with the OLD generator (both run1 and run3) 
    
Old_gen_HNL_scalings = {20:1228.625, #Scalings Pawel sent to use for the HNL samples created with the old generator version
                        50:199.4745, #These will not direclty work because the different BR will also affect the number which decay away before reaching uboone
                        100:49.9434,
                        150:22.19445,
                        180:15.43885,
                        200:12.49565}