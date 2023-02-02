#My variables, used the _v versions (i.e should be all tracks/showers included)

event_vars = ['run', 'sub', 'evt'] #Needed to identify events

shower_vars = ['shr_energy_tot', 'shr_energy', 'shr_theta_v', 'shr_pca_0','shr_pca_1','shr_pca_2', 'shr_phi_v', 'shr_px_v','shr_py_v','shr_pz_v',
              'shr_openangle_v', 'shr_tkfit_start_x_v', 'shr_tkfit_start_y_v', 'shr_tkfit_start_z_v','shr_tkfit_theta_v', 'shr_tkfit_phi_v',
              'shr_start_x_v', 'shr_start_y_v', 'shr_start_z_v', 'shr_dedx_y_v', 'shr_dedx_v_v','shr_dedx_u_v',  
              'shr_tkfit_dedx_y_v','shr_tkfit_dedx_v_v', 'shr_tkfit_dedx_u_v', 'shrmoliereavg', 'shr_distance', 'shr_score', 'shr_hits_max']
              #'all_shr_energies']

track_vars = ['trk_len_v', 'trk_theta_v', 'trk_phi_v', 'trk_energy', 'trk_energy_tot', 'trk_energy_hits_tot', 'trk_distance_v', 'trk_score_v',
             'trk_chipr_best','trk_bragg_p_v', 'trk_bragg_mip_v', 'trk_hits_max', 'trk_start_x_v', 'trk_start_y_v', 'trk_start_z_v', 
             'trk_end_x_v', 'trk_end_y_v', 'trk_end_z_v', 'trk_dir_x_v', 'trk_dir_y_v', 'trk_dir_z_v', 
             'trk_calo_energy_y_v','trk_calo_energy_u_v','trk_calo_energy_v_v']
             #'all_trk_energies'] #To include this I think I need to resize the dataframe somehow so it can be flattened

multiplicity_vars = ['nslice', 'n_pfps', 'n_tracks', 'n_showers']

all_other_vars = ['merge_bestdot', 'merge_bestdist', 'merge_vtx_x', 'merge_vtx_y', 'merge_vtx_z', 'shrclusdir0','shrclusdir1','shrclusdir2',
                 'shrclusfrac0','shrclusfrac1','shrclusfrac2', 'reco_nu_vtx_x', 'reco_nu_vtx_y', 'reco_nu_vtx_z',
                 'crtveto', 'topological_score', 'flash_time', 'nu_flashmatch_score', 'shrPCA_1Cr', 'shrPCA_2Cr','shrPCA_3Cr',
                 'shrMCSMom', 'shrStart_5cm', 'DeltaMed', 'CylFrac_1cm', 'NeutrinoEnergy2', 'SliceCaloEnergy2', 'pi0_radlen1','pi0_radlen2', 
                 'pi0_dot1', 'pi0_dot2', 'pi0_shrscore1', 'pi0_shrscore2', 'pi0_dir1_x', 'pi0_dir1_y', 'pi0_dir1_z', 'pi0_dir2_x', 'pi0_dir2_y',
                 'pi0_dir2_z', 'pi0_gammadot', 'secondshower_Y_charge', 'secondshower_Y_vtxdist', 
                 'secondshower_Y_eigenratio', 'secondshower_Y_dot', 'secondshower_Y_dir', 'bdt_cosmic', 'bdt_ext', 'contained_fraction', 
                 'contained_sps_ratio', 'extra_energy_y', 'slclustfrac', 'slnhits', 'pt','p_assume_muon', 'CosmicIPAll3D',
                 'pfnplanehits_U', 'pfnplanehits_V','pfnplanehits_Y']

weight_related = ['weightSplineTimesTune', 'ppfx_cv', 'npi0'] #For MC only, includes true npi0 for scaling

sys_vars = ["weightsPPFX","weightsGenie","weightsReint"]

swtrig_vars = ["swtrig_pre","swtrig_post"] #Think this is only for flash time plots

Fiducial_variables = ["trk_sce_start_x_v","trk_sce_start_y_v","trk_sce_start_z_v","trk_sce_end_x_v","trk_sce_end_y_v","trk_sce_end_z_v" ]

First_pass_vars = event_vars + shower_vars + track_vars + multiplicity_vars + all_other_vars
First_pass_vars_MC = event_vars + shower_vars + track_vars + multiplicity_vars + all_other_vars + weight_related
First_pass_vars_for_BDT = shower_vars + track_vars + multiplicity_vars + all_other_vars

#Reduced variables for BDT

New_feature_list = ['shrclusdir2', 'n_tracks', 'trk_energy', 'shr_theta_v', 'contained_sps_ratio', 'trk_chipr_best', 'shr_px_v',
                    'trk_end_x_v', 'n_pfps', 'pfnplanehits_V', 'pfnplanehits_U', 'trk_calo_energy_u_v', 'nu_flashmatch_score', 'trk_score_v',
                    'NeutrinoEnergy2', 'shr_phi_v', 'pfnplanehits_Y', 'shr_pz_v', 'trk_theta_v', 'trk_phi_v', 'trk_energy_hits_tot',
                    'trk_dir_z_v', 'SliceCaloEnergy2'] #Common most important features across all BDT models

Rest_for_preselection = ['nslice', 'flash_time', 'contained_fraction', 'trk_score', 'crtveto']

Additional_vars = ['shr_energy_tot','trk_energy_tot', 'n_showers']

New_variables = event_vars + New_feature_list + Rest_for_preselection + Additional_vars + swtrig_vars
New_variables_MC = New_variables + weight_related


Final_features = ['shr_theta_v', 'shr_phi_v', 'shr_px_v', 'shr_py_v', 'shr_pz_v', 'shrclusdir0', 'shrclusdir1', 'shrclusdir2', 'shr_energy_tot', #shr related
                  'trk_theta_v', 'trk_phi_v', 'trk_dir_x_v', 'trk_dir_y_v', 'trk_dir_z_v', 'trk_energy', 'trk_energy_hits_tot', 'trk_energy_tot', #track related
                  'trk_score_v', 'trk_calo_energy_u_v', 'trk_end_x_v', 'trk_chipr_best', #rest track related
                  'pfnplanehits_U', 'pfnplanehits_V', 'pfnplanehits_Y', 'NeutrinoEnergy2','SliceCaloEnergy2', 'nu_flashmatch_score', 'contained_sps_ratio', #energy/misc
                  'flash_time', 'contained_fraction', 'trk_score', 'crtveto'] #rest for preselection


#-------------FINAL VARIABLES------------------------#
Final_variable_list = event_vars + multiplicity_vars + swtrig_vars + Fiducial_variables + Final_features
Final_variable_list_MC = Final_variable_list + weight_related

#Variables used in Pre-selection
Preselection_vars = ["nslice","flash_time","nu_flashmatch_score","NeutrinoEnergy2","contained_fraction","trk_score","trk_score_v","n_pfps"] + swtrig_vars
Preselection_vars_MC = Preselection_vars + weight_related
Preselection_vars_CRT = Preselection_vars+["crtveto"]
Preselection_vars_CRT_MC = Preselection_vars_CRT+weight_related

Truth_vars = ['mc_pdg','mc_E','mc_vx','mc_vy','mc_vz','mc_endx','mc_endy','mc_endz','mc_px','mc_py','mc_pz',
              'mc_primary_pdg','mc_primary_px','mc_primary_py','mc_primary_pz','mc_completeness','mc_purity','npi0']