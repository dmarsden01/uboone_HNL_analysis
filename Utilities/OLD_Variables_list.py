general_var = ['DeltaMed', 'shr_pca_2', 'merge_vtx_x', 'shr_dedx_v_v', 'shr_dedx_u_v', 'topological_score', 
     'pfnplanehits_V', 'merge_vtx_y', 'pt_assume_muon', 'shrclusfrac1', 
     'CosmicDirAll3D', 'DeltaMed1h', 'p', 'shrMCSMom', 'dtrk', 'shr_phi_v', 'slpdg', 'DeltaMed2h', 
     'shr_energy_tot', 'extra_energy_y', 'merge_bestdist', 'secondshower_Y_vtxdist', 'merge_bestdot', 
     'shr_score', 'n_tracks_contained', 'hits_v', 'shrclusdir0', 'pfnplanehits_Y', 
     'CosmicIPAll2DEnds', 'n_tracks', 'shr_pca_0', 'pfnplanehits_U', 'shrclusdir2', 'shr_start_x', 
     'contained_sps_ratio', 'hits_ratio', 'shr_py', 'dvtx', 'CosmicDirAll2DEnds', 'trk_len_v', 
     'trk_chipr_best', 'pt', 'shr_px_v', 'CosmicIP', 'trk_energy_hits_tot', 'shr_pz', 'shr_pz_v', 
     'shr_dedx_U', 'shr_dedx_y_v', 'shr_start_y', 'slclustfrac', 'nu_flashmatch_score', 'trk_energy', 
     'nhits_pl1', 'shr_moliere_rms_v', 'secondshower_Y_eigenratio', 'shr_distance', 'shrclusdir1',
     'shr_px', 'total_hits_y', 'shr_pca_1', 'nslhits_pl1', 'trk_hits_u_tot', 'evnhits', 'flash_time', 
     'secondshower_Y_dir', 'shrPCALen', 'shr_py_v', 'shr_pitch_y_v', 'shr_moliere_avg_v', 'NeutrinoEnergy2',
     'shr_start_x_v', 'trk_mcs_muon_mom_v', 'trk_score_v', 'hits_u', 'shrclusfrac0', 'secondshower_Y_charge',
     'shr_dedx_Y', 'shr_start_y_v', 'shr_dist_v', 'shr_theta', 'shr_phi', 'n_showers', 'secondshower_Y_nhit',
     'trk_energy_muon_mcs', 'shr_start_z', 'pfnhits', 'shrPCA1CAS', 'CosmicIPAll2DOvlp', 
     'shrclusfrac2', 'trk_dir_y_v', 'shr_dedx_V', 'p_assume_muon', 'contained_fraction', 'reco_nu_vtx_y', 
     'trk_dir_z_v', 'shr_pitch_u_v', 'secondshower_Y_dot'] #contained_sps_ratio was duplicated so I took it out.

Cr_variable = ['shrPCA_1Cr', 'shrPCA_1Cr2h', 'shrPCA_3Cr2h', 'shrPCA_3Cr', 'shrPCA_1Cr1h', 'shrPCA_3Cr1h']

cm_cariable = ['shrPCA1CMed_5cm', 'CylFrac_2cm', 'shrStartMCS_2_5cm', 'shrStartMCS_5cm', 'shrStart_2_5cm', 
     'shrStart_5cm', 'shrPCA1CMed_2_5cm', 'CylFrac_1cm', 'shrMCSAS_5cm', 'shrMCSAS_2_5cm']

two_shr_var = ['pi0_dedx1_U', 'pi0_radlen1', 'pi0_shrscore2', 'pi0_dot2', 'pi0_dir1_y', 'pi0_radlen2', 
     'pi0_gammadot', 'pi0_dir2_x', 'pi0_energy2_V', 'pi0_mass_Y', 'pi0_dedx2_fit_V', 'pi0_dedx1_V', 
     'pi0_dedx1_fit_V', 'pi0_shrscore1', 'pi0_dir1_z', 'pi0_dir2_z', 'pi0_dedx2_fit_Y', 'pi0_dedx1_fit_Y', 
     'pi0_dedx2_U', 'pi0_energy1_Y', 'pi0_dedx1_fit_U', 'pi0_dedx2_Y', 'pi0_energy2_Y', 'pi0_energy2_U', 
     'pi0_dot1', 'pi0_dir1_x', 'pi0_dedx2_V', 'pi0_energy1_V', 'pi0_dir2_y', 'pi0_dedx2_fit_U', 
     'pi0_dedx1_Y'] #taken some out of this, which I believe are variables Aditya constructed

weight_related = ['weightSplineTimesTune', 'ppfx_cv', 'npi0']


for_presel = ['nslice', 'trk_score']

event_vars = ['run', 'sub', 'evt']

extras_Owen = ['trk_distance_v', 'trk_theta_v', 'trk_phi_v', 'trk_llr_pid_score_v', 'trk_range_muon_mom_v', 'trk_energy_proton_v', 'trk_calo_energy_y_v', 
               'trk_calo_energy_u_v', 'trk_calo_energy_v_v', 'trk_sce_end_x_v', 'trk_sce_end_y_v', 'trk_sce_end_z_v', 'trk_sce_start_x_v', 'trk_sce_start_y_v', 
               'trk_sce_start_z_v', 'trk_dir_x_v', 'shr_tkfit_dedx_y_v', 'trk_bragg_mip_v', 
               'trk_bragg_mip_v_v', 'pfpplanesubclusters_U', 'pfpplanesubclusters_V', 'pfpplanesubclusters_Y', 
               'reco_nu_vtx_sce_x', 'reco_nu_vtx_sce_y', 'reco_nu_vtx_sce_z', 
               'reco_nu_vtx_x', 'reco_nu_vtx_z', 'SliceCaloEnergy2', 'CosmicIPAll3D', 'slnhits', 
               'bdt_cosmic', 'bdt_ext']
    
Aditya_pre_BDT_vars_ALL_MC = general_var+Cr_variable+cm_cariable+two_shr_var+weight_related+for_presel #For MC
Aditya_pre_BDT_vars_ALL = general_var+Cr_variable+cm_cariable+two_shr_var+for_presel

my_full_vars_MC = general_var+Cr_variable+cm_cariable+two_shr_var+weight_related+for_presel+event_vars+extras_Owen
my_full_vars = general_var+Cr_variable+cm_cariable+two_shr_var+for_presel+event_vars+extras_Owen


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
             #'all_trk_energies']

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

First_pass_vars = event_vars + shower_vars + track_vars + multiplicity_vars + all_other_vars
First_pass_vars_MC = event_vars + shower_vars + track_vars + multiplicity_vars + all_other_vars + weight_related
First_pass_vars_for_BDT = shower_vars + track_vars + multiplicity_vars + all_other_vars


#Without _v variables

# event_vars = ['run', 'sub', 'evt'] #Needed to identify events

# shower_vars = ['shr_energy_tot', 'shr_energy', 'shr_theta', 'shr_pca_0','shr_pca_1','shr_pca_2', 'shr_phi', 'shr_px_v','shr_py_v','shr_pz_v',
#               'shr_openangle', 'shr_tkfit_start_x', 'shr_tkfit_start_y', 'shr_tkfit_start_z','shr_tkfit_theta', 'shr_tkfit_phi',
#               'shr_start_x', 'shr_start_y', 'shr_start_z', 'shr_dedx_Y', 'shr_dedx_V','shr_dedx_U',  
#               'shr_tkfit_dedx_Y','shr_tkfit_dedx_V', 'shr_tkfit_dedx_U', 'shrmoliereavg', 'shr_distance', 'shr_score', 'shr_hits_max']

# track_vars = ['trk_len', 'trk_theta', 'trk_phi', 'trk_energy', 'trk_energy_tot', 'trk_distance', 'trk_score', 'trk_chipr_best',
#              'trk_bragg_p', 'trk_bragg_mip', 'trk_hits_max', 'trk_start_x_v', 'trk_start_y_v', 'trk_start_z_v', 
#              'trk_end_x_v', 'trk_end_y_v', 'trk_end_z_v', 'trk_dir_x_v', 'trk_dir_y_v', 'trk_dir_z_v']

# multiplicity_vars = ['nslice', 'n_pfps', 'n_tracks', 'n_showers']

# all_other_vars = ['merge_bestdot', 'merge_bestdist', 'merge_vtx_x', 'merge_vtx_y', 'merge_vtx_z', 'shrclusdir0','shrclusdir1','shrclusdir2',
#                  'shrclusfrac0','shrclusfrac1','shrclusfrac2', 'reco_nu_vtx_x', 'reco_nu_vtx_y', 'reco_nu_vtx_z',
#                  'crtveto', 'topological_score', 'flash_time', 'nu_flashmatch_score', 'shrPCA_1Cr', 'shrPCA_2Cr','shrPCA_3Cr',
#                  'shrMCSMom', 'shrStart_5cm', 'DeltaMed', 'CylFrac_1cm', 'NeutrinoEnergy2', 'SliceCaloEnergy2', 'pi0_radlen1','pi0_radlen2', 
#                  'pi0_dot1', 'pi0_dot2', 'pi0_shrscore1', 'pi0_shrscore2', 'pi0_dir1_x', 'pi0_dir1_y', 'pi0_dir1_z', 'pi0_dir2_x', 'pi0_dir2_y',
#                  'pi0_dir2_z', 'pi0_shrscore1', 'pi0_shrscore2', 'pi0_gammadot', 'secondshower_Y_charge', 'secondshower_Y_vtxdist', 
#                  'secondshower_Y_eigenratio', 'secondshower_Y_dot', 'secondshower_Y_dir', 'bdt_cosmic', 'bdt_ext', 'contained_fraction', 
#                  'contained_sps_ratio', 'extra_energy_y', 'slclustfrac', 'slnhits', 'pt','p_assume_muon']

# weight_related = ['weightSplineTimesTune', 'ppfx_cv', 'npi0'] #For MC only, includes true npi0 for scaling

# general_var = [ 
#      'pfnplanehits_V', 'pt_assume_muon', 
#      'CosmicDirAll3D', 'DeltaMed1h', 'p', 'shrMCSMom', 'dtrk', 'shr_phi_v', 'slpdg', 'DeltaMed2h', 
#      'n_tracks_contained', 'hits_v', 'pfnplanehits_Y', 
#      'CosmicIPAll2DEnds', 'pfnplanehits_U', 
#      'hits_ratio', 'dvtx', 'CosmicDirAll2DEnds', 'trk_len_v', 
#      'CosmicIP', 'trk_energy_hits_tot',  
#      'nhits_pl1', 'shr_moliere_rms_v', 
#      'total_hits_y', 'nslhits_pl1', 'trk_hits_u_tot', 'evnhits',
#      'shrPCALen', 'shr_pitch_y_v', 'shr_moliere_avg_v', 
#      'trk_mcs_muon_mom_v', 'trk_score_v', 'hits_u', 
#      'shr_dist_v', 'n_showers', 'secondshower_Y_nhit',
#      'trk_energy_muon_mcs', 'pfnhits', 'shrPCA1CAS', 'CosmicIPAll2DOvlp', 
#      'trk_dir_y_v',
#      'trk_dir_z_v', 'shr_pitch_u_v']

# Cr_variable = ['shrPCA_1Cr', 'shrPCA_1Cr2h', 'shrPCA_3Cr2h', 'shrPCA_3Cr', 'shrPCA_1Cr1h', 'shrPCA_3Cr1h']

# cm_cariable = ['shrPCA1CMed_5cm', 'CylFrac_2cm', 'shrStartMCS_2_5cm', 'shrStartMCS_5cm', 'shrStart_2_5cm', 
#      'shrStart_5cm', 'shrPCA1CMed_2_5cm', 'CylFrac_1cm', 'shrMCSAS_5cm', 'shrMCSAS_2_5cm']

# two_shr_var = ['pi0_dedx1_U',  'pi0_dot2', 
#      'pi0_energy2_V', 'pi0_mass_Y', 'pi0_dedx2_fit_V', 'pi0_dedx1_V', 
#      'pi0_dedx1_fit_V',  'pi0_dedx2_fit_Y', 'pi0_dedx1_fit_Y', 
#      'pi0_dedx2_U', 'pi0_energy1_Y', 'pi0_dedx1_fit_U', 'pi0_dedx2_Y', 'pi0_energy2_Y', 'pi0_energy2_U', 
#      'pi0_dot1', 'pi0_dedx2_V', 'pi0_energy1_V', 'pi0_dedx2_fit_U', 
#      'pi0_dedx1_Y'] #taken some out of this, which I believe are variables Aditya constructed

# weight_related = ['weightSplineTimesTune', 'ppfx_cv', 'npi0']


# for_presel = ['nslice', 'trk_score', 'n_pfps']

# extras_Owen = ['trk_distance_v', 'trk_theta_v', 'trk_phi_v', 'trk_llr_pid_score_v', 'trk_range_muon_mom_v', 'trk_energy_proton_v', 'trk_calo_energy_y_v', 
#                'trk_calo_energy_u_v', 'trk_calo_energy_v_v', 'trk_sce_end_x_v', 'trk_sce_end_y_v', 'trk_sce_end_z_v', 'trk_sce_start_x_v', 'trk_sce_start_y_v', 
#                'trk_sce_start_z_v', 'trk_dir_x_v', 'shr_tkfit_dedx_y_v', 'trk_bragg_mip_v', 
#                'trk_bragg_mip_v_v', 'pfpplanesubclusters_U', 'pfpplanesubclusters_V', 'pfpplanesubclusters_Y', 
#                'reco_nu_vtx_sce_x', 'reco_nu_vtx_sce_y', 'reco_nu_vtx_sce_z', 
#                'reco_nu_vtx_x', 'reco_nu_vtx_z', 'SliceCaloEnergy2', 'CosmicIPAll3D', 'slnhits', 
#                'bdt_cosmic', 'bdt_ext']
