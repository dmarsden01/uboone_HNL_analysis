import numpy as np

bins_var = {'shrclusdir2':np.linspace(0,360, 21),
            'n_tracks':np.linspace(-0.5,5.5, 7),
            'n_pfps':np.linspace(0.5,5.5, 6),
            'trk_energy':np.linspace(0.0,0.7,16),
            'shr_theta_v':np.linspace(0,3.2,21), 
            'shr_px_v':np.linspace(-1.0,1.0,21), 
            'trk_end_x_v':np.linspace(0,260,21),
            'pfnplanehits_V':np.linspace(0,300,21),
            'pfnplanehits_U':np.linspace(0,300,21),
            'pfnplanehits_Y':np.linspace(0,300,21),
            'shr_phi_v':np.linspace(-3.2,3.2,21),
            'shr_pz_v':np.linspace(-1.0,1.0,21),
            'trk_theta_v':np.linspace(0,3.2,21),
            'trk_phi_v':np.linspace(-3.2,3.2,21),
            'trk_energy_hits_tot':np.linspace(0,0.5,21),
            'trk_energy_tot':np.linspace(0,0.75,21),
            'trk_dir_z_v':np.linspace(-1.0,1.0,21),
            'SliceCaloEnergy2':np.linspace(0,150,21),
            'nu_flashmatch_score':np.linspace(0,15,21),
            'shr_energy_tot':np.linspace(0,0.25,21),
            'trk_chipr_best':np.linspace(0,500,21),
            'contained_sps_ratio':np.linspace(0,1.0,11),
            'contained_fraction':np.linspace(0,1.0,11),
            'trk_calo_energy_u_v':np.linspace(0,500,21),
            'trk_score_v':np.linspace(0,1.0,21),
            'nslice':np.linspace(-0.5,1.5,2),
            'flash_time':np.linspace(0,25,40),
            'crtveto':np.linspace(-0.5,1.5,2),
            'n_showers':np.linspace(-0.5,5.5,7),
            'NeutrinoEnergy2':np.linspace(0,500,21),
            'trk_score':np.linspace(0.5,1.0,11),
            'min_x':np.linspace(9,253,41),
            'max_x':np.linspace(9,253,41),
            'min_y':np.linspace(-112,112,41),
            'max_y':np.linspace(-112,112,41),
            'min_z':np.linspace(14,1020,41),
            'max_z':np.linspace(14,1020,41)}