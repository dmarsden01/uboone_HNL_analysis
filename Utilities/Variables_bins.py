import numpy as np

bins_var = {'shrclusdir2':np.linspace(0,3.2, 21),
            'n_tracks':np.linspace(0.5,5.5, 6),
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
            'trk_dir_z_v':np.linspace(-1.0,1.5,21),
            'SliceCaloEnergy2':np.linspace(0,150,21)}