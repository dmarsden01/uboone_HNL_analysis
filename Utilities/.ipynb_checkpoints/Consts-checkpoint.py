

cosmic="(trk1_backtracked_pdg==0 or trk2_backtracked_pdg==0)"
notcosmic="~"+cosmic
diff="(trk1_backtracked_pdg!=trk2_backtracked_pdg)"

q_wellreco=diff+"and"+notcosmic


flash_time_cutmax=16.44
flash_time_cutmin=5.64
flash_time_cut=f"flash_time>{flash_time_cutmin} and flash_time<{flash_time_cutmax}"

nu_flashmatch_score_cutmax=20
nu_flashmatch_score_cut=f"nu_flashmatch_score<{nu_flashmatch_score_cutmax}"


##TPC
TPCxlo=0
TPCxhi=256
TPCylo=-115
TPCyhi=117
TPCzlo=0
TPCzhi=1037


# ##Krish

# TPCxlo=8.45
# TPCxhi=244.8
# TPCylo=-106.5
# TPCyhi=106.5
# TPCzlo=5
# TPCzhi=1031.8


verinfid_cut=f"VertLoc_x>{TPCxlo} and VertLoc_x<{TPCxhi} and VertLoc_y>{TPCylo} and VertLoc_y<{TPCyhi} and VertLoc_z>{TPCzlo} and VertLoc_z<{TPCzhi}"  #TPC boundry from pawels internal note


max_x_cut=253
min_x_cut=9#9
max_y_cut=112
min_y_cut=-112
max_z_cut= 1020#1020
min_z_cut=14

# max_x_cut=244.8
# min_x_cut=8.45#9
# max_y_cut=106.5
# min_y_cut=-106.5
# max_z_cut= 1031.8#1020
# min_z_cut=5

#f"min_x>{min_x_cut} and max_x<{max_x_cut} 

max_extent_cut=f"min_y>{min_y_cut} and max_y<{max_y_cut} and min_z>{min_z_cut} and max_z<{max_z_cut}" + " and "+f"min_x>{min_x_cut} and max_x<{max_x_cut}"

closestcos_cutmax=20
closestcos_cut=f"_closestNuCosmicDist>{closestcos_cutmax}"

npfp_cutmax=4
npfp_cut=f"n_pfps<={npfp_cutmax}"


ntrack_cutmax=3
ntrack_cut=f"n_tracks<={ntrack_cutmax}"

nshower_cutmax=2
nshower_cut=f"n_showers<={nshower_cutmax}"

q0="n_pfps>-99999"

cosangle_cutmin=-0.94
cosangle_cut=f"cosangle>{cosangle_cutmin}"

proton_cutmin=-0.5
proton_cut=f"trk1_llr_pid_score_v>{proton_cutmin} and  trk2_llr_pid_score_v>{proton_cutmin}"

nueng_cutmax=500
nueng_cut=f"NeutrinoEnergy2<{nueng_cutmax}"

trk1len_cutmax=50
trk1len_cut=f"trk1_len_v<{trk1len_cutmax}"


# q1=q0+ " & "  + verinfid_cut
# q2=q1+" & "+ flash_time_cut
# q3=q2+" & "+ nu_flashmatch_score_cut
# q4=q3+ " & "  +npfp_cut+ " & " + ntrack_cut + " & " + nshower_cut
# q5=q4+ " & " + max_extent_cut
# q6=q5+ " & " + cosangle_cut
# q7=q6+ " & " + proton_cut
# q8=q7+ " & " + nueng_cut
# q9=q8+ " & " + trk1len_cut
# q_preselection_run1=q9


# q_crtveto=q_preselection_run1 +  ' & crtveto==0'
# q_closecut=q_crtveto + " & " + closestcos_cut
# q_preselection_run3=q_closecut

# q_preselection_list_run1=[q0,q1,q2,q3,q4,q5,q6,q7,q8,q9]
# q_preselection_list_run3=[q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q_crtveto,q_closecut]

# q_preselection_names_run1=["Has Vertex","Fid vol","Flash Time cut",f"Flash Match score<{nu_flashmatch_score_cutmax}","PFP mult cut","Max extent (x,y,z)",f"Opening Angle>-{cosangle_cutmin}","Proton Cut","Energy Reco Cut","Longest Track Lenght"]
# q_preselection_names_run3=["Has Vertex","Fid vol","Flash Time cut",f"Flash Match score<{nu_flashmatch_score_cutmax}","PFP mult cut","Max extent (x,y,z)",f"Opening Angle>-{cosangle_cutmin}","Proton Cut","Energy Reco Cut","Longest Track Lenght","CRT veto","Closest Cosmic Dist (CRT)"]




q1=q0+" & "+ flash_time_cut
q2=q1+" & "+ nu_flashmatch_score_cut
q3=q2+ " & "  +npfp_cut+ " & " + ntrack_cut + " & " + nshower_cut
q4=q3+ " & " + max_extent_cut
q5=q4+ " & " + cosangle_cut
q6=q5+ " & " + proton_cut
q7=q6+ " & " + nueng_cut
q8=q7+ " & " + trk1len_cut
q_preselection_run1=q8


q_crtveto=q_preselection_run1 +  ' & crtveto==0'
q_closecut=q_crtveto + " & " + closestcos_cut
q_preselection_run3=q_closecut

q_preselection_list_run1=[q0,q1,q2,q3,q4,q5,q6,q7,q8]
q_preselection_list_run3=q_preselection_list_run1+[q_crtveto,q_closecut]

q_preselection_names_run1=["Has Vertex","Flash Time",f"Flash Match Score","Object Multiplicity","Containment (x,y,z)",f"Opening Angle","Proton Likelihood","Maximum Energy","Longest Track Length"]
q_preselection_names_run3=q_preselection_names_run1+["CRT Veto","Closest CRT Cosmic"]

q_sigtraining= "nu_purity_from_pfp>0.90"






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



