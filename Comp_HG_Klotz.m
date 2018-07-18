close all;
clear all;

%% Initialise
HV = '25'; %'02'
config_file = ['Sim_DirConfigForwardComputation_HV',HV,'_NoGlobal'];
eval(config_file);
%Sim_DirConfigForwardComputation_HV25;
%Sim_DirConfigForwardComputation_HV02;

postprocess_only = false;

% strain data in 24 segs, which will depend on the images we have
% if there are 5 images before apical region, then we will skip the first
strain_index_1 = 1;
strain_index_24 = 24;

pres = 8.0; %%mmHg, end-diastolic pressure loading

load(['Results/Sim_FS_LV_Ca_Cb_Klotz_HV',HV,'.mat'],'mpara')


A = mpara.A;
B = mpara.B;
Af = mpara.Af; 
Bf = mpara.Bf; 
As = mpara.As;
Bs = mpara.Bs; 
Afs = mpara.Afs;
Bfs = mpara.Bfs;

tic
[LVVolumeAba, strainAbaTotal, strainComparison, SuccessB] = ...
                Sim_LVPassiveForwardSimulationMatPres_NoGlobal(A,   B, ...
                                                  Af,  Bf, ...
                                                  As,  Bs, ...
                                                  Afs, Bfs, ...
                                                  pres, ...
                                                  false,...
                                                  options); 
                                              
%% save the volume and 24 strains
ResLVVol  = LVVolumeAba;
ResLVStrain  = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24);
Times = toc;
save(['Results/HG_Klotz_HV',HV,'.mat'],'ResLVVol','ResLVStrain','mpara','Times','HV')