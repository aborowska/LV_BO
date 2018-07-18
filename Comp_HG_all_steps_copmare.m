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

load(['Results/HG_all_steps_HV',HV,'.mat'],'y_HG','x_HG')
ResLVVol_HG = zeros(4,1);
ResLVStrain_HG = zeros(4,24);
Times_HG = zeros(4,1);

for ii = 1:4
    A = x_HG(ii,1);
    B = x_HG(ii,2);
    Af = x_HG(ii,3); 
    Bf = x_HG(ii,4); 
    As = x_HG(ii,5);
    Bs = x_HG(ii,6); 
    Afs = x_HG(ii,7);
    Bfs = x_HG(ii,8);

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
    ResLVVol_HG(ii,1)  = LVVolumeAba;
    ResLVStrain_HG(ii,:)  = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24);
    Times_HG(ii,1) = toc;
end

load(['Results/DataMRI_HV',HV,'.mat']);
y_check = Obj_fun(2,ResLVStrain_HG,ResLVVol_HG,strainData,LVEDVMRI);

save(['Results/HG_all_steps_HV',HV,'check.mat'],'ResLVVol_HG','ResLVStrain_HG',...
    'x_HG','y_HG','y_check','Times_HG','HV')