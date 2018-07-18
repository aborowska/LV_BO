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

%% Parameter ranges
parnames = {'A','B','Af','Bf','As','Bs','Afs','Bfs'};

Param_norm_approx;
N = 200; 
Theta = [];
while size(Theta,1) < 50
    % Theta_fO1 = mvnrnd(parmean,parcov_fO1,N);
    Theta_fO2 = mvnrnd(parmean,parcov_fO2,N);

    ind_ok = bsxfun(@gt, Theta_fO2, Bounds(1,:))  &    bsxfun(@lt, Theta_fO2, Bounds(2,:)) ;
    ind = all(ind_ok,2);
    sum(ind_ok,1)
    Theta = Theta_fO2(ind,:);
end
N = size(Theta,1);

%% Prelocation for the results
ResLVVol = zeros(N,1);
ResLVStrain = zeros(N,strain_index_24+1-strain_index_1);
      
 
Times = zeros(N,1);
ii_start = 1;

if exist(['Init_MVN_ALL_HV',HV,'.mat'],'file')
    load(['Init_MVN_ALL_HV',HV,'.mat'],'ResLVVol','ResLVStrain','Theta','Times')
    ii_start = find(Times, 1, 'last' ) + 1;
end

for ii = ii_start : N
    fprintf('Init_MVN_ALL_v2 >>>>> iter = %ii\n',ii)
    tic
    A = Theta(ii,1);
    B = Theta(ii,2);
    Af = Theta(ii,3); 
    Bf = Theta(ii,4); 
    As = Theta(ii,5);
    Bs = Theta(ii,6); 
    Afs = Theta(ii,7);
    Bfs = Theta(ii,8);


    [LVVolumeAba, strainAbaTotal, strainComparison, SuccessB] = ...
                    Sim_LVPassiveForwardSimulationMatPres_NoGlobal(A,   B, ...
                                                      Af,  Bf, ...
                                                      As,  Bs, ...
                                                      Afs, Bfs, ...
                                                      pres, ...
                                                      postprocess_only,...
                                                      options); 
    %% save the volume and 24 strains
    ResLVVol(ii,1) = LVVolumeAba;
    ResLVStrain(ii,:) = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24);
%     SimRes(ii) = strainComparison;
    Times(ii,1) = toc;
    save(['Init_MVN_ALL_HV',HV,'.mat'],'ResLVVol','ResLVStrain','Theta','mpara','Times','HV')
end

save(['Init_MVN_ALL_HV',HV,'.mat'],'ResLVVol','ResLVStrain','Theta','mpara','Times','HV')
