close all;
clear all;

%% Initialise
HV = '25'; %'02'
config_file = ['Sim_DirConfigForwardComputation_HV',HV];
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

% When designing the parameter space, a range could be potentially chosen as
% a [0.05 10]
% b [0.5 30]
% af [0.05 40]
% bf [0.5 30]
% as [0.05 40]
% bs [0.05 30]
% afs [0.05 10]
% bfs [0.5 30]
Bounds = [0.05 10; % the lower bound might be too low? at least a=0.06 and b=21 crushed abaqus
        0.5 30;
        0.05 40;
        0.5 30;
        0.05 40;
        0.05 30;
        0.05 10;
        0.5 30];

Bounds_scale = Bounds(:,2) - Bounds(:,1);

%% Latin Hypercube Design for the ALL parameters
P = 8; % we have 8 parmetersb
N = P*10; % the rule for the LH design
Theta = lhsdesign(N,P);
Theta = bsxfun(@times, Theta, Bounds_scale');
Theta = bsxfun(@plus, Theta, Bounds(:,1)');


%% Prelocation for the results
ResLVVol = zeros(N,1);
ResLVStrain = zeros(N,strain_index_24+1-strain_index_1);

% SimRes.strainAbaTotal = {}; 
% SimRes.strainAbaTotal_fsn = {}; 
% SimRes.strainAbaTotalOnEachSlice = {};
% SimRes.strainAbaTotalOnEachSlice_fsn= {}; 
% SimRes.LVVolumeAba = {}; 
% SimRes.LVVolumeOri = {}; 
% SimRes.dis = {}; 
% SimRes.fiberStrain_cra = {}; 
% SimRes.fiberStrain_fsn = {}; 
% SimRes.node_ori = {}; 
% SimRes.elem = {}; 
% SimRes.nodeIndex_endo = {}; 
% 
% SimRes(N).LVVolumeOri = NaN;            
 
Times = zeros(N,1);
ii_start = 1;

if exist(['Init_LHS_ALL_HV',HV,'.mat'],'file')
    load(['Init_LHS_ALL_HV',HV,'.mat'],'ResLVVol','ResLVStrain','Theta','Times')
    ii_start = find(Times, 1, 'last' ) + 1;
end

for ii = ii_start : N
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
                    Sim_LVPassiveForwardSimulationMatPres(A,   B, ...
                                                      Af,  Bf, ...
                                                      As,  Bs, ...
                                                      Afs, Bfs, ...
                                                      pres, ...
                                                      postprocess_only); 
    %% save the volume and 24 strains
    ResLVVol(ii,1) = LVVolumeAba;
    ResLVStrain(ii,:) = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24);
%     SimRes(ii) = strainComparison;
    Times(ii,1) = toc;
    save(['Init_LHS_ALL_HV',HV,'.mat'],'ResLVVol','ResLVStrain','Theta','mpara','Times','HV')
end

save(['Init_LHS_ALL_HV',HV,'.mat'],'ResLVVol','ResLVStrain','Theta','mpara','Times','HV')
