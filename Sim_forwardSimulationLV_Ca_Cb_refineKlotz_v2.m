%% Initially prepared by Hao Gao
%% here we will optimize the volume according to the klotz curve

clear all; close all; clc;
 
%% load gloabal options for HV1. Useful data will be saved in globla option structure
HV = '25';
config_file = ['Sim_DirConfigForwardComputation_HV',HV];
% Sim_DirConfigForwardComputation_HV02;
eval(config_file);


%% read the strain data from images and LV measured end-diastolic volume
load(['DataMRI_HV',HV,'.mat'], 'HV','LVEDVMRI','strainData','strainDataFlag',...
    'strain_index_1','strain_index_24');

%% using the right parameter sets
%% A, B, Af, Bf, As, Bs, Afs, Bfs are the 8 unknown parameters for HGO
%% myocardial material model
% A = 0.224487;
% B = 1.621500;
% Af = 2.426717;
% Bf = 1.826862;
% As = 0.556238;
% Bs = 0.774678;
% Afs = 0.390516;
% Bfs = 1.695000;


pres = 8.0; %%mmHg, end-diastolic pressure loading


%% here we need to provide the eigen values, which is calculated
%% (LV-LV_mean)*eigVecs, LV-LV_mean needs to be a row vector
%% wc_pca = [1.0 1.0 1.0 1.0 1.0]; %%only use five eigen vectors for reconstruction the mesh
%  wc_pca = [-18.2695   -1.4049    2.8309  -14.6374  -25.4179]; %% that is from HV03

%% using how many cpus for parameter computing, 4 - 6 will be good choice
%% we might have only 48 tokens currently, 16 more are being purchased
%% so please use only maximum 40 tokens, which means if 4 token per matlab session,
%% then you can run up to 10 matlab sessions in parallel. 
%% for the current setting
%% options.cpunumber = 6;  %%default value
  

%% that part is for debugging, should not be changed once the model is setup properly
postprocess_only = false;
regenerate_fibre = false;


%% before calling the optimization function, we will need to update the global options
load(['Sim_FS_LV_Ca_Cb_HV',HV,'.mat'],'ResLVVol','ResLVStrain','SimRes','Ca','Cb','mpara','Times')

% use the Ca-Cb pair leading to the smallest value of the objective
% function
fO1_HG = Obj_fun(1,ResLVStrain,ResLVVol,strainData,LVEDVMRI);
fO2_HG = Obj_fun(2,ResLVStrain,ResLVVol,strainData,LVEDVMRI);
[fO1_min, ind_min] = min(fO1_HG);
ind_b = mod(ind_min,length(Cb));
ind_a = (ind_min-ind_b)/length(Cb);
Ca = Ca(ind_a+1); 
Cb = Cb(ind_b);


% Ca = 0.4;
% Cb = 0.3;
mpara.A = 0.23619*Ca;
mpara.B = 10.810*Cb;
mpara.Af = 20.037*Ca;
mpara.Bf = 14.154*Cb;
mpara.As = 3.7245*Ca;
mpara.Bs = 5.1645*Cb;
mpara.Afs = 0.41088*Ca;
mpara.Bfs = 11.300*Cb;
options.mpara = mpara; %%remember to update here


options.strainData= strainData;
options.LVEDVMRI = LVEDVMRI;
options.LVEDP_high = 30;
options.LVEDP_norm = 8;


Lb_aafs = [0.1 0.05];
Ub_aafs = [10.0 10.0];
options_Opt = optimset('Algorithm', 'sqp', 'TolFun', 1e-3, ...
                       'TolX',0.001,'Diffminchange',1e-3);
ca_cb_0 = [1.0 1.0];
tic
[ca_cb,fval,exitflag,output] = fmincon(@LVPassive_model_optimization_fmincon_ca_cb_klotz,...
    ca_cb_0,[],[], [], [], Lb_aafs,Ub_aafs,[],options_Opt);
Times = toc;

save(['Sim_FS_LV_Ca_Cb_Klotz_HV',HV,'.mat'],...
    'ca_cb','fval','exitflag','output','Ca','Cb','mpara','Times')
