close all; clear all; clc;
%% Initially prepared by Hao Gao
%% here we will optimize af bf
 
clear all; close all; clc;
 
%% load gloabal options for HV1. Useful data will be saved in globla option structure
HV = '25';
config_file = ['Sim_DirConfigForwardComputation_HV',HV];
% Sim_DirConfigForwardComputation_HV02;
eval(config_file);


%% read the strain data from images and LV measured end-diastolic volume
load(['DataMRI_HV',HV,'.mat'], 'HV','LVEDVMRI','strainData','strainDataFlag',...
    'strain_index_1','strain_index_24');

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


%% before calling the optimization function, we will need to update the paramters from the second step	
load(['Sim_FS_LV_Ca_Cb_Klotz_HV',HV,'.mat'],'mpara')
%       A: 0.0218
%       B: 4.1625
%      Af: 1.8485
%      Bf: 5.4501
%      As: 0.3436
%      Bs: 1.9886
%     Afs: 0.0379
%     Bfs: 4.3512
% % mpara.A = 0.100000;
% % mpara.B = 7.545964;
% % mpara.Af = 3.886494;
% % mpara.Bf = 9.880257;
% % mpara.As = 0.722426;
% % mpara.Bs = 3.605100;
% % mpara.Afs = 0.100000;
% % mpara.Bfs = 7.888011;
options.mpara = mpara; %%remember to update here


options.strainData= strainData;
options.strainDataFlag = strainDataFlag;
options.LVEDVMRI = LVEDVMRI; %LVVolumeMeasured.edv;
options.LVEDP_high = 30;
options.LVEDP_norm = 8;
options.strain_index_1 = 1;
options.strain_index_24 = 24;


cd(abaqusDir);
fid_log = fopen(opt_log_filename,'a');
cd(workingDir);
fprintf(fid_log, '\n \n beginning of step 3 for af bf\n \n');

options_Opt = optimset('Algorithm', 'sqp', 'TolFun', 1e-6, ...
                           'TolX',0.000001,'Diffminchange',1e-4,...
                           'Diffmaxchange',1e-3, 'MaxIter', 100);
afmin = 0.1;
afmax = 5;
bfmin = 0.1;    
bfmax = 5;

Lb = [afmin bfmin ];
Ub = [afmax bfmax ];
x0 = [1.0  1.0];

tic
[x,fval,exitflag,output] = fmincon(@LVPassive_model_optimization_fmincon_afbf,...
    x0,[],[],[],[],Lb,Ub,[],options_Opt);
Times = toc;


save(['Sim_FS_LV_af_bf_HV',HV,'.mat'],...
    'x','fval','exitflag','output','mpara','mpara_old','Times')

mpara_old = mpara;
mpara.Af = mpara.Af*x(1);
mpara.Bf = mpara.Bf*x(2);


save(['Sim_FS_LV_af_bf_HV',HV,'.mat'],...
    'x','fval','exitflag','output','mpara','mpara_old','Times')
