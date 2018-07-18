% clear all;
close all;

HV = '02';% '25';

%% Load LV measured end-diastolic volume and strains
load(['DataMRI_HV',HV,'.mat'], 'HV','LVEDVMRI','strainData','strainDataFlag',...
    'strain_index_1','strain_index_24')
% load(['DataMRI_HV',HV,'.mat']);


%% Load the results and compute the objective functions

%% For HG step 1 grid selection
load(['Sim_FS_LV_Ca_Cb_HV',HV,'.mat'],'ResLVVol','ResLVStrain','SimRes','Ca','Cb','mpara','Times')

ResLVStrain_HG = ResLVStrain;
ResLVVol_HG = ResLVVol;
strainData_HG = strainData;
LVEDVMRI_HG = LVEDVMRI;

fO1_HG = Obj_fun(1,ResLVStrain_HG,ResLVVol_HG,strainData_HG,LVEDVMRI_HG);
fO2_HG = Obj_fun(2,ResLVStrain_HG,ResLVVol_HG,strainData_HG,LVEDVMRI_HG);

%% For BO init based on LHD
% load(['Init_LHS_Ca_Cb_HV',HV,'.mat'])
% N = size(SimRes,2);
% 
% fO1 = Obj_fun(1,ResLVStrain,ResLVVol,strainData,LVEDVMRI);
% fO2 = Obj_fun(2,ResLVStrain,ResLVVol,strainData,LVEDVMRI);
% 
% % [fO1_max, ind_max] = max(fO1);
% % Cab_max = Cab(ind_max,:);
% % 
% % Cab(ind_max,:) = [];
% % fO1(ind_max) = [];
% % fO2(ind_max) = [];

 
 
%% Plot the results
% Reshape vector to matrices for visulatisation
% vertical axis: Ca 
% horizontal axis: Cb
% reshape(Times,10,10)';

[X,Y] = meshgrid(Ca, Cb);
Z1 = reshape(fO1_HG,10,10)'; % Rows: varying Ca; Columns: varying Cb 
Z2 = reshape(fO2_HG,10,10)';

figure(60)
subplot(1,2,1)
hold on
surf(X,Y,Z1,'edgecolor','none')
% scatter3(Cab(:,2),Cab(:,1),fO1,'MarkerEdgeColor','red')
hold off
title(['HV',HV,': fO1 (unnormalized volume)'])
xlabel('Cb')
ylabel('Ca')
subplot(1,2,2)
hold on
surf(X,Y,Z2,'edgecolor','none')
% scatter3(Cab(:,2),Cab(:,1),fO2,'MarkerEdgeColor','red')
hold off
title(['HV',HV,': fO2 (normalized volume)'])
xlabel('Cb')
ylabel('Ca')
% legend('HG Step 1','BO Init LHS')
