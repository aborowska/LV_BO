% clear all; close all
HV = '25';
load(['Results/DataMRI_HV',HV,'.mat']);

x_HG = zeros(4,8);
y_HG = zeros(4,1);

load(['Results/Sim_FS_LV_Ca_Cb_HV',HV,'.mat'],'ResLVVol','ResLVStrain','SimRes','Ca','Cb','mpara','Times')
fO1_HG = Obj_fun(1,ResLVStrain,ResLVVol,strainData,LVEDVMRI);
fO2_HG = Obj_fun(2,ResLVStrain,ResLVVol,strainData,LVEDVMRI);
[fO1_min, ind_min] = min(fO1_HG);
ind_b = mod(ind_min,length(Cb));
ind_a = (ind_min-ind_b)/length(Cb);
Ca = Ca(ind_a+1); 
Cb = Cb(ind_b);
x_HG(1,:) = [0.23619*Ca, 10.810*Cb, 20.037*Ca,14.154*Cb, 3.7245*Ca, 5.1645*Cb, 0.41088*Ca, 11.300*Cb];
y_HG(1) = fO2_HG(ind_min);


% load(['Results/Sim_FS_LV_Ca_Cb_Klotz_HV',HV,'.mat'])
load(['Results/HG_Klotz_HV',HV,'.mat'],'mpara','ResLVVol','ResLVStrain')
x_HG(2,:) = cell2mat(struct2cell(mpara))';
y_HG(2,:) = Obj_fun(2,ResLVStrain',ResLVVol,strainData,LVEDVMRI);

load(['Results/Sim_FS_LV_af_bf_HV',HV,'.mat'],'mpara','fval')
x_HG(3,:) = cell2mat(struct2cell(mpara))';
y_HG(3,:) = fval;


load(['Results/Sim_FS_LV_a_b_HV',HV,'.mat'],'mpara','fval');
x_HG(4,:) = cell2mat(struct2cell(mpara))';
y_HG(4,:) = fval;

save(['Results/HG_all_steps_HV',HV,'.mat'],'y_HG','x_HG')
% Comp_HG_all_steps;
load(['Results/HG_all_steps_HV',HV,'.mat'],'y_HG','x_HG','ResLVVol_HG','ResLVStrain_HG');

%%
load(['Results/BayOpt_ALL_afterCaCb_HV',HV,'.mat'],'i1','ResLVVol','ResLVStrain','x','y','Cab')
T0 = 10; %46; % 39; %10;%20; %79;
T1 = 29; %146; %126; %39; %98;
T = 29; % 146; %126; %29; %39; %98; % 117;

figure(1)
subplot(3,1,1)
plot(1:T0,ResLVVol(1:T0))
hold on
plot((T0+1):T1,ResLVVol((T0+1):T1),'r')
plot((T1+1):T,ResLVVol((T1+1):T),'y')
plot(1:T,LVEDVMRI+0*(1:T),'m')
hold off
title('Volume')
xlim([1,T])

subplot(3,1,2)
plot(1:T0,mean(ResLVStrain(1:T0,:),2))
hold on
plot((T0+1):T1,mean(ResLVStrain((T0+1):T1,:),2),'r')
plot((T1+1):T,mean(ResLVStrain((T1+1):T,:),2),'y')
plot(1:T,mean(strainData(1:24))+0*(1:T),'m')
hold off
title('Mean strain')
xlim([1,T])
% legend('Init LHS','1st BO run','2nd BO run','data')
% legend('Init LHS','1st BO run','data')
% legend('Init 50% min after Ca Cb','1st BO run','data')
legend('Init w/ BO 4 params','BO 8 params','data')

[y_min, ind_min] = min(y);

subplot(3,1,3)
plot(1:T0,y(1:T0))
hold on
plot((T0+1):T1,y((T0+1):T1),'r')
plot((T1+1):T,y((T1+1):T),'y')
plot(1:T,y_min+0*(1:T),'g')
YL = ylim;
line([ind_min ind_min],YL,'Color','green')
hold off
title('f02 value (normalised volume)')
xlim([1,T])
legend('Init w/ BO 4 params','BO 8 params','min fO2')



%%
[y_min, ind_min] = min(y);

figure(2)
subplot(3,1,1)
plot(1:T0,ResLVVol(1:T0))
hold on
plot((T0+1):T1,ResLVVol((T0+1):T1),'r')
% plot((T1+1):T,ResLVVol((T1+1):T),'y')
plot(1:T,LVEDVMRI+0*(1:T),'m')

plot(1:T,ResLVVol(4,1)+0*(1:T),'k')
plot(1:T,ResLVVol(1,1)+0*(1:T),'k:')

hold off
title('Volume')
xlim([1,T])
% legend('Init w/ BO 4 params','BO 8 params',...
%     'min fO2','HGO final','HGO step1')
legend('Init w/ BO Ca Cb (50% best)','BO 8 params',...
    'min fO2','HGO final','HGO step1')



subplot(3,1,2)
plot(1:T0,mean(ResLVStrain(1:T0,:),2))
hold on
plot((T0+1):T1,mean(ResLVStrain((T0+1):T1,:),2),'r')
plot((T1+1):T,mean(ResLVStrain((T1+1):T,:),2),'y')
plot(1:T,mean(strainData(1:24))+0*(1:T),'m')

plot(1:T,mean(ResLVStrain_HG(4,:))+0*(1:T),'k')
plot(1:T,mean(ResLVStrain_HG(1,:))+0*(1:T),'k:')

hold off
title('Mean strain')
xlim([1,T])
% legend('Init w/ BO 4 params','BO 8 params',...
%     'min fO2','HGO final','HGO step1')
legend('Init w/ BO Ca Cb (50% best)','BO 8 params',...
    'data','HGO final','HGO step1')


subplot(3,1,3)
plot(1:T0,y(1:T0))
hold on
plot((T0+1):T1,y((T0+1):T1),'r')
% plot((T1+1):T,y((T1+1):T),'y')
plot(1:T,y_min+0*(1:T),'g')

plot(1:T,y_HG(4,1)+0*(1:T),'k')
plot(1:T,y_HG(1,1)+0*(1:T),'k:')

YL = ylim;
line([ind_min ind_min],YL,'Color','green')
hold off
title('f02 value (normalised volume)')
xlim([1,T])
% legend('Init w/ BO 4 params','BO 8 params',...
%     'min fO2','HGO final','HGO step1')
legend('Init w/ BO Ca Cb (50% best)','BO 8 params',...
    'min fO2','HGO final','HGO step1')



 

%%
parnames = {'A','B','Af','Bf','As','Bs','Afs','Bfs'};

figure(3)
for ii = 1:8
    subplot(2,4,ii)

    plot(1:T0,x(1:T0,ii))
    hold on
    plot((T0+1):T1,x((T0+1):T1,ii),'r')
    % plot((T1+1):T,y((T1+1):T),'y')
    plot(1:T,x(ind_min,ii) +0*(1:T),'g')

    plot(1:T,x_HG(4,ii)+0*(1:T),'k')
    plot(1:T,x_HG(1,ii)+0*(1:T),'k:')
    title(parnames(ii))
end

legend('Init w/ BO Ca Cb (50% best)','BO 8 params',...
    'min fO2','HGO final','HGO step1')
