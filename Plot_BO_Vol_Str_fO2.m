% clear all; close all
HV = '25';
load(['Results/DataMRI_HV',HV,'.mat']);
figure(1)
T0 = 46; % 39; %10;%20; %79;
T1 = 146; %126; %39; %98;
T = 146; %126; %29; %39; %98; % 117;

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

