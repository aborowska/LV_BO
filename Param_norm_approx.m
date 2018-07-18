% mpara.A = 0.23619;
% mpara.B = 10.810;
% mpara.Af = 20.037;
% mpara.Bf = 14.154;
% mpara.As = 3.7245;
% mpara.Bs = 5.1645;
% mpara.Afs = 0.41088;
% mpara.Bfs = 11.300;
parmean = [ 0.23619; 10.810; 20.037; 14.154; 3.7245; 5.1645; 0.41088; 11.300]';


% ~~ Three standard deviations account for 99.7% of the sample population being studied,
Bounds = [0.1 10; % the lower bound might be too low? at least a=0.06 and b=21 crushed abaqus
        0.5 30;
        0.1 40;
        0.5 30;
        0.1 40;
        0.1 30;
        0.05 10;
        0.1 30]';

parstd = (Bounds(2,:) - Bounds(1,:))./[6 6 6 6 6 6 6 6];
% parstd = (Bounds(2,:) - Bounds(1,:))./[10 6 6 6 10 6 10 6];

% Table 2 Correlation coefficient SCM with LV volume
% a b af bf as bs afs bfs
parcorr_fO1 = ...
   [1.0000 1.0000 0.9999 0.9998 0.9570 0.0000 0.9999 1.0000
    0.0000 1.0000 0.9999 0.9998 0.9571 0.0000 0.9999 1.0000
    0.0000 0.0000 1.0000 1.0000 0.9572 0.0000 0.9997 0.9999
    0.0000 0.0000 0.0000 1.0000 0.9572 0.0000 0.9997 0.9999
    0.0000 0.0000 0.0000 0.0000 1.0000 0.1000 0.9568 0.9570
    0.0000 0.0000 0.0000 0.0000 0.0000 1.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 1.0000 0.9999
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 1.0000];
 
parcorr_fO1 = triu(parcorr_fO1,1) +  triu(parcorr_fO1,1)' + eye(8);
% Table 3 Correlation coefficient SCM with normalized LV volume
% a b af bf as bs afs bfs
parcorr_fO2 = ...
   [1.0000  0.9918 -0.0505 -0.1202 -0.0100  0.0200  0.9155  0.8590
    0.0000  1.0000 -0.0940 -0.1600 -0.0020  0.0200  0.8865  0.7985
    0.0000  0.0000  1.0000  0.9952  0.0600 -0.0060 -0.2722  0.2759
    0.0000  0.0000  0.0000  1.0000  0.0400  0.0000 -0.3442  0.2026
    0.0000  0.0000  0.0000  0.0000  1.0000  0.0040 -0.0200 -0.0040
    0.0000  0.0000  0.0000  0.0000  0.0000  1.0000  0.0040  0.0000
    0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  1.0000  0.8296
    0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  1.0000];
  
parcorr_fO2 = triu(parcorr_fO2,1) +  triu(parcorr_fO2,1)' + eye(8);

% /** convert correlation matrix to covariance matrix **/
% R = {1.00 0.25 0.90,
%      0.25 1.00 0.50,
%      0.90 0.50 1.00 };
%  
% /** standard deviations of each variable **/
% c = {1  4  9};
% D = diag(c);
%  
% S = D*R*D; /** covariance matrix **/
D = diag(parstd);
parcov_fO1 = D*parcorr_fO1*D;
parcov_fO2 = D*parcorr_fO2*D;
