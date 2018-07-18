%%  Ca and Cb optimization

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


%% Initial parameter values (from Wang's paper)
% They are from ex-vivo experimental data
% provided from original Holzapfel and Ogden paper, 2009.
% The pressure is chosen to be 8 mmHg for healthy volunteers.

mpara.A = 0.23619;
mpara.B = 10.810;
mpara.Af = 20.037;
mpara.Bf = 14.154;
mpara.As = 3.7245;
mpara.Bs = 5.1645;
mpara.Afs = 0.41088;
mpara.Bfs = 11.300;

pres = 8.0; %%mmHg, end-diastolic pressure loading

% inital params for step af, bf (from the  second step)
% before calling the optimization function, we will need to update the paramters from the second step	
% A = 0.100000;
% B = 7.545964;
% Af = 3.886494;
% Bf = 9.880257;
% As = 0.722426;
% Bs = 3.605100;
% Afs = 0.100000;
% Bfs = 7.888011; 


% materialParam_startLine = 19;
% pressure_loadingLine = 14;
% eplsion = 0.01;



%% Latin Hypercube Design for the range of Ca and Cb
% Ca = 0.1 : 0.1 : 1.0;
% Cb = 0.1 : 0.1 : 1.0;
P = 2; % we have 2 parmeters, Ca and Cb
N = P*10; % the rule for the LH design
Cab = lhsdesign(N,P);
Cab = Cab*0.9 + 0.1;


%% Prelocation for the results
ResLVVol = zeros(N,1);
ResLVStrain = zeros(N,strain_index_24+1-strain_index_1);

SimRes.strainAbaTotal = {}; 
SimRes.strainAbaTotal_fsn = {}; 
SimRes.strainAbaTotalOnEachSlice = {};
SimRes.strainAbaTotalOnEachSlice_fsn= {}; 
SimRes.LVVolumeAba = {}; 
SimRes.LVVolumeOri = {}; 
SimRes.dis = {}; 
SimRes.fiberStrain_cra = {}; 
SimRes.fiberStrain_fsn = {}; 
SimRes.node_ori = {}; 
SimRes.elem = {}; 
SimRes.nodeIndex_endo = {}; 

SimRes(N).LVVolumeOri = NaN;            
 
Times = zeros(N,1);

for ii = 1 : N
    tic
    A = Cab(ii,1)*mpara.A;
    B = Cab(ii,2)*mpara.B;
    Af = Cab(ii,1)*mpara.Af; 
    Bf = Cab(ii,2)*mpara.Bf; 
    As = Cab(ii,1)*mpara.As;
    Bs = Cab(ii,2)*mpara.Bs; 
    Afs = Cab(ii,1)*mpara.Afs;
    Bfs = Cab(ii,2)*mpara.Bfs;


    [LVVolumeAba, strainAbaTotal,strainComparison, SuccessB] = ...
                    Sim_LVPassiveForwardSimulationMatPres(A,   B, ...
                                                      Af,  Bf, ...
                                                      As,  Bs, ...
                                                      Afs, Bfs, ...
                                                      pres, ...
                                                      postprocess_only); 
    %%to save the volume and 24 strains
    ResLVVol(ii,1) = LVVolumeAba;
    ResLVStrain(ii,:) = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24);
    SimRes(ii) = strainComparison;
    Times(ii,1) = toc;
end

save(['Init_LHS_Ca_Cb_HV',HV,'.mat'],'ResLVVol','ResLVStrain','SimRes','Cab','mpara','Times','HV')
