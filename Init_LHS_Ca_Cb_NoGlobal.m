%%  Ca and Cb optimization

%% Data and initialisation
% The data from the in vivo measurements consist of 24 regional
% circumferential strains and LV cavity volume at ED (end-diastolic).
% Because the ventricular pressure recording is not available, a
% population-based ED pressure (8 mmHg) is assumed.
Sim_DirConfigForwardComputation_HV02_NoGlobal;

postprocess_only = false;

% strain data in 24 segs, which will depend on the images we have
% if there are 5 images before apical region, then we will skip the first
strain_index_1 = 1;
strain_index_24 = 24;

% The pressure is chosen to be 8 mmHg for healthy volunteers.
pres = 8.0; %%mmHg, end-diastolic pressure loading


%% Initial parameter values (from Wang's paper)
% They are from ex-vivo experimental data
% provided from original Holzapfel and Ogden paper, 2009.

mpara.A = 0.23619;
mpara.B = 10.810;
mpara.Af = 20.037;
mpara.Bf = 14.154;
mpara.As = 3.7245;
mpara.Bs = 5.1645;
mpara.Afs = 0.41088;
mpara.Bfs = 11.300;


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

% load('try_single_call.mat','LVVolumeAba','strainAbaTotal','strainComparison');


%% 
parfor ii = 1 : N
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
                    Sim_LVPassiveForwardSimulationMatPres_NoGlobal(A,   B, ...
                                                      Af,  Bf, ...
                                                      As,  Bs, ...
                                                      Afs, Bfs, ...
                                                      pres, ...
                                                      postprocess_only, ...
                                                      options); 
    %%to save the volume and 24 strains
    ResLVVol(ii,1) = LVVolumeAba;
    ResLVStrain(ii,:) = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24);
    SimRes(ii) = strainComparison;
    Times(ii,1) = toc;
end

save('Init_LHS_Ca_Cb_NoGlobal.mat','ResLVVol','ResLVStrain','SimRes','Cab','mpara','Times')



% strainComparison = strainComparisonAbaqusVsMRI();
% 
% %feval_total = strainComparison.strainDiff;  %%this is not useful
% %feval_total(end+1) = (strainComparison.LVVolumeAba-strainComparison.LVVolumeMRI)/strainComparison.LVVolumeMRI;
% 
% %feval_total_fmincon = sum(feval_total.^2);


%% Some constrains when deciding whether the parameter set is reasonable or not.
% (1) First based on ejection fraction
% the possible ejection fraction range will be from 20% to 75%, so that it covers most of situation. 
% EF is calculated as
% EF = (LV_EDV –LV_ESV)/LV_EDV
% In the simulator, it can be approximated as
% EF = (strainComparison.LVVolumeAba - strainComparison.LVVolumeOri)/ strainComparison.LVVolumeAba
% (2) Secondly, there is also a constraint on the LV_EDV, which should be too large, but this depends the original size of the heart.
% From the same wiki page, stroke volume is around 95 mL (± 14 mL), 
% that means the maximum is around 120mL, 
% so we can limit the maximum LVEDV for that subject should be always less than 120 mL, or
% strainComparison.LVVolumeAba - strainComparison.LVVolumeOri < 120

xxx = [SimRes.node_ori];
xxx1 = reshape(xxx,[31856,3,20]);

LVVolumeAba = [SimRes.LVVolumeAba]';
LVVolumeOri = [SimRes.LVVolumeOri]';
EF =  (LVVolumeAba - LVVolumeOri)./LVVolumeAba;
LV_EDV = LVVolumeAba - LVVolumeOri;
ind_ok = ((EF >= 0.2) & (EF <= 0.75) & (LV_EDV < 120));
Cab(~ind_ok,:) % unreasonable parameters
Cab(ind_ok,:) % reasonable parameters

%% Step 2: optimize Ca Cb with fO1 = sum((eps_i-esp_i^*).^2) + (V-V0)^2
% eps_i^* is the measure regional circumferenctial strain
% V0 is the measure LV cavity volume from cine images
% fO1 is volume dominated
fO1 =  ([SimRes.LVVolumeAba]' -  [SimRes.LVVolumeOri]').^2 + ...
    sum(().^2);