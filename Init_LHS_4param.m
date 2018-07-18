%%  Ca and Cb optimization

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

%% Latin Hypercube Design for the range of theta_1, theta_2, theta_3, theta_4
P = 4; % we have 4 parmeters 
N = P*10; % the rule for the LH design
Theta = lhsdesign(N,P);
Theta = Theta*0.9 + 0.1; %%     WHAT RANGES FOR THE SCALINGS????? CURRENTLY: [0.1-1]

if false
    load(['BayOpt_4param_HV',HV,'_oldBounds.mat'],'x','y','i1')
%     load(['BayOpt_4param_HV',HV,'.mat'],'x','y','i1')

    y = y((end-i1+1+1):end);
    x = x((end-i1+1+1):end,:); 
    [y_best, ind_best] = sort(y);
    x = x(ind_best,:); 
    N = round(length(y)/2);
    y = y_best(1:N,:);
    x = x(1:N,:);
        
    for ii = 1:4
        subplot(2,2,ii)
        histogram(x(:,ii), 20)
        title(['Theta_',num2str(ii)]) 
        % CONCLUSIONS: 
        % * THETA_1 AND THETA_3 WITHING 0.1-1 BOUNDS 
        % * THETA_2 HITS THE LOWER BOUND 0.1
        % * THETA_4 HITS THE UPPER BOUND 1
    end
end

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

for ii = 1 : N
    tic
    A = Theta(ii,1)*mpara.A;
    B = Theta(ii,1)*mpara.B;
    Af = Theta(ii,2)*mpara.Af; 
    Bf = Theta(ii,3)*mpara.Bf; 
    As = Theta(ii,2)*mpara.As;
    Bs = Theta(ii,3)*mpara.Bs; 
    Afs = Theta(ii,4)*mpara.Afs;
    Bfs = Theta(ii,4)*mpara.Bfs;

    try
        [LVVolumeAba, strainAbaTotal,strainComparison, SuccessB] = ...
                        Sim_LVPassiveForwardSimulationMatPres_NoGlobal(A,   B, ...
                                                          Af,  Bf, ...
                                                          As,  Bs, ...
                                                          Afs, Bfs, ...
                                                          pres, ...
                                                          postprocess_only,...
                                                          options); 
        %%to save the volume and 24 strains
        ResLVVol(ii,1) = LVVolumeAba;
        ResLVStrain(ii,:) = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24);
    %     SimRes(ii) = strainComparison;
    catch
        fprintf('Abaqus crashed at iteration %i \n',ii)
        fprintf('For theta = %6.4f, %6.4f, %6.4f, %6.4f\n',Theta(ii,1),Theta(ii,2),Theta(ii,3),Theta(ii,4))
    end
    Times(ii,1) = toc;
end
save(['Init_LHS_4param_HV',HV,'.mat'],'ResLVVol','ResLVStrain','Theta','mpara','Times','HV') % ,'SimRes');
