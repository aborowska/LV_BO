%% Initially prepared by Hao Gao
%% here we will optimize the volume according to the klotz curve
% clear all; close all; %clc;
clear all; 
close all;

addpath(genpath('../../../GPstuff'));

%% load gloabal options for HV1. Useful data will be saved in globla option structure
HV = '25';
config_file = ['Sim_DirConfigForwardComputation_HV',HV,'_NoGlobal'];
% Sim_DirConfigForwardComputation_HV02;
eval(config_file);


%% read the strain data from images and LV measured end-diastolic volume
load(['DataMRI_HV',HV,'.mat'], 'HV','LVEDVMRI','strainData','strainDataFlag',...
    'strain_index_1','strain_index_24');

%% using the right parameter sets
%% A, B, Af, Bf, As, Bs, Afs, Bfs are the 8 unknown parameters for HGO
%% myocardial material model
mpara.A = 0.23619;
mpara.B = 10.810;
mpara.Af = 20.037;
mpara.Bf = 14.154;
mpara.As = 3.7245;
mpara.Bs = 5.1645;
mpara.Afs = 0.41088;
mpara.Bfs = 11.300;

vpara = [0.23619, 10.810, 20.037,  14.154, 3.7245, 5.1645, 0.41088, 11.300];

pres = 8.0; %%mmHg, end-diastolic pressure loading
 
%% Load the initial values based on BO for Ca and Cb
fprintf('*** Using BO for Ca and Cb *** \n')
load(['BayOpt_Ca_Cb_HV',HV,'.mat'],...
    'ResLVVol','ResLVStrain','x','y','i1')

% the loaded x is the scalining vector [Ca,Cb]
ResLVVol = ResLVVol((end-i1+1+1):end,:);
ResLVStrain = ResLVStrain((end-i1+1+1):end,:);
y = y((end-i1+1+1):end);
x = x((end-i1+1+1):end,:); 

[y_best, ind_best] = sort(y);
ResLVVol = ResLVVol(ind_best,:);
ResLVStrain = ResLVStrain(ind_best,:);
x = x(ind_best,:); 
% x_best = x(1,:); % optimal scaling

%% retrieve the 8 parameters to set the updated bounds for af and bf
N = round(length(y)/2);
y = y_best(1:N,:);
x = x(1:N,:);
ResLVVol = ResLVVol(1:N,:);
ResLVStrain = ResLVStrain(1:N,:);

Cab = x;
x = repmat(vpara,N,1);
x = x .* repmat(Cab,1,4); % now x is the set of 50% (=10 here) best 8-dim parameters
Bounds = [    0.1000    0.5000    0.1000    0.5000    0.1000    0.0500    0.1000    0.5000
             10.0000   30.0000   40.0000   30.0000   40.0000   30.0000   10.0000   30.0000];
Bounds(1,:) = min(x)*0.9;
Bounds(2,:) = max(x)*1.1;
Bounds_scale = Bounds(2,:) - Bounds(1,:);


%% choose the best scaling to fix the remaining parameters
% rescale the parameters
% mpara.A = Ca_Cb_opt(1)*mpara.A;
% mpara.B = Ca_Cb_opt(2)*mpara.B;
% mpara.Af = Ca_Cb_opt(1)*mpara.Af; 
% mpara.Bf = Ca_Cb_opt(2)*mpara.Bf; 
% mpara.As = Ca_Cb_opt(1)*mpara.As;
% mpara.Bs = Ca_Cb_opt(2)*mpara.Bs; 
% mpara.Afs = Ca_Cb_opt(1)*mpara.Afs;
% mpara.Bfs = Ca_Cb_opt(2)*mpara.Bfs;
mpara.A = x(1,1);
mpara.B = x(1,2);
mpara.Af = x(1,3);
mpara.Bf = x(1,4); 
mpara.As = x(1,5);
mpara.Bs = x(1,6);
mpara.Afs = x(1,7);
mpara.Bfs = x(1,8);    


%% set the initial point [af, bf] and the corresponding objective function for the current optimisation
x = x(1,3:4);
y = y(1);
ResLVVol = ResLVVol(1,:);
ResLVStrain = ResLVStrain(1,:);
P = size(x,2);

%% construct GP to model the (ABAQUS) function
len_sc = Bounds_scale(3:4)/10;
cfse = gpcf_sexp('lengthScale',len_sc,'lengthScale_prior',prior_t('s2',4),'magnSigma2',.1,'magnSigma2_prior',prior_sqrtt('s2',10^2));
% cfse = gpcf_sexp('lengthScale',1,'magnSigma2',1,'magnSigma2_prior',prior_sqrtt('s2',10^2)); %GPCF_SEXP  Create a squared exponential (exponentiated quadratic) covariance function
lik = lik_gaussian('sigma2', 0.001, 'sigma2_prior', prior_fixed);
% gp = gp_set('cf', {cfc, cfl, cfl2, cfse}, 'lik', lik);
gp = gp_set('cf', {cfse}, 'lik', lik);
 

%% ----- conduct Bayesian optimization -----

% Set the options for optimizer of the acquisition function
optimf = @fmincon;
optdefault = struct('GradObj','on','LargeScale','off','Algorithm','trust-region-reflective','TolFun',1e-9,'TolX',1e-6);
opt = optimset(optdefault);
% lb = [0.1 0.1];     % lower bound of the input space
% ub = [5 5];   % upper bound of the input space
lb = Bounds(1,3:4); %[0.1 0.05];
ub = Bounds(2,3:4); %[10.0 10.0];
 


Times = [];
controls.lb = lb;
controls.ub = ub;
controls.nstarts = 100;
controls.maxiter = 20;
controls.improv = 1e-6;

i1 = 1;
improv = inf;   % improvement between two successive query points
%demo_bayesoptimization
while ((i1 < controls.maxiter) && (improv > controls.improv))
    fprintf('**** BayOpt iter = %i***\n',i1)
    tic
    % Train the GP model for objective function and calculate variables
    % that are needed when calculating the Expected improvement
    % (Acquisition function) 
    if i1>1
        gp = gp_optim(gp,x,y);              %   GP_OPTIM  Optimize paramaters of a Gaussian process
                                            % MAP estimate for the parameters
                                            % gp_optim works as a wrapper for usual gradient based optimization functions
        [gpia,pth,th]=gp_ia(gp,x,y);        %   [GP_ARRAY, P_TH, TH, EF, VARF, PF, FF] = GP_IA(GP, X, Y, XT, OPTIONS)
                                            %    takes a GP structure GP with covariates X and observations Y
                                            %    and returns an array of GPs GP_ARRAY and corresponding weights
                                            %    P_TH and hyperparameter values.
        gp = gp_unpak(gp,sum(bsxfun(@times,pth,th)));
    end
    [K, C] = gp_trcov(gp,x);
    invC = inv(C);
    a = C\y;
    fmin = min( y );
%     fmin = min( fx(x) );
    
%     % Calculate EI and the posterior of the function for visualization
%     [Ef,Varf] = gp_pred(gp, x, y, xl);
%     EI = expectedimprovement_eg(xl, gp, x, a, invC, fmin);

    % optimize acquisition function
    %  * Note! Opposite to the standard notation we minimize negative Expected
    %    Improvement since Matlab optimizers seek for functions minimum
    %  * Note! We alternate the acquisition function between Expected
    %    Improvement and expected variance. The latter helps the
    %    optimization so that it does not get stuck in local mode
    % Here we use multiple starting points for the optimization so that we
    % don't crash into suboptimal mode of acquisition function
    if mod(i1,5)==0  % Do just exploration by finding the maximum variance location
        fh_eg = @(x_new) expectedvariance_eg(x_new, gp, x, [], invC); 
                        % expectedVariance_eg    Calculate the negative expected variance and its
                        %                        gradient  
    else
        fh_eg = @(x_new) expectedimprovement_eg(x_new, gp, x, a, invC, fmin);
    end
    indbest = find(y == fmin);
    xstart = [repmat(controls.lb,controls.nstarts,1) + repmat(controls.ub-controls.lb,controls.nstarts,1).*rand(controls.nstarts,P) ]; 
    x_news = zeros(controls.nstarts,P);
    for s1 = 1:controls.nstarts
        x_news(s1,:) = optimf(fh_eg, xstart(s1,:), [], [], [], [], controls.lb, controls.ub, [], opt);
    end
    EIs = expectedimprovement_eg(x_news, gp, x, a, invC, fmin);
    x_new = x_news( find(EIs==min(EIs),1), : ); % pick up the point where Expected Improvement is maximized


    [EFT, VARFT, LPYT, EYT, VARYT] = gp_pred(gp, x, y, x_new);
    % put new sample point to the list of evaluation points
    x(end+1,:) = x_new;
    
    %% CALL ABAQUS
    A = mpara.A;
    B = mpara.B;
%     Af = x_new(1)*mpara.Af; 
%     Bf = x_new(2)*mpara.Bf; 
    Af = x_new(1); 
    Bf = x_new(2); 
    As = mpara.As;
    Bs = mpara.Bs; 
    Afs = mpara.Afs;
    Bfs = mpara.Bfs;


    [LVVolumeAba, strainAbaTotal, strainComparison, SuccessB] = ...
                    Sim_LVPassiveForwardSimulationMatPres_NoGlobal(A,   B, ...
                                                      Af,  Bf, ...
                                                      As,  Bs, ...
                                                      Afs, Bfs, ...
                                                      pres, ...
                                                      false,...
                                                      options);

    ResLVVol(end+1,1) = LVVolumeAba;
    ResLVStrain(end+1,:) = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24);
    
    fval = [bsxfun(@minus, ResLVStrain(end,(strain_index_1:1:strain_index_24)), strainData(strain_index_1:1:strain_index_24)'),...
        (ResLVVol(end)-LVEDVMRI)/LVEDVMRI];
    fval = sum(fval.^2,2); 
    
    y(end+1,:) = fval;
    
    Times = [Times; toc];
    
    improv = abs(y(end) - y(end-1));
    i1 = i1 + 1;

    save(['BayOpt_af_bf_HV',HV,'.mat'],...
        'x','y','ResLVVol','ResLVStrain','i1','Times','gp','lik','cfse','controls')
end 

save(['BayOpt_af_bf_HV',HV,'.mat'],...
    'x','y','ResLVVol','ResLVStrain','i1','Times','gp','lik','cfse','controls')