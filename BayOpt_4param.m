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
% postprocess_only = false;
% regenerate_fibre = false;

%% Load the initial values based on LHSampling
fprintf('*** Using initial LHSampling *** \n')
load(['Init_LHS_4param_HV',HV,'.mat'],'ResLVVol','ResLVStrain','Theta','Times')
ind = (ResLVVol ~= 0);
ResLVVol = ResLVVol(ind,:);
ResLVStrain = ResLVStrain(ind,:);
x = Theta(ind,:);

y = Obj_fun(2,ResLVStrain,ResLVVol,strainData,LVEDVMRI,...
    strain_index_1,strain_index_24);
 
P = size(x,2);
%% Parameters and bounds
parnames = {'A','B','Af','Bf','As','Bs','Afs','Bfs'};
Bounds = [    0.1000    0.5000    0.1000    0.5000    0.1000    0.0500    0.1000    0.5000
             10.0000   30.0000   40.0000   30.0000   40.0000   30.0000   10.0000   30.0000];
Bounds_scale = Bounds(2,:) - Bounds(1,:);

%% construct GP to model the (ABAQUS) function
% cfc = gpcf_constant('constSigma2',10,'constSigma2_prior', prior_fixed);
% cfl = gpcf_linear('coeffSigma2', .01, 'coeffSigma2_prior', prior_sqrtt()); 
% cfl2 = gpcf_squared('coeffSigma2', .01, 'coeffSigma2_prior', prior_sqrtt(), 'interactions', 'on');
% cfse = gpcf_sexp('lengthScale',[5 5],'lengthScale_prior',prior_t('s2',4),'magnSigma2',.1,'magnSigma2_prior',prior_sqrtt('s2',10^2));

len_sc = [0.01 0.01 0.01 0.01];
cfse = gpcf_sexp('lengthScale',len_sc,'lengthScale_prior',prior_t('s2',4),'magnSigma2',.1,'magnSigma2_prior',prior_sqrtt('s2',10^2));
% cfse = gpcf_sexp('lengthScale',1,'magnSigma2',1,'magnSigma2_prior',prior_sqrtt('s2',10^2)); %GPCF_SEXP  Create a squared exponential (exponentiated quadratic) covariance function
lik = lik_gaussian('sigma2', 0.001, 'sigma2_prior', prior_fixed);
% gp = gp_set('cf', {cfc, cfl, cfl2, cfse}, 'lik', lik);
gp = gp_set('cf', {cfse}, 'lik', lik);
 

 
%% --- MCMC ---
disp(' MCMC integration over the parameters')
[gp_rec,gp_mcmc,opt_mcmc] = gp_mc(gp, x, y, 'nsamples', 1000,'display',20);
%    [RECORD, GP, OPT] = GP_MC(GP, X, Y, OPTIONS) Takes the Gaussian 
%    process structure GP, inputs X and outputs Y. Returns record 
%    structure RECORD with parameter samples, the Gaussian process GP
%    at current state of the sampler and an options structure OPT 
%    containing all the options in OPTIONS and information of the
%    current state of the sampler (e.g. the random number seed)

% After sampling we delete the burn-in and thin the sample chain
gp_rec = thin(gp_rec, 21, 2);
mn = zeros(1,5);

% subplot(3,3,1)
% hold on
% histogram(gp_rec.cf{1, 1}.magnSigma2,20)
mn(1) = mean(gp_rec.cf{1, 1}.magnSigma2);
% line([0.1 0.1],[0, 100],'color','red')
% line([mn(1), mn(1)],[0, 100],'color','blue')
% hold off
% title('magnSigma2')
for jj = 1:4
%     subplot(3,3,1+jj)
%     hold on
%     histogram(gp_rec.cf{1, 1}.lengthScale(:,jj),20)
    mn(jj+1) = mean(gp_rec.cf{1, 1}.lengthScale(:,jj));
%     line([len_sc(jj),len_sc(jj)],[0, 100],'color','red')
%     line([mn(jj+1), mn(jj+1)],[0, 100],'color','blue')
%     hold off
%     title(['lengthScale(',num2str(jj),')'])
end
% suptitle('100 draws, red - init, blue - post. mean')

% Make the predictions
% [Eft_mc, Varft_mc] = gp_pred(gp_rec, x, y, xt);

%% GP with "optimised parameters"
cfse = gpcf_sexp('lengthScale',mn(2:5),...
    'lengthScale_prior',prior_t('nu',4),... %prior_t('s2',4),...
    'magnSigma2',mn(1),...
    'magnSigma2_prior',prior_sqrtt('s2',10^2));
% cfse = gpcf_sexp('lengthScale',1,'magnSigma2',1,'magnSigma2_prior',prior_sqrtt('s2',10^2)); %GPCF_SEXP  Create a squared exponential (exponentiated quadratic) covariance function
lik = lik_gaussian('sigma2', 0.001, 'sigma2_prior', prior_fixed);
% gp = gp_set('cf', {cfc, cfl, cfl2, cfse}, 'lik', lik);
gp = gp_set('cf', {cfse}, 'lik', lik);


%% ----- conduct Bayesian optimization -----

% Set the options for optimizer of the acquisition function
optimf = @fmincon;
optdefault = struct('GradObj','on','LargeScale','off',...
    'Algorithm','trust-region-reflective','Display','Off',...
    'TolFun',1e-9,'TolX',1e-6);
opt = optimset(optdefault);
% lb = [0.1 0.1 0.1 0.1];     % lower bound of the input space <--------------????
% ub = [1 1 1 1];   % upper bound of the input space <--------------- ????
lb = [0.1 0.05 0.1 0.1];  
ub = [1 1 1 2]; 

if false
    HV = '25';
    load(['BayOpt_4param_HV',HV,'_oldBounds.mat'],'x','y','i1')
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


Times = [];
controls.lb = lb;
controls.ub = ub;
controls.nstarts = 500;
controls.maxiter = 101;
controls.improv = 1e-6;

i1 = 1;
improv = inf;   % improvement between two successive query points

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

%     xstart = [repmat(lb,controls.nstarts,1) + repmat(ub-lb,controls.nstarts,1).*rand(controls.nstarts,P) ]; 
        % Latin Hypercube Design for the new parameters parameters
    xstart = lhsdesign(controls.nstarts,P);
    xstart = bsxfun(@times, xstart, controls.ub - controls.lb);
    xstart = bsxfun(@plus, xstart, controls.lb);

    x_news = zeros(controls.nstarts,P);
    for s1 = 1:controls.nstarts
        x_news(s1,:) = optimf(fh_eg, xstart(s1,:), [], [], [], [], lb, ub, [], opt);
    end
    EIs = expectedimprovement_eg(x_news, gp, x, a, invC, fmin);
    x_new = x_news( find(EIs==min(EIs),1), : ); % pick up the point where Expected Improvement is maximized

    [EFT, VARFT, LPYT, EYT, VARYT] = gp_pred(gp, x, y, x_new);

    % put new sample point to the list of evaluation points
    x(end+1,:) = x_new;
    
    %% CALL ABAQUS
    A = x_new(1)*mpara.A;
    B = x_new(1)*mpara.B;
    Af = x_new(2)*mpara.Af; 
    Bf = x_new(3)*mpara.Bf; 
    As = x_new(2)*mpara.As;
    Bs = x_new(3)*mpara.Bs; 
    Afs = x_new(4)*mpara.Afs;
    Bfs = x_new(4)*mpara.Bfs;

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

    save(['BayOpt_4param_HV',HV,'.mat'],'Theta',...
        'x','y','ResLVVol','ResLVStrain','i1','Times','gp','lik','cfse','controls')
end 

save(['BayOpt_4param_HV',HV,'.mat'],'Theta',...
    'x','y','ResLVVol','ResLVStrain','i1','Times','gp','lik','cfse','controls')
