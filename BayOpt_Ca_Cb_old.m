%% Initially prepared by Hao Gao
%% here we will optimize the volume according to the klotz curve
clear all; close all; clc;

addpath(genpath('../../../GPstuff'));

%% load gloabal options for HV1. Useful data will be saved in globla option structure
HV = '25';
config_file = ['Sim_DirConfigForwardComputation_HV',HV];
% Sim_DirConfigForwardComputation_HV02;
eval(config_file);


%% read the strain data from images and LV measured end-diastolic volume
load(['DataMRI_HV',HV,'.mat'], 'HV','LVEDVMRI','strainData','strainDataFlag',...
    'strain_index_1','strain_index_24');

%% using the right parameter sets
%% A, B, Af, Bf, As, Bs, Afs, Bfs are the 8 unknown parameters for HGO
%% myocardial material model
% A = 0.224487;
% B = 1.621500;
% Af = 2.426717;
% Bf = 1.826862;
% As = 0.556238;
% Bs = 0.774678;
% Afs = 0.390516;
% Bfs = 1.695000;


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
postprocess_only = false;
regenerate_fibre = false;


%% before calling the optimization function, we will need to update the global options
load(['Sim_FS_LV_Ca_Cb_HV',HV,'.mat'],'ResLVVol','ResLVStrain','SimRes','Ca','Cb','mpara','Times')

% use the Ca-Cb pair leading to the smallest value of the objective
% function
fO1_HG = Obj_fun(1,ResLVStrain,ResLVVol,strainData,LVEDVMRI);
fO2_HG = Obj_fun(2,ResLVStrain,ResLVVol,strainData,LVEDVMRI);
[fO1_min, ind_min] = min(fO1_HG);
ind_b = mod(ind_min,length(Cb));
ind_a = (ind_min-ind_b)/length(Cb);
Ca = Ca(ind_a+1); 
Cb = Cb(ind_b);


% Ca = 0.4;
% Cb = 0.3;
mpara.A = 0.23619*Ca;
mpara.B = 10.810*Cb;
mpara.Af = 20.037*Ca;
mpara.Bf = 14.154*Cb;
mpara.As = 3.7245*Ca;
mpara.Bs = 5.1645*Cb;
mpara.Afs = 0.41088*Ca;
mpara.Bfs = 11.300*Cb;
options.mpara = mpara; %%remember to update here


options.strainData= strainData;
options.LVEDVMRI = LVEDVMRI;
options.LVEDP_high = 30;
options.LVEDP_norm = 8;


Lb_aafs = [0.1 0.05];
Ub_aafs = [10.0 10.0];
options_Opt = optimset('Algorithm', 'sqp', 'TolFun', 1e-3, ...
                       'TolX',0.001,'Diffminchange',1e-3);
ca_cb_0 = [1.0 1.0];
tic
[ca_cb,fval,exitflag,output] = fmincon(@LVPassive_model_optimization_fmincon_ca_cb_klotz,...
    ca_cb_0,[],[], [], [], Lb_aafs,Ub_aafs,[],options_Opt);
Times = toc;

save(['BayOpt_Ca_Cb_Klotz_HV',HV,'.mat'],'Ca','Cb','mpara','Times')



%% The objective function
% fx = @(x) -log( (mvnpdf([x(:,1) x(:,2)],[-1.5 -2.5], [1 0.3; 0.3 1]) + 0.3*mvnpdf([x(:,1) x(:,2)],[2 3], [3 0.5; 0.5 4])).*...
%     mvnpdf([x(:,1) x(:,2)],[0 0], [100 0; 0 100])) ./15 -1;
 
fx = @(xx)  LVPassive_model_optimization_fmincon_ca_cb_klotz(xx);
 
%% Load the initial values based on LHSampling
load('Init_LHS_Ca_Cb_HV25.mat')
scatter3(Cab(:,1),Cab(:,2),y)


%% construct GP to model the (ABAQUS) function
% cfc = gpcf_constant('constSigma2',10,'constSigma2_prior', prior_fixed);
% cfl = gpcf_linear('coeffSigma2', .01, 'coeffSigma2_prior', prior_sqrtt()); 
% cfl2 = gpcf_squared('coeffSigma2', .01, 'coeffSigma2_prior', prior_sqrtt(), 'interactions', 'on');
% cfse = gpcf_sexp('lengthScale',[5 5],'lengthScale_prior',prior_t('s2',4),'magnSigma2',.1,'magnSigma2_prior',prior_sqrtt('s2',10^2));
cfse = gpcf_sexp('lengthScale',[0.01 0.01],'lengthScale_prior',prior_t('s2',4),'magnSigma2',.1,'magnSigma2_prior',prior_sqrtt('s2',10^2));
% cfse = gpcf_sexp('lengthScale',1,'magnSigma2',1,'magnSigma2_prior',prior_sqrtt('s2',10^2)); %GPCF_SEXP  Create a squared exponential (exponentiated quadratic) covariance function
lik = lik_gaussian('sigma2', 0.001, 'sigma2_prior', prior_fixed);
% gp = gp_set('cf', {cfc, cfl, cfl2, cfse}, 'lik', lik);
gp = gp_set('cf', {cfse}, 'lik', lik);
 

%% ----- conduct Bayesian optimization -----

% Set the options for optimizer of the acquisition function
optimf = @fmincon;
optdefault = struct('GradObj','on','LargeScale','off','Algorithm','trust-region-reflective','TolFun',1e-9,'TolX',1e-6);
opt = optimset(optdefault);
% lb=[0 0];     % lower bound of the input space
% ub=[1 1];   % upper bound of the input space
% Lb_aafs = [0.1 0.05];
% Ub_aafs = [10.0 10.0];
lb = [0.1 0.05];
up = [10.0 10.0];

[X,Y] = meshgrid(linspace(lb(1),ub(1),100),linspace(lb(2),ub(2),100));
xl = [X(:) Y(:)];
% draw initial points
% x = [-4 -4;-4 4;4 -4;4 4;0 0];
% y = fx(x);
x = Cab;

% figure, % figure for visualization
i1 = 1;
maxiter = 20;
improv = inf;   % improvement between two successive query points
while i1 < maxiter && improv>1e-6
    % Train the GP model for objective function and calculate variables
    % that are needed when calculating the Expected improvement
    % (Acquisition function) 
    if i1>1
        gp = gp_optim(gp,x,y);              %   GP_OPTIM  Optimize paramaters of a Gaussian process 
        [gpia,pth,th]=gp_ia(gp,x,y);        %   [GP_ARRAY, P_TH, TH, EF, VARF, PF, FF] = GP_IA(GP, X, Y, XT, OPTIONS)
                                            %    takes a GP structure GP with covariates X and observations Y
                                            %    and returns an array of GPs GP_ARRAY and corresponding weights
                                            %    P_TH and hyperparameter values.
        gp = gp_unpak(gp,sum(bsxfun(@times,pth,th)));
    end
    [K, C] = gp_trcov(gp,x);
    invC = inv(C);
    a = C\y;
%     fmin = min( y );
    fmin = min( fx(x) );
    
    % Calculate EI and the posterior of the function for visualization
    [Ef,Varf] = gp_pred(gp, x, y, xl);
    EI = expectedimprovement_eg(xl, gp, x, a, invC, fmin);

    % optimize acquisition function
    %  * Note! Opposite to the standard notation we minimize negative Expected
    %    Improvement since Matlab optimizers seek for functions minimum
    %  * Note! We alternate the acquisition function between Expected
    %    Improvement and expected variance. The latter helps the
    %    optimization so that it does not get stuck in local mode
    % Here we use multiple starting points for the optimization so that we
    % don't crash into suboptimal mode of acquisition function
    if mod(i1,5)==0  % Do just exploration by finding the maimum variance location
        fh_eg = @(x_new) expectedvariance_eg(x_new, gp, x, [], invC); 
                        % expectedVariance_eg    Calculate the negative expected variance and its
                        %                        gradient  
    else
        fh_eg = @(x_new) expectedimprovement_eg(x_new, gp, x, a, invC, fmin);
    end
    indbest = find(y == fmin);
    nstarts = 20;
    xstart = [repmat(lb,nstarts,1) + repmat(ub-lb,nstarts,1).*rand(nstarts,2) ]; 
    for s1=1:length(xstart)
        x_new(s1,:) = optimf(fh_eg, xstart(s1,:), [], [], [], [], lb, ub, [], opt);
    end
    xnews = x_new;
    EIs = expectedimprovement_eg(x_new, gp, x, a, invC, fmin);
    x_new = x_new( find(EIs==min(EIs),1), : ); % pick up the point where Expected Improvement is maximized
    % x_new = [0 1] for lengthsacel = [5 5] or [1 1]  
    % x_new = [0.0757 0.1954] --> EYT =    -2.9858
    % x_new = [0.1067 0.8530] --> EYT =    4.1979e-19

    [EFT, VARFT, LPYT, EYT, VARYT] = gp_pred(gp, x, y, x_new);
%     scatter3(x(:,1),x(:,2),y)
%     hold on
%     scatter3(x_new(:,1),x_new(:,2),EYT,'MarkerEdgeColor','r')    
%     hold off
    
    % put new sample point to the list of evaluation points
    x(end+1,:) = x_new;
    
    % CALL ABAQUS
    fx = 
    y(end+1,:) = fx(x_new);     % calculate the function value at query point

    if false
        % visualize
        clf
        % Plot the objective function
        subplot(2,2,1),hold on, title('Objective, query points')
        box on
        pcolor(X,Y,Z),shading flat
        clim = caxis;
        l1=plot(x(1:end-1,1),x(1:end-1,2), 'rx', 'MarkerSize', 10);
        %plot(x(end,1),x(end,2), 'ro', 'MarkerSize', 10)
        l2=plot(xnews(:,1),xnews(:,2), 'ro', 'MarkerSize', 10);
        l3=plot(x(end,1),x(end,2), 'ro', 'MarkerSize', 10, 'linewidth', 3);
        legend([l1,l2,l3], {'function evaluation points','local modes of acquisition function','The next query point'})
        % Plot the posterior mean of the GP model for the objective function
        subplot(2,2,2),hold on, title(sprintf('GP prediction, mean, iter: %d',i1))
        box on
        pcolor(X,Y,reshape(Ef,100,100)),shading flat
        caxis(clim)
        % Plot the posterior variance of GP model
        subplot(2,2,4),hold on, title('GP prediction, variance')
        box on
        pcolor(X,Y,reshape(Varf,100,100)),shading flat
        l2=plot(xnews(:,1),xnews(:,2), 'ro', 'MarkerSize', 10);
        l3=plot(x(end,1),x(end,2), 'ro', 'MarkerSize', 10, 'linewidth', 3);
        % Plot the expected improvement 
        subplot(2,2,3), hold on, title(sprintf('Expected improvement %.2e', min(EIs)))
        box on
        pcolor(X,Y,reshape(EI,100,100)),shading flat
        plot(xnews(:,1),xnews(:,2), 'ro', 'MarkerSize', 10);
        plot(x(end,1),x(end,2), 'ro', 'MarkerSize', 10, 'linewidth', 3);


        improv = abs(y(end) - y(end-1));
        i1=i1+1;
        pause
    end
end 