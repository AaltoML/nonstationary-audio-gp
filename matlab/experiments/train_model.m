clear; close all;
addpath('../');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals
addpath(genpath('gppad/'));
addpath(genpath('nmf/'));
addpath(genpath('toolsGP/'));
soundPath = '../../audio/speech/'; % Specify where to load the data from
% soundPath = '../../audio/textures/';
% soundPath = '../../audio/source_sep/training_data/';

  % load signal
  % Name of file to load:
  File = 'speech0_female';
%   File = 'speech1_male';
%   File = 'speech2_male';
%   File = 'stim312_wind';
%   File = '011PFNOM_M60_train';
%   File = '011PFNOM_M64_train';
%   File = '011PFNOM_M67_train';
%   File = '131EGLPM_M60_train';
%   File = '131EGLPM_M64_train';
%   File = '131EGLPM_M67_train';
%   File = '311CLNOM_M60_train';
%   File = '311CLNOM_M64_train';
%   File = '311CLNOM_M67_train';
%   File = 'ALVARADO_M60_train';
%   File = 'ALVARADO_M64_train';
%   File = 'ALVARADO_M67_train';
  fs_ = 16000; % sampling rate of file
  
  D = 16; % number of frequency channels
  N = 3; % number of NMF components
  dt = 1;
  kernel1 = 'matern32'; % kernel for subbands
  kernel2 = 'matern52'; % kernel for amplitude envelopes
%   mod_sparsity_ = 1.; % 0 = softplus
  link = @(g) log(1+exp(g-1)); % link function
%   link_ = @(g, mod_sparsity_) 0.5*log(1+exp(g-mod_sparsity_)); % link function
%   link = @(g) 0.5*exp(g);
  rng(12345,'twister')
  w_lik_train = 1e-3; % observation noise
  w_lik_test = 1e-3;
  
  inf_type = 'EP'; % 'EP', 'IHGP_EP', 'EKF'
  
  %%% tune the hyperparameters? %%%
  optimise = 1;
    % training params
    p_cubature = 7; % order of cubature for Gaussian integral
    max_iters = 15; % maximum number of iterations
    
  %   likfunc = @likModulatorNMFPower;
%   likfunc = @likModulatorNMFPowerSq; % models the square of the amplitudes (i.e the spectrogram)
  likfunc = @likModulatorPreCalcwn; % models the square, and takes pre-calculated sigma pts
  
  ep_fraction = 0.1; % Power EP fraction (0.8-1 seems best)
  
  %%% training params %%%
  ep_itts_train = 1; % EP iterations during training
  ep_damping_train = linspace(0.01, 0.1, ep_itts_train); % EP damping during training
  
  %%% testing params %%%
  ep_itts_test = 20; % EP iterations during testing
  ep_damping_test = linspace(0.01, 0.1, ep_itts_test); % EP damping during testing (0.1 seems best)
  
  l_iter = 1; % local EKF (inner loop) iterations
  
  train_duration = 16000;
  
  
%% Load signal and pre-process
  [y_,fs] = audioread([soundPath,File,'.wav']); % reads in the file
  y = resample(y_, fs_, fs); % downsample the input
  fs = fs_;
  normaliser = sqrt(var(y));
  y_norm = y/normaliser; % rescale the input to unit variance
  T = length(y_norm);
  t = linspace(1,T,T)';
  yTrain = y_norm(1:train_duration);
  TTrain = length(yTrain);
  tTrain = linspace(1,TTrain,TTrain)';
  
  
%% Pre-optimise the filter bank with exponential kernel (probabilistic phase vocoder)
  
  % Learn properties of the filters (centre frequency and width)
  opts = struct;
  opts.verbose = 0; % view plots of the fitting process
  opts.minT = 100;
  opts.maxT = 1000;
  opts.numIts = 10;
  opts.numLevels = 30;
  opts.bet = 750;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together
  
  [varx,lamx,om] = fit_probSTFT_SD(y_norm,D,'exp',opts); % trains filters to match the spectrum
  [omega, om_ind] = sort(om);
  lenx = 1 ./ lamx(om_ind);
  varx = varx(om_ind);
  
  [A1,Q1,H1,Pinf1,K1,tau1] = get_disc_model(lamx,varx,omega,D,'exp',6);
  
  
%% Run GP prob amp demodulation
  ZTest1 = kernel_ss_probFB(y_norm,A1,Q1,H1,Pinf1,K1,0,tau1);
  ATest1 = abs(ZTest1');
  
  disp('Running GPPAD to get initialisation for smoothed envelopes')
  % We're interested in NMF lengthscales on the order of ~100ms and longer,
  % anything shorter than this can be soaked up in the subband variation.
  len = fs/10;
  [mods,cars] = GPPAD(real(ZTest1)',len);
  
  
%% Initialise with normal NMF
  HInit = exp(randn(T,N));
  ks = ceil(T*rand(N,1));
  WInit = mods(ks,:);
  vary = zeros(T,D);
  Opts.restarts = 20;
  Opts.numIts = 500;
  %[WEst1,HEst1,info1] = nmf(ATrain,WInit,HInit,[],[],[],vary,Opts);
  [WEst,HEst,info1] = nmf_fp(mods,WInit,HInit,vary,Opts);

  % order according to slowness (mean square derivative)
  fastness = mean(diff(HEst).^2)./var(HEst);
  [val,ind] = sort(fastness,'descend');
  HEst = HEst(:,ind); 
  WEst = WEst(ind,:);
%%
  lenx_nmf = zeros(N,1);
  mux = zeros(N,1);
  varx_nmf = zeros(N,1);

  for k=1:N
    % threshold
    logHthresh = log(exp(HEst(:,k)+1e-8)-1);
  
    % filter
    filt = log(1+exp(-1/2*([-100:100].^2)/(1^2)));
    filt = filt/sum(filt);
    logHsm = conv(logHthresh,filt,'same');
  
    % fit GP
    mux(k) = mean(logHsm);
    [lenx_nmf(k),varx_nmf(k),info] = trainSEGP_RS(logHsm-mux(k));
  
  end
  
  % shift the link function to alter sparsity
%   mod_sparsity = -max(mux);
%   link = @(g) link_(g, mod_sparsity);
  
  
%% Now set up the full model  
  % define constraints
  range_var_fast = [0.01,  0.1]; % subband variance should be fixed
  range_len_fast = [ 100, 2000]; % subband lengthscales should be short
  range_omega    = [   0, 2*pi]; % subband centre frequencies
  range_var_slow = [   1,   10]; % NMF modulator variances
  range_len_slow = [ 200, 5000]; % NMF modulator lengthscales should be long
  range_W        = [   0,  1.5]; % NMF weights
  
  constraints = [range_var_fast; range_len_fast; range_omega; ...
                 range_var_slow; range_len_slow; range_W];
  init_param = @(constraint,D_,N_) constraint(1) + (constraint(end)-constraint(1))*rand(D_,N_);
  
  % which parameters to optimise:
  %              w_lik, var_f, len_f, omega, var_s, len_s, W(nmf)
  tune_hypers = [0,     0,     1,     0,     1,     1,     1     ];
  
  
  % pick model parameters
  var_fast = mean(range_var_fast).*ones(D,1);
  len_fast = min(max(1.*lenx,range_len_fast(1)+1),range_len_fast(end)-1);
%   omega = linspace(pi/4,pi/50,D)';
  var_slow = min(max(1.*varx_nmf,range_var_slow(1)+1),range_var_slow(end)-1);
  len_slow = min(max(5.*lenx_nmf,range_len_slow(1)+1),range_len_slow(end)-1);
  W = WEst';
  
  
%% setup and optimisation
  
  % pre-calculate the sigma-point weights
  [wn,xn_unscaled,~] = utp_ws(p_cubature,N);

  % Moments
%   mom = @(hyp,mu,s2,nmfW,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,p_cubature,ep_fraction,'infEP');
  mom = @(hyp,mu,s2,nmfW,ep_frac,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,ep_frac,wn,xn_unscaled,'infEP');
  
  
  % State space model
  ss = @(x,p1,p2,kern1,kern2) ss_modulators_nmf(p1,p2,kern1,kern2);
  
  num_lik_params = length(w_lik_train);
  
  % Hyperparameters initial
  w_all = {log(w_lik_train); inv_sigmoid(var_fast,range_var_fast); inv_sigmoid(len_fast,range_len_fast);...
           inv_sigmoid(omega,range_omega); inv_sigmoid(var_slow,range_var_slow); ...
           inv_sigmoid(len_slow,range_len_slow); inv_sigmoid(W(:),range_W)};
       
  w_init = w_all(tune_hypers>0); w_fixed = w_all(tune_hypers==0);
  w_init = cell2mat(w_init);
  w_fixed = cell2mat(w_fixed);
  
  yTrain_ = yTrain;
  
  %{
  trainGapLen = 100;
  trainGapPos = 500:1000:train_duration-500;%[500,1500,4500,6000,7500];
  NTrainGaps = length(trainGapPos);
  trainGaps = ceil(linspace(gapLim(1),gapLim(2),numgaps));
  ind = [];
  for ng=1:NTrainGaps
    ind = [ind,trainGapPos(ng)+[-ceil(trainGapLen/2):+ceil(trainGapLen/2)]];
  end  
%   yTrain_(ind) = NaN;
  %}
  
  
  if optimise
    disp('Now we finally tune the parameters with our model:')
    tic
    % Optimization options
    opts = optimset('GradObj','off','display','iter','DerivativeCheck','off','MaxIter',max_iters);
    
    % Optimize hyperparameters w.r.t. log marginal likelihood
    switch inf_type
      case 'IHGP_EP'
        [w_opt,ll] = fminunc(@(w_init) ihgp_adf_modulator_nmf_constraints(w_init,tTrain,yTrain_,ss,mom,[],kernel1,kernel2,num_lik_params,D,N, ...
                                                              ep_fraction,ep_damping_train,ep_itts_train,constraints,w_fixed,tune_hypers), ...
                         w_init,opts);
      case 'EP'
        [w_opt,ll] = fminunc(@(w_init) gf_ep_modulator_nmf_constraints(w_init,tTrain,yTrain_,ss,mom,[],kernel1,kernel2,num_lik_params,D,N, ...
                                                              ep_fraction,ep_damping_train,ep_itts_train,constraints,w_fixed,tune_hypers), ...
                         w_init,opts);
      case 'EKF'
        [w_opt,ll] = fminunc(@(w_init) gf_giekf_modulator_nmf_constraints(w_init,tTrain,yTrain_,ss,mom,[],kernel1,kernel2,num_lik_params,D,N, ...
                                                              ep_itts_train,l_iter,constraints,w_fixed,tune_hypers,opts.GradObj), ...
                         w_init,opts);
      otherwise
        error('inf_type not recognised');
    end
    toc              
  end
  
  
  
%% Saving
       
  % unwrap the optimisable and fixed parameters with constraints
  w_ind = 1; wf_ind = 1;  
  if tune_hypers(1)
    lik_param = exp(w_opt(1:num_lik_params));
    w_ind = w_ind + num_lik_params;
  else
    lik_param = exp(w_fixed(1:num_lik_params));
    wf_ind = wf_ind + num_lik_params;
  end 
  param1 = [];
  param2 = [];
  for i=2:6
    if tune_hypers(i)
      if i<=4
        param1 = [param1; sigmoid(w_opt(w_ind:w_ind+D-1),constraints(i-1,:))];
        w_ind = w_ind + D;
      else
        param2 = [param2; sigmoid(w_opt(w_ind:w_ind+N-1),constraints(i-1,:))];
        w_ind = w_ind + N;
      end
    else
      if i<=4
        param1 = [param1; sigmoid(w_fixed(wf_ind:wf_ind+D-1),constraints(i-1,:))];
        wf_ind = wf_ind + D;
      else
        param2 = [param2; sigmoid(w_fixed(wf_ind:wf_ind+N-1),constraints(i-1,:))];
        wf_ind = wf_ind + N;
      end
    end
  end
  
  if tune_hypers(7)
    Wnmf = reshape(sigmoid(w_opt(w_ind(1):end),constraints(6,:)),[D,N]);
  else
    Wnmf = reshape(sigmoid(w_fixed(wf_ind(1):end),constraints(6,:)),[D,N]);
  end
  
  results = struct;
  results.File = File;
  results.soundPath = soundPath;
  results.fs = fs;
  results.y_norm = y_norm;
  results.normaliser = normaliser;
  results.D = D; % number of frequency channels
  results.N = N; % number of NMF components
  results.kernel1 = kernel1; % kernel for subbands
  results.kernel2 = kernel2; % kernel for amplitude envelopes
  results.link = link; % link function
  results.w_lik = lik_param;
  results.param1 = param1;
  results.param2 = param2;
  results.Wnmf = Wnmf;
  results.inf_type = inf_type; % 'IHGP_EP', 'EP' or 'EKF'
  results.ep_fraction = ep_fraction; % Power EP fraction
  results.ep_damping_train = ep_damping_train; % EP damping
  results.ep_damping_test = ep_damping_test; % EP damping
  results.ep_itts_train = ep_itts_train;
  results.ep_itts_test = ep_itts_test;
  results.train_duration = train_duration;
  results.max_iters = max_iters;
  results.w_init = w_init;
  results.w_fixed = w_fixed;
  results.w_opt = w_opt;
  results.constraints = constraints;
  results.tune_hypers = tune_hypers;
  matName = strcat('results_',File,'_',inf_type,'.mat');
  save(matName,'-struct','results');
  