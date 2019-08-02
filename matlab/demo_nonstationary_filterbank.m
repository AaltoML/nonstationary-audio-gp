clear; close all;
addpath('prob_filterbank/'); % filterbank code
addpath('symmetric-cubature-rules/'); % code for approximate Gaussian integrals
soundPath = '../audio/speech/'; % Specify where to load the data from
% soundPath = '../audio/textures/';

  % load signal
  File = 'speech0_female'; % Name of file to load
%   File = 'stim312_wind'; % Name of file to load
  fs_ = 24000; % sampling rate of file
  
  D = 12; % number of frequency channels
  N = 3; % number of NMF components
  dt = 1;
  kernel1 = 'exp'; % kernel for subbands
  kernel2 = 'matern52'; % kernel for amplitude envelopes
  T = 24000; % max length of signal
  t = linspace(1,T,T)';
  mod_sparsity = 0.; % 0 = softplus
  link = @(g) log(1+exp(g-mod_sparsity)); % link function
%   link = @(g) exp(g);
  rng(12345)
  w_lik = 1e-3; % observation noise
  p_cubature = 9; % order of cubature for Gaussian integral
  max_iters = 5; % maximum number of iterations
  
  %%% tune the hyperparameters? %%% (optimisation still quite unstable)
  optimise = 0;
  
  ep_fraction = 0.5; % Power EP fraction
  ep_itts = 3; % EP iterations
  ep_damping = linspace(0.5, 0.5, ep_itts); % EP damping
  
  
%% Load signal and pre-process
  [y,fs] = audioread([soundPath,File,'.wav']); % reads in the file
  yTest = resample(y, fs_, fs); % downsample the input
  fs = fs_;
  normaliser = sqrt(var(yTest));
  yTest = yTest/normaliser; % rescale the input to unit variance
  
  yTest = yTest(4501:4501+T-1);
  yTrain = yTest;
  
%% Pre-optimise the filter bank with exponential kernel (probabilistic phase vocoder)
  % Learn properties of the filters (centre frequency and width)
  opts = struct;
  opts.verbose = 1; % view plots of the fitting process
  opts.minT = 100;
  opts.maxT = 1000;
  opts.numIts = 10;
  opts.numLevels = 30;
  opts.bet = 750;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together
  
  [varx,lamx,om] = fit_probSTFT_SD(yTrain,D,'exp',opts); % trains filters to match the spectrum
  [omega, om_ind] = sort(om);
  lenx = 1 ./ lamx(om_ind);
  varx = varx(om_ind);
  
  [A1,Q1,H1,Pinf1,K1,tau1] = get_disc_model(lamx,varx,omega,D,'exp',6);

if ~strcmp(kernel1,'exp')
% Then tune parameters with new kernel
  % Learn properties of the filters (centre frequency and width)
  opts = struct;
  opts.verbose = 1; % view plots of the fitting process
  opts.minT = 500;
  opts.maxT = 1000;
  opts.numIts = 10;
  opts.numLevels = 20;
  opts.bet = 750;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together
  opts.bandwidth_lim = 5;  % limits the bandwidths to n times the initial values
  
  opts.theta_init = [varx;lamx;omega];
  [varx,lamx,om] = fit_probSTFT_SD(yTrain,D,kernel1,opts); % trains filters to match the spectrum
  [omega, om_ind] = sort(om);
  if strcmp(kernel1,'matern32')
    lenx = sqrt(3) ./ lamx(om_ind);
  elseif strcmp(kernel1,'matern52')
    lenx = sqrt(3) ./ lamx(om_ind);
  else
    error('cant calculate lengthscales')
  end
  varx = varx(om_ind);
  
  [A1,Q1,H1,Pinf1,K1,tau1] = get_disc_model(lamx,varx,omega,D,'exp',6);
end
  
  
%%
  ZTest1 = kernel_ss_probFB(yTest,A1,Q1,H1,Pinf1,K1,0,tau1);
  ATest1 = abs(ZTest1');
  [Hnmf,Wnmf] = nnmf(ATest1,N);
%   Wnmf = Wnmf' / (D/10 * max(Wnmf(:)));
  Wnmf = Wnmf' / max(Wnmf(:));
%   Wnmf = max(Wnmf,0.1);
  
  
%% Now set up the full model  
  % define some constraints
  range_var_fast = [0.01,  0.1]; % subband variance should be fixed
  range_len_fast = [  20, 2000]; % subband lengthscales should be short
  range_omega    = [1e-3, 2*pi]; % subband centre frequencies
  range_var_slow = [   2,    5]; % NMF modulator variances
  range_len_slow = [ 300, 5000]; % NMF modulator lengthscales should be long
  range_W        = [1e-3, 1.25]; % NMF weights
  
  constraints = [range_var_fast; range_len_fast; range_omega; ...
                 range_var_slow; range_len_slow; range_W];
  init_param = @(constraint,D_,N_) constraint(1) + (constraint(end)-constraint(1))*rand(D_,N_);
  
  % which parameters to optimise:
  %              w_lik, var_f, len_f, omega, var_s, len_s, W(nmf)
  tune_hypers = [1,     1,     1,     1,     1,     1,     1     ];
  
  
  % pick some random model parameters
  var_fast = mean(range_var_fast)*ones(D,1);
  len_fast = min(max(lenx,range_len_fast(1)),range_len_fast(end));
%   var_fast = varx;
%   omega = linspace(pi/4,pi/50,D)';
  var_slow = mean(range_var_slow)*ones(N,1);
  len_slow = linspace(range_len_slow(1)+2,range_len_slow(end)-2,N)';
  W = Wnmf;%init_param(constraints(6,:),D,N);
  
  
%% Prior predictive checks
  % State space sampling
  w_sub = [var_fast;len_fast;omega];
  w_mod = [var_slow;len_slow];
  [F,L,Qc,H,Pinf] = ss_modulators_nmf(w_sub,w_mod,kernel1,kernel2);
  [A,Q] = lti_disc(F,L,Qc,dt);
  
  z = zeros(length(L),T);
  y_samp = zeros(T,1); y_f = zeros(T,D); y_s = zeros(T,N);
  z(:,t(1)) = chol(Pinf,'lower')' * randn(length(L),1);
  y_samp(t(1)) = (H(1:D,:) * z(:,t(1)))' * W * H(D+1:D+N,:) * link(z(:,t(1)));
  y_f(t(1),:) = H(1:D,:) * z(:,t(1)); y_s(t(1),:) = H(D+1:D+N,:) * z(:,t(1));
  for k=2:T
    z(:,t(k)) = A * z(:,t(k-1)) + chol(Q)' * randn(length(L),1); % prior
    y_f(t(k),:) = H(1:D,:) * z(:,t(k)); y_s(t(k),:) = H(D+1:D+N,:) * z(:,t(k)); % component ground truth
    y_samp(t(k)) = (H(1:D,:) * z(:,t(k)))' * W * H(D+1:D+N,:) * link(z(:,t(k))); % likelihood
  end
  
  mods = link(y_s) * W';
  figure(1);clf
  for i=1:D
    subplot(D+1,2,2*i-1)
    plot(t/fs,y_f(:,i))
    hold on
    plot(t/fs,mods(:,i))
    if i==1
      legend('fast varying signal','slow varying modulator','fontsize',10)
      title('Prior predictive checks')
    end
    if i<=N
      subplot(D+1,2,2*i)
      plot(t/fs,link(y_s(:,i)),'r')
    end
  end
  subplot(D+1,1,D+1)
  plot(t/fs,y_samp,'g')
  legend('sum product','fontsize',10)
  
  
%% setup and optimisation
  
  likfunc = @likModulatorNMFPower;
  
  l_iter = 1; % local EKF (inner loop) iterations
  
  % Moments
  mom = @(hyp,mu,s2,nmfW,ep_frac,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,p_cubature,ep_frac,'infEP');
  
  % State space model
  ss = @(x,p1,p2,kern1,kern2) ss_modulators_nmf(p1,p2,kern1,kern2);
  
  num_lik_params = length(w_lik);
  
  % Hyperparameters initial
  w = log([w_lik;...
             var_fast;len_fast;omega;...
             var_slow;len_slow;...
             W(:)]);
  
  if optimise
    tic
    ep_itts_opt = 1;
    % Optimization options
    opts = optimset('GradObj','off','display','iter','DerivativeCheck','off','MaxIter',max_iters);
    
    % Optimize hyperparameters w.r.t. log marginal likelihood
    [w_opt,ll] = fminunc(@(w) gf_ep_modulator_nmf(w,t,yTest,ss,mom,[],kernel1,kernel2,num_lik_params,D,N,ep_fraction,ep_damping,ep_itts_opt),w,opts);
    % other inference methods:
%     [w_opt,ll] = fminunc(@(w) gf_giekf_modulator_nmf(w,t,yTest,ss,mom,[],kernel1,kernel2,num_lik_params,D,N,ep_itts_opt,l_iter),w,opts);
%     [w_opt,ll] = fminunc(@(w) ihgp_ep_modulator_nmf(w,t,yTest,ss,mom,[],kernel1,kernel2,num_lik_params,D,N,ep_fraction,ep_damping,ep_itts_opt),w,opts);
    toc
  end
  
  
%% inference and plotting
  
  if exist('w_opt','var'), w_ = w_opt; else, w_ = w; end
  % ADF filtering
  tic
   [Eft,Varft,Covft,lb,ub,out] = gf_ep_modulator_nmf(w_,t,yTest,ss,mom,t,kernel1,kernel2,num_lik_params,D,N,ep_fraction,ep_damping,ep_itts);
   % other inference methods:
%    [Eft,Varft,Covft,lb,ub,out] = gf_giekf_modulator_nmf(w_,t,yTest,ss,mom,t,kernel1,kernel2,num_lik_params,D,N,ep_itts,l_iter);
%    [Eft,Varft,Covft,lb,ub,out] = ihgp_ep_modulator_nmf(w_,t,yTest,ss,mom,t,kernel1,kernel2,num_lik_params,D,N,ep_fraction,ep_damping,ep_itts);
   
  toc
  
  if tune_hypers(end)
    W_ = reshape(exp(w_(end-D*N+1:end)),[D,N]);
  else
    W_ = reshape(exp(w_fixed(end-D*N+1:end)),[D,N]);
  end
  
%   W_ = W; % fixed weights
  
  Eft_mod = zeros(size(Eft(D+1:end,:)));
  Varft_mod = zeros(size(Varft(D+1:end,:)));
  s = 250;
  sub_samp = zeros([size(Eft(1:D,:)),s]);
  mod_samp = zeros([size(Eft(D+1:end,:)),s]);
  
  figure(2);clf
  for i=1:D
    sub_samp(i,:,:) = bsxfun(@plus,bsxfun(@times,randn(T,s),sqrt(Varft(i,:))'),Eft(i,:)');
    subplot(D+1,2,2*i-1)
    hold on
    plot(Eft(i,:),'r')
    plot(lb(i,:),'g--')
    plot(ub(i,:),'g--')
    if i==1
        legend('posterior mean','95% confidence','fontsize',10)
        title('subbands')
    end
    if i<=N
      mod_samp(i,:,:) = bsxfun(@plus,bsxfun(@times,randn(T,s),sqrt(Varft(D+i,:))'),Eft(D+i,:)');
      Eft_mod(i,:) = mean(link(mod_samp(i,:,:)),3);
      Varft_mod(i,:) = var(link(mod_samp(i,:,:)),0,3);
      subplot(D+1,2,2*i)
      hold on
      plot(Eft_mod(i,:),'r')
      plot(Eft_mod(i,:)+1.96*sqrt(Varft_mod(i,:)),'g--')
      plot(Eft_mod(i,:)-1.96*sqrt(Varft_mod(i,:)),'g--')
      if i==1, title('NMF components'), end
    end
  end
  sig_samp = zeros(T,s);
  for v=1:s
    sig_samp(:,v) = sum((W_ * link(mod_samp(:,:,v))) .* sub_samp(:,:,v),1);
  end
  Esig = mean(sig_samp,2);
  Vsig = var(sig_samp,0,2);
  subplot(D+1,1,D+1)
  plot(yTest)
  hold on
  plot(Esig,'r')
  plot(Esig + 1.96*sqrt(Vsig),'g--')
  plot(Esig - 1.96*sqrt(Vsig),'g--')
  title('signal (via NMF)')
  
  figure(3);clf
  plot(yTest)
  hold on
  plot(Esig,'r')
  plot(Esig + 1.96*sqrt(Vsig),'g--')
  plot(Esig - 1.96*sqrt(Vsig),'g--')
  legend('ground truth','posterior mean','95% confidence','fontsize',10)
  title('signal (via NMF)')
  
  fprintf('RMSE: %d \n',sqrt(mean((yTest-Esig).^2)))
  fprintf('lZ: %d \n',sum(out.lZ))
  
  var_fast_opt = exp(w_(2:D+1));
  len_fast_opt = exp(w_(D+2:2*D+1));
  omega_opt = exp(w_(2*D+2:3*D+1));
  switch kernel1
    case 'exp'
      lam_fast_opt = 1./len_fast_opt;
    case 'matern32'
      lam_fast_opt = sqrt(3)./len_fast_opt;
    case 'matern52'
      lam_fast_opt = sqrt(5)./len_fast_opt;
    otherwise
      error('mapping not implemented')
  end
  %%
  plot_pSTFT_kern_cts(var_fast_opt,omega_opt,lam_fast_opt,kernel1,1,1,4);
