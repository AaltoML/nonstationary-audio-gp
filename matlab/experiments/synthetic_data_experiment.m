function synthetic_data_experiment(inf_type,ep_itts,ep_damping,save_results)

if nargin < 4
    save_results = 0;
end
% clear; close all;
addpath('../');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals
soundPath = '../../audio/speech/'; % Specify where to load the data from
% soundPath = '../audio/textures/';
% soundPath = '../audio/music/';

  % load signal
  File = 'speech0_female'; % Name of file to load
  fs_ = 16000; % sampling rate of file
  
  D = 5; % number of frequency channels
  N = 2; % number of NMF components
  dt = 1;
  kernel1 = 'matern32'; % kernel for subbands
  kernel2 = 'matern52'; % kernel for amplitude envelopes
  T = 5000;
  t = linspace(1,T,T)';
%   mod_sparsity = 1.; % 0 = softplus
  link = @(g) log(1+exp(g-1)); % link function
  rng(12345,'twister')
  w_lik = 1e-4;
  p_cubature = 9; % order of cubature for Gaussian integral

%   likfunc = @likModulatorNMFPower;
  likfunc = @likModulatorNMFPowerSq;
  
  ep_fraction = 0.75; % Power EP fraction (0.8-1 seems best)
  
  %%% testing params %%%
  if nargin < 3
    ep_damping = 0.1; % EP damping during testing (0.1 seems best)
  end
%   ep_itts = 1; % EP iterations during testing
  
  l_iter = 1; % local EKF (inner loop) iterations
  
%   inf_type = 'EKF'; % 'EP', 'EKF', 'IHGP_EP'
  
  
%% Load signal and pre-process
  [y,fs] = audioread([soundPath,File,'.wav']); % reads in the file
  yTest = resample(y, fs_, fs); % downsample the input
  fs = fs_;
  normaliser = sqrt(var(yTest));
  yTest = yTest/normaliser; % rescale the input to unit variance
  yTest = yTest(4501:4501+T-1);
  
  
%% Pre-optimise the filter bank with exponential kernel (probabilistic phase vocoder)
  yTrain = yTest;
  % Learn properties of the filters (centre frequency and width)
  opts = struct;
  opts.verbose = 0; % view plots of the fitting process
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
  
  
%%
  ZTest1 = kernel_ss_probFB(yTest,A1,Q1,H1,Pinf1,K1,0,tau1);
  ATest1 = abs(ZTest1');
  [Hnmf,Wnmf] = nnmf(ATest1,N);
%   Wnmf = Wnmf' / (D/10 * max(Wnmf(:)));
  Wnmf = Wnmf' / max(Wnmf(:));
  Wnmf = max(Wnmf,0.1);
  
  
%% Now set up the full model  
  % define some constraints
  range_var_fast = [0.01,  0.1]; % subband variance should be fixed
  range_len_fast = [ 100, 1000]; % subband lengthscales should be short
  range_omega    = [   0, 2*pi]; % subband centre frequencies
  range_var_slow = [   1,   10]; % NMF modulator variances
  range_len_slow = [ 300, 5000]; % NMF modulator lengthscales should be long
  range_W        = [   0, 1.25]; % NMF weights
  
  
  % pick some random model parameters
  var_fast = mean(range_var_fast)*ones(D,1);
  len_fast = min(max(lenx,range_len_fast(1)),range_len_fast(end));
%   var_fast = varx;
%   omega = linspace(pi/4,pi/50,D)';
  var_slow = mean(range_var_slow)*ones(N,1);
  len_slow = linspace(range_len_slow(1)+2,range_len_slow(end)/10-2,N)';
  W = Wnmf; %init_param(constraints(6,:),D,N);
  
  
%% Prior predictive checks
  % State space sampling
  w_sub = [var_fast;len_fast;omega];
  w_mod = [var_slow;len_slow];
  
  [F,L,Qc,H,Pinf] = ss_modulators_nmf(w_sub,w_mod,kernel1,kernel2);
  [A,Q] = lti_disc(F,L,Qc,dt);
  
  z = zeros(length(L),T);
  y_samp = zeros(T,1); y_f = zeros(T,D); y_s = zeros(T,N);
  z(:,t(1)) = chol(Pinf,'lower')' * randn(length(L),1);
  y_samp(t(1)) = (H(1:D,:) * z(:,t(1)))' * sqrt(W * H(D+1:D+N,:) * link(z(:,t(1))));
  y_f(t(1),:) = H(1:D,:) * z(:,t(1)); y_s(t(1),:) = H(D+1:D+N,:) * z(:,t(1));
  for k=2:T
    z(:,t(k)) = A * z(:,t(k-1)) + chol(Q)' * randn(length(L),1); % prior
    y_f(t(k),:) = H(1:D,:) * z(:,t(k)); y_s(t(k),:) = H(D+1:D+N,:) * z(:,t(k)); % component ground truth
    y_samp(t(k)) = (H(1:D,:) * z(:,t(k)))' * sqrt(W * H(D+1:D+N,:) * link(z(:,t(k)))); % likelihood
  end
  
  mods = link(y_s) * W';
  figure(1);clf
  for i=1:D
    subplot(D+1,2,2*i-1)
    plot(y_f(:,i))
    hold on
    plot(mods(:,i))
    if i==1
      legend('fast varying signal','slow varying modulator','fontsize',10)
      title('Prior predictive checks')
    end
    if i<=N
      subplot(D+1,2,2*i)
      plot(link(y_s(:,i)),'r')
    end
  end
  subplot(D+1,1,D+1)
  plot(y_samp,'g')
  legend('sum product','fontsize',10)
  
  
%% setup and optimisation
  
  % Moments
  mom = @(hyp,mu,s2,nmfW,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,p_cubature,ep_fraction,'infEP');
  
  % State space model
  ss = @(x,p1,p2,kern1,kern2) ss_modulators_nmf(p1,p2,kern1,kern2);
  
  num_lik_params = length(w_lik);
  
  w_true = log([w_lik;...
                var_fast;...
                len_fast;...
                omega;...
                var_slow;...
                len_slow;...
                W(:)]);             
  
  
%% inference and plotting
  
  % EP filtering
  tic
  switch inf_type
      case 'EP'
        [Eft,Varft,~,lb,ub] = gf_ep_modulator_nmf(w_true,t,y_samp,ss,mom,t,kernel1,kernel2,num_lik_params,D,N, ...
                                                                ep_fraction,ep_damping,ep_itts);
      case 'IHGP_EP'
        [Eft,Varft,~,lb,ub] = ihgp_ep_modulator_nmf(w_true,t,y_samp,ss,mom,t,kernel1,kernel2,num_lik_params,D,N, ...
                                                                ep_fraction,ep_damping,ep_itts);
      case 'EKF'
        [Eft,Varft,~,lb,ub] = gf_giekf_modulator_nmf(w_true,t,y_samp,ss,mom,t,kernel1,kernel2,num_lik_params,D,N, ...
                                                                ep_itts,l_iter);
  end
  toc
  
  
 
  

 
 %% plotting
 
  Eft_mod = zeros(size(Eft(D+1:end,:)));
  Varft_mod = zeros(size(Varft(D+1:end,:)));
  s = 10000;
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
%     title(sprintf('channel %d - subband',i))
    if i==1, legend('posterior mean','95% confidence','fontsize',10);
        title('Frequency channels'), end
    if i<=N
      mod_samp(i,:,:) = bsxfun(@plus,bsxfun(@times,randn(T,s),sqrt(Varft(D+i,:))'),Eft(D+i,:)');
      Eft_mod(i,:) = mean(link(mod_samp(i,:,:)),3);
      Varft_mod(i,:) = var(link(mod_samp(i,:,:)),0,3);
      subplot(D+1,2,2*i)
      hold on
      plot(Eft_mod(i,:),'r')
      plot(Eft_mod(i,:)+1.96*sqrt(Varft_mod(i,:)),'g--')
      plot(Eft_mod(i,:)-1.96*sqrt(Varft_mod(i,:)),'g--')
%       title(sprintf('NMF component %d',i))
      if i==1, title('NMF components'), end
    end
  end
  sig_samp = zeros(T,s);
  for v=1:s
%     sig_samp(:,v) = sum((W * link(mod_samp(:,:,v))) .* sub_samp(:,:,v),1);
    sig_samp(:,v) = sum(sqrt(W * link(mod_samp(:,:,v))) .* sub_samp(:,:,v),1); % if using sqrt model
  end
  Esig = mean(sig_samp,2);
  Vsig = var(sig_samp,0,2);
  subplot(D+1,1,D+1)
  plot(y_samp)
  hold on
  plot(Esig,'r')
  plot(Esig + 1.96*sqrt(Vsig),'g--')
  plot(Esig - 1.96*sqrt(Vsig),'g--')
  title('signal (via NMF)')
  
  figure(3);clf
  plot(y_samp)
  hold on
  plot(Esig,'r')
  plot(Esig + 1.96*sqrt(Vsig),'g--')
  plot(Esig - 1.96*sqrt(Vsig),'g--')
  legend('ground truth','posterior mean','95% confidence','fontsize',10)
  title('signal (via NMF)')
  
  RMSE_sig = sqrt(mean((y_samp-Esig).^2));
  fprintf('signal RMSE: %d \n',RMSE_sig)
  
  Esub = mean(sub_samp,3)';
  Vsub = var(sub_samp,[],3)';
  RMSE_sub = sqrt(mean((y_f-Esub).^2,[1,2]));
  fprintf('subband RMSE: %d \n',RMSE_sub)
  
  Emod = mean(mod_samp,3)';
  Vmod = var(mod_samp,[],3)';
  RMSE_mod = sqrt(mean((link(y_s)-link(Emod)).^2,[1,2]));
  fprintf('modulator RMSE: %d \n',RMSE_mod)
  
%% saving
if save_results
  results = struct;
  results.y_samp = y_samp;
  results.D = D; % number of frequency channels
  results.N = N; % number of NMF components
  results.kernel1 = kernel1; % kernel for subbands
  results.kernel2 = kernel2; % kernel for amplitude envelopes
  results.link = link; % link function
  results.w_true = w_true;
  results.inf_type = inf_type;
  results.ep_itts = ep_itts;
  results.Eft = Eft;
  results.Varft = Varft;
  results.lb = lb;
  results.ub = ub;
  results.Esig = Esig;
  results.Vsig = Vsig;
  results.Esub = Esub;
  results.Vsub = Vsub;
  results.Emod = Emod;
  results.Vmod = Vmod;
  results.y_f = y_f;
  results.y_s = y_s;
  results.RMSE_sig = RMSE_sig;
  results.RMSE_sub = RMSE_sub;
  results.RMSE_mod = RMSE_mod;
  matName = strcat('synthetic_results_',inf_type,'_',string(ep_itts),'.mat');
  save(matName,'-struct','results');
end  
 
end