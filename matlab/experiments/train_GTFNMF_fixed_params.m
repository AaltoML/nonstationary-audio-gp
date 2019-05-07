function train_GTFNMF_fixed_params(soundPath,File,opts)

rng(12345,'twister')

fs_ = opts.fs;
D = opts.D;
N = opts.N;
kernel1 = opts.kernel1;

%% Load signal and pre-process
  [y_,fs] = audioread([soundPath,File,'.wav']); % reads in the file
  y = resample(y_, fs_, fs); % downsample the input
  fs = fs_;
  normaliser = sqrt(var(y));
  y_norm = y/normaliser; % rescale the input to unit variance
  T = length(y_norm);
  
  
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
  
  
  
if ~strcmp(kernel1,'exp')
  % Learn properties of the filters (centre frequency and width)
  opts = struct;
  opts.verbose = 1; % view plots of the fitting process
  opts.minT = 500;
  opts.maxT = 1000;
  opts.numIts = 6;
  opts.numLevels = 15;
  opts.bet = 5000;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together
  
  opts.theta_init = [varx;lamx(om_ind);omega];
  [varx,lamx,om] = fit_probSTFT_SD(y_norm,D,kernel1,opts); % trains filters to match the spectrum
  [omega, om_ind] = sort(om);
  if strcmp(kernel1,'matern32')
    lenx = sqrt(3) ./ lamx(om_ind);
  elseif strcmp(kernel1,'matern52')
    lenx = sqrt(5) ./ lamx(om_ind);
  end
  varx = varx(om_ind);
end

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
    logHthresh = log(exp(HEst(:,k)+1e-2)-1);
  
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
  %{
  range_var_fast = [0.01,  0.1]; % subband variance should be fixed
  range_len_fast = [ 100, 2000]; % subband lengthscales should be short
  range_omega    = [   0, 2*pi]; % subband centre frequencies
  range_var_slow = [   1,   10]; % NMF modulator variances
  range_len_slow = [ 200, 5000]; % NMF modulator lengthscales should be long
  range_W        = [   0,  1.5]; % NMF weights
  %}
  W = WEst';
  
  param1 = [0.06.*ones(D,1); lenx; omega];
  param2 = [1.5*varx_nmf; 1.5.*lenx_nmf];
  
  
%% Saving  
  
  results = struct;
  results.File = File;
  results.soundPath = soundPath;
  results.fs = fs;
  results.y_norm = y_norm;
  results.normaliser = normaliser;
  results.D = D; % number of frequency channels
  results.N = N; % number of NMF components
  results.param1 = param1;
  results.param2 = param2;
  results.Wnmf = W;
%   results.inf_type = inf_type; % 'IHGP_EP', 'EP' or 'EKF'
%   results.suffix = suffix;
  matName = strcat('trained_',File,'_',kernel1,'.mat');
  save(matName,'-struct','results');
  
end
  