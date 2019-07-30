clear
addpath('symmetric-cubature-rules/'); % code for approximate Gaussian integrals

  D = 12; % number of frequency channels
  N = 3; % number of NMF components
  dt = 1;
  kernel1 = 'matern32'; % kernel for subbands
  kernel2 = 'matern52'; % kernel for amplitude envelopes
  T = 5000;
  t = linspace(1,T,T)';
  mod_sparsity = 0.; % 0 = softplus
  link = @(g) log(1+exp(g-mod_sparsity)); % link function
%   link = @(g) exp(g);
  rng(12345)
  w_lik = 1e-4; % observation noise
  p_cubature = 7; % order of cubature for Gaussian integral
  max_iters = 5; % maximum number of iterations
  
  %%% tune the hyperparameters? %%% (optimisation still quite unstable)
  optimise = 0;
  
  ep_fraction = 0.5; % Power EP fraction
  ep_itts = 3; % EP iterations
  ep_damping = linspace(0.5, 0.5, ep_itts); % EP damping
  
  % define some constraints
  range_var_fast = [0.01,  0.1]; % subband variance should be fixed
  range_len_fast = [  20,  500]; % subband lengthscales should be short
  range_omega    = [   0, 2*pi]; % subband centre frequencies
  range_var_slow = [   2,    5]; % NMF modulator variances
  range_len_slow = [ 200, 2000]; % NMF modulator lengthscales should be long
  range_W        = [   0,  2/D]; % NMF weights
  
  constraints = [range_var_fast; range_len_fast; range_omega; ...
                 range_var_slow; range_len_slow; range_W];
  init_param = @(constraint,D_,N_) constraint(1) + (constraint(end)-constraint(1))*rand(D_,N_);
  
  % which parameters to optimise:
  %              w_lik, var_f, len_f, omega, var_s, len_s, W(nmf)
  tune_hypers = [0,     0,     1,     0,     1,     1,     0     ];
  
  
  % pick some random model parameters
  var_fast = mean(range_var_fast)*ones(D,1);
  len_fast = init_param(constraints(2,:),D,1);
  omega = linspace(pi/4,pi/50,D)';
  var_slow = mean(range_var_slow)*ones(N,1);
  len_slow = linspace(range_len_slow(1)+2,range_len_slow(end)-2,N)';
  W = init_param(constraints(6,:),D,N);
  
  

% State space sampling

  w_sub = [var_fast;len_fast;omega];
  w_mod = [var_slow;len_slow];
  [F,L,Qc,H,Pinf] = ss_modulators_nmf(w_sub,w_mod,kernel1,kernel2);
  [A,Q] = lti_disc(F,L,Qc,dt);

  z = zeros(length(L),T);
  y = zeros(T,1); y_f = zeros(T,D); y_s = zeros(T,N);
%   rng('default');
  z(:,t(1)) = chol(Pinf,'lower')' * randn(length(L),1);
  y(t(1)) = (H(1:D,:) * z(:,t(1)))' * W * H(D+1:D+N,:) * link(z(:,t(1)));
  y_f(t(1),:) = H(1:D,:) * z(:,t(1)); y_s(t(1),:) = H(D+1:D+N,:) * z(:,t(1));
  for k=2:T
    z(:,t(k)) = A * z(:,t(k-1)) + chol(Q)' * randn(length(L),1); % prior
    y_f(t(k),:) = H(1:D,:) * z(:,t(k)); y_s(t(k),:) = H(D+1:D+N,:) * z(:,t(k)); % component ground truth
    y(t(k)) = (H(1:D,:) * z(:,t(k)))' * W * H(D+1:D+N,:) * link(z(:,t(k))); % likelihood
  end
  
  mods = link(y_s) * W';
  figure(1);clf
  for i=1:D
    subplot(D+1,2,2*i-1)
    plot(y_f(:,i))
    hold on
    plot(mods(:,i))
    title(sprintf('channel %d',i))
    if i==1, legend('fast varying signal','slow varying modulator','fontsize',10), end
    if i<=N
      subplot(D+1,2,2*i)
      plot(link(y_s(:,i)),'r')
      title(sprintf('NMF component %d',i))
    end
  end
  subplot(D+1,1,D+1)
  plot(y,'g')
  legend('sum product','fontsize',10)
  title('Observations')
  

%% setup and optimisation

  likfunc = @likModulatorNMFPower;
  
  % Moments
  mom = @(hyp,mu,s2,nmfW,ep_frac,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,p_cubature,ep_frac,'infEP');
  
  % State space model
  ss = @(x,p1,p2,kern1,kern2) ss_modulators_nmf(p1,p2,kern1,kern2);
  
  num_lik_params = length(w_lik);
  
  % Hyperparameters actual
  w_all = {log(w_lik); inv_sigmoid(var_fast,range_var_fast); inv_sigmoid(len_fast,range_len_fast);...
           inv_sigmoid(omega,range_omega);inv_sigmoid(var_slow,range_var_slow); ...
           inv_sigmoid(len_slow,range_len_slow); inv_sigmoid(W(:),range_W)};
       
  w = w_all(tune_hypers>0); w_fixed = w_all(tune_hypers==0);
  w = cell2mat(w);
  w_fixed = cell2mat(w_fixed);

  if optimise
    % Hyperparameters initial
    init_var_fast = mean(range_var_fast)*ones(D,1);
    init_len_fast = init_param(constraints(2,:),D,1);
    init_omega = init_param(constraints(3,:),D,1);
    init_var_slow = init_param(constraints(4,:),N,1);
    init_len_slow = init_param(constraints(5,:),N,1);
    init_W = init_param(constraints(6,:),D,N);
    w_init = {log(w_lik); inv_sigmoid(init_var_fast,range_var_fast); inv_sigmoid(init_len_fast,range_len_fast);...
       inv_sigmoid(init_omega,range_omega); inv_sigmoid(init_var_slow,range_var_slow); ...
       inv_sigmoid(init_len_slow,range_len_slow); inv_sigmoid(init_W(:),range_W)};
    
    w_init = w_init(tune_hypers>0);
    w_init = cell2mat(w_init);
  
    % Optimization options
    opts = optimset('GradObj','off','display','iter','DerivativeCheck','on','MaxIter',max_iters);
  
    % Optimize hyperparameters w.r.t. log marginal likelihood
    [w_opt,ll] = fminunc(@(w_init) gf_ep_modulator_nmf_constraints(w_init,t,y,ss,mom,[],kernel1,kernel2,num_lik_params,D,N,...
                                                                   ep_fraction,ep_damping,ep_itts,constraints,w_fixed,tune_hypers), ...
                           w_init,opts);
  end
  

%% inference and plotting
  
  if exist('w_opt','var'), w_ = w_opt; else, w_ = w; end
  % ADF filtering
  tic
  [Eft,Varft,Covft,lb,ub,out] = gf_ep_modulator_nmf_constraints(w_,t,y,ss,mom,t,kernel1,kernel2,num_lik_params,D,N,...
                                                                ep_fraction,ep_damping,ep_itts,constraints,w_fixed,tune_hypers);
  toc
  
  if tune_hypers(end)
    W_ = reshape(sigmoid(w_(end-D*N+1:end),range_W),[D,N]);
  else
    W_ = reshape(sigmoid(w_fixed(end-D*N+1:end),range_W),[D,N]);
  end
  
  Eft_mod = zeros(size(Eft(D+1:end,:)));
  Varft_mod = zeros(size(Varft(D+1:end,:)));
  s = 1000;
  sub_samp = zeros([size(Eft(1:D,:)),s]);
  mod_samp = zeros([size(Eft(D+1:end,:)),s]);
  
  figure(2);clf
  for i=1:D
    sub_samp(i,:,:) = bsxfun(@plus,bsxfun(@times,randn(T,s),sqrt(Varft(i,:))'),Eft(i,:)');
    subplot(D+1,2,2*i-1)
    plot(y_f(:,i))
    hold on
    plot(Eft(i,:))
    plot(lb(i,:),'g--')
    plot(ub(i,:),'g--')
    title(sprintf('channel %d - subband',i))
    if i==1, legend('ground truth','posterior mean','95% confidence','fontsize',10), end
    if i<=N
      mod_samp(i,:,:) = bsxfun(@plus,bsxfun(@times,randn(T,s),sqrt(Varft(D+i,:))'),Eft(D+i,:)');
      Eft_mod(i,:) = mean(link(mod_samp(i,:,:)),3);
      Varft_mod(i,:) = var(link(mod_samp(i,:,:)),0,3);
      subplot(D+1,2,2*i)
      plot(link(y_s(:,i)))
      hold on
      plot(Eft_mod(i,:))
      plot(Eft_mod(i,:)+1.96*sqrt(Varft_mod(i,:)),'g--')
      plot(Eft_mod(i,:)-1.96*sqrt(Varft_mod(i,:)),'g--')
      title(sprintf('NMF component %d',i))
    end
  end
  sig_samp = zeros(T,s);
  for v=1:s
    sig_samp(:,v) = sum((W_ * link(mod_samp(:,:,v))) .* sub_samp(:,:,v),1);
  end
  Esig = mean(sig_samp,2);
  Vsig = var(sig_samp,0,2);
  subplot(D+1,1,D+1)
  plot(y)
  hold on
  plot(Esig)
  plot(Esig + 1.96*sqrt(Vsig),'g--')
  plot(Esig - 1.96*sqrt(Vsig),'g--')
  title('signal (via NMF)')
  
  figure(3);clf
  plot(y)
  hold on
  plot(Esig)
  plot(Esig + 1.96*sqrt(Vsig),'g--')
  plot(Esig - 1.96*sqrt(Vsig),'g--')
  legend('ground truth','posterior mean','95% confidence','fontsize',10)
  title('signal (via NMF)')
  
  fprintf('RMSE: %d \n',sqrt(mean((y-Esig).^2)))
  fprintf('lZ: %d \n',sum(out.lZ))
  