clear
addpath('symmetric-cubature-rules/'); % code for approximate Gaussian integrals

  D = 10; % number of frequency channels
  N = 2; % number of NMF components
  dt = 1;
  kernel1 = 'matern32'; % kernel for subbands
  kernel2 = 'matern52'; % kernel for amplitude envelopes
  T = 5000;
  t = linspace(1,T,T)';
  link = @(g) log(1+exp(g)); % link function
%   link = @(g) exp(g);
  rng('default')
  rng(100,'twister')
  w_lik = 1e-4; % observation noise
  p_cubature = 7; % order of cubature for Gaussian integral
  
  %%% tune the hyperparameters? %%%
  optimise = 0;
  
  % pick some random model parameters
  len_fast = 150 + 400*rand(D,1);
  var_fast = 0.01*ones(D,1);
  omega = linspace(pi/3,pi/50,D)';
  len_slow = linspace(200,1500,N)';
  var_slow = 5 + 5*rand(N,1);
  W = 0.1*abs((2.*rand(D,N)).^2-0.2);
  

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
%   likfunc = @likModulatorNMF;
  
  ep_fraction = 1; % Power EP fraction
  
  ep_itts = 10; % EP iterations
  
  ep_damping = linspace(0.1, 0.5, ep_itts);%0.1; % EP damping
  
  % Moments
  mom = @(hyp,mu,s2,nmfW,ep_frac,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,p_cubature,ep_frac,'infEP');
%   mom = @(hyp,mu,s2,nmfW,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,p_cubature,'infEP');
  
  % State space model
  ss = @(x,p1,p2,kern1,kern2) ss_modulators_nmf(p1,p2,kern1,kern2);
  
  num_lik_params = length(w_lik);
  
  % Hyperparameters actual
  w = log([w_lik;var_fast;len_fast;omega;var_slow;len_slow;W(:)]);

  if optimise
    W_init = rand(D,N);
    % Hyperparameters initial
    w_init = log([w_lik;...
             0.01*ones(size(var_fast));8*ones(size(len_fast));omega;...
             ones(size(var_slow));600.*ones(size(len_slow));...
             W_init(:)]);
  
    % Optimization options
    opts = optimset('GradObj','off','display','iter','DerivativeCheck','on','MaxIter',8);
  
    % Optimize hyperparameters w.r.t. log marginal likelihood
    [w_opt,ll] = fminunc(@(w_init) gf_ep_modulator_nmf(w,t,y,ss,mom,[],kernel1,kernel2,num_lik_params,D,N,ep_fraction,ep_damping,ep_itts), ...
                           w_init,opts);
  end
  

%% inference and plotting
  
  if exist('w_opt','var'), w_ = w_opt; else, w_ = w; end
  % EP filtering
  tic
  [Eft,Varft,Covft,lb,ub,out] = gf_ep_modulator_nmf(w_,t,y,ss,mom,t,kernel1,kernel2,num_lik_params,D,N,ep_fraction,ep_damping,ep_itts);
  toc
  
  
  W_ = reshape(exp(w_(num_lik_params+3*D+2*N+1:end)),[D,N]);
  
  Eft_mod = zeros(size(Eft(D+1:end,:)));
  Varft_mod = zeros(size(Varft(D+1:end,:)));
  s = 500;
  sub_samp = zeros([size(Eft(1:D,:)),s]);
  mod_samp = zeros([size(Eft(D+1:end,:)),s]);
  
  mod_mean = W_*link(Eft(D+1:end,:));
  sub_mean = Eft(1:D,:);
  sig_mean = sum(mod_mean .* sub_mean);
  
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
%   plot(sig_mean,'k--')
  legend('ground truth','posterior mean','95% confidence','fontsize',10)
  title('signal (via NMF)')
  
  fprintf('RMSE: %d \n',sqrt(mean((y-Esig).^2)))
