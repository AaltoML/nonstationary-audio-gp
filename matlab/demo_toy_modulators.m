clear
addpath('symmetric-cubature-rules/'); % code for approximate Gaussian integrals

  var_fast = [ 0.1;   0.1];%   0.4];
  len_fast = [ 50.;   40.];%    20];
  omega    = [pi/4;  pi/6];%  pi/8];
  var_slow = [  2.;    3.];%    4.];
  len_slow = [500.;  700.];%    90];
  D = length(omega);
  dt = 1;
  kernel1 = 'matern32';
  kernel2 = 'matern52';
  T = 2000;
  t = linspace(1,T,T)';
  link = @(g) log(1+exp(g)); % link function
%   link = @(g) exp(g);
  rng('default')
  rng(123)
  w_lik = 1e-8; % observation noise
  p_cubature = 9; % order of cubature for Gaussian integral
  
  %%% tune the hyperparameters? %%%
  optimise = 0;
  

% State space sampling

  w = [var_fast;len_fast;omega;var_slow;len_slow];
  [F,L,Qc,H,Pinf] = ss_modulators(w,kernel1,kernel2);
  [A,Q] = lti_disc(F,L,Qc,dt);

  z = zeros(length(L),T);
  y = zeros(T,1); y_f = zeros(T,D); y_s = zeros(T,D);
  z(:,t(1)) = chol(Pinf,'lower')' * randn(length(L),1);
  y(t(1)) = (H(1:D,:) * z(:,t(1)))' * H(D+1:2*D,:) * link(z(:,t(1)));
  y_f(t(1),:) = H(1:D,:) * z(:,t(1)); y_s(t(1),:) = H(D+1:2*D,:) * z(:,t(1));
  for k=2:T
    z(:,t(k)) = A * z(:,t(k-1)) + chol(Q)' * randn(length(L),1); % prior
    y_f(t(k),:) = H(1:D,:) * z(:,t(k)); y_s(t(k),:) = H(D+1:2*D,:) * z(:,t(k)); % component ground truth
    y(t(k)) = (H(1:D,:) * z(:,t(k)))' * H(D+1:2*D,:) * link(z(:,t(k))); % likelihood
  end
  
  figure(1);clf
  subplot(4,1,1)
  plot(y_f(:,1))
  hold on
  plot(link(y_s(:,1)))
  legend('fast varying signal 1','slow varying modulator 1','fontsize',10)
  title('Channel 1')
  subplot(4,1,2)
  plot(y_f(:,2))
  hold on
  plot(link(y_s(:,2)))
  legend('fast varying signal 2','slow varying modulator 2','fontsize',10)
  title('Channel 2')
  %{
  subplot(4,1,3)
  plot(y_f(:,3))
  hold on
  plot(link(y_s(:,3)))
  legend('fast varying signal 3','slow varying modulator 3','fontsize',10)
  title('Channel 3')
  %}
  subplot(4,1,4)
  plot(y,'g')
  legend('sum product','fontsize',10)
  title('Observations')
  

%% setup and optimisation

  likfunc = @likModulator;
  
  % Moments
  mom = @(hyp,mu,s2,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,p_cubature,'infEP');
  
  % State space model
  ss = @(x,p,kern1,kern2) ss_modulators(p,kern1,kern2);
  
  num_lik_params = length(w_lik);
  
  % Hyperparameters actual
  w = log([w_lik;var_fast;len_fast;omega;var_slow;len_slow]);

  if optimise
    % Hyperparameters initial
    w = log([w_lik;ones(size(omega));50*ones(size(omega));ones(size(omega));omega;5.*ones(size(omega))]);
  
    % Optimization options
    opts = optimset('GradObj','off','display','iter','DerivativeCheck','on','MaxIter',10);
  
    % Optimize hyperparameters w.r.t. log marginal likelihood
    [w2,ll] = fminunc(@(w) gf_adf_modulator(w,t,y,ss,mom,[],kernel1,kernel2,num_lik_params), ...
                           w,opts);
  end
  

%% inference and plotting
  
  if exist('w2','var'), w_ = w2; else, w_ = w; end
  % ADF filtering
  tic
  [Eft,Varft,Covft,lb,ub,out] = gf_ep_modulator(w_,t,y,ss,mom,t,kernel1,kernel2,num_lik_params);
  toc
  
  % GHKF filtering (also feat EKF filter, but commented out)
  tic
  %[Eft,Varft,Covft,lb,ub,out] = ckf_modulator(w_,t,y,ss,mom,t,kernel1,kernel2,num_lik_params,link);
  toc
  
  
  
  Eft_mod = zeros(size(Eft(D+1:end,:)));
  Varft_mod = zeros(size(Varft(D+1:end,:)));
  s = 1000;
  sub_samp = zeros([size(Eft(1:D,:)),s]);
  mod_samp = zeros([size(Eft(D+1:end,:)),s]);
  
  figure(2);clf
  for i=1:D
    sub_samp(i,:,:) = bsxfun(@plus,bsxfun(@times,randn(T,s),sqrt(Varft(i,:))'),Eft(i,:)');
    subplot(2*D+1,1,i)
    plot(y_f(:,i))
    hold on
    plot(Eft(i,:))
    plot(lb(i,:),'g--')
    plot(ub(i,:),'g--')
    title(sprintf('channel %d - subband',i))
    if i==1, legend('ground truth','posterior mean','95% confidence','fontsize',10), end
    mod_samp(i,:,:) = bsxfun(@plus,bsxfun(@times,randn(T,s),sqrt(Varft(D+i,:))'),Eft(D+i,:)');
    Eft_mod(i,:) = mean(link(mod_samp(i,:,:)),3);
    Varft_mod(i,:) = var(link(mod_samp(i,:,:)),0,3);
    subplot(2*D+1,1,D+i)
    plot(link(y_s(:,i)))
    hold on
    plot(Eft_mod(i,:))
    plot(Eft_mod(i,:)+1.96*sqrt(Varft_mod(i,:)),'g--')
    plot(Eft_mod(i,:)-1.96*sqrt(Varft_mod(i,:)),'g--')
    title(sprintf('channel %d - modulator',i))
  end
  sig_samp = zeros(T,s);
  for v=1:s
    sig_samp(:,v) = sum(link(mod_samp(:,:,v)) .* sub_samp(:,:,v),1);
  end
  Esig = mean(sig_samp,2);
  Vsig = var(sig_samp,0,2);
  subplot(2*D+1,1,2*D+1)
  plot(y)
  hold on
  plot(Esig)
  plot(Esig + 1.96*sqrt(Vsig),'g--')
  plot(Esig - 1.96*sqrt(Vsig),'g--')
  title('sum product')  
  %%
  fprintf('RMSE: %d \n',sqrt(mean((y-Esig).^2)))
  