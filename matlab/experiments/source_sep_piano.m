clear; close all;
addpath('../');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals
addpath(genpath('gppad/'));
addpath(genpath('nmf/'));
addpath(genpath('toolsGP/'));
% soundPath = '../../audio/speech/'; % Specify where to load the data from
% soundPath = '../../audio/textures/';
soundPath = '../../audio/source_sep/training_data/';
soundPathTest = '../../audio/source_sep/test_data/';

instrument = 'piano';

  % piano
  File1 = '011PFNOM_M60_train';
  File2 = '011PFNOM_M64_train';
  File3 = '011PFNOM_M67_train';
  FileM = '011PFNOM_mixture';
  FileG1 = '011PFNOM_C_part';
  FileG2 = '011PFNOM_E_part';
  FileG3 = '011PFNOM_G_part';
  %{
  % electric guitar
  File1 = '131EGLPM_M60_train';
  File2 = '131EGLPM_M64_train';
  File3 = '131EGLPM_M67_train';
  FileM = '131EGLPM_mixture';
  FileG1 = '131EGLPM_C_part';
  FileG2 = '131EGLPM_E_part';
  FileG3 = '131EGLPM_G_part';
  % clarinet
  File1 = '311CLNOM_M60_train';
  File2 = '311CLNOM_M64_train';
  File3 = '311CLNOM_M67_train';
  FileM = '311CLNOM_mixture';
  FileG1 = '311CLNOM_C_part';
  FileG2 = '311CLNOM_E_part';
  FileG3 = '311CLNOM_G_part';
  % ?
  File1 = 'ALVARADO_M60_train';
  File2 = 'ALVARADO_M64_train';
  File3 = 'ALVARADO_M67_train';
  FileM = 'ALVARADO_mixture';
  FileG1 = 'ALVARADO_C_part';
  FileG2 = 'ALVARADO_E_part';
  FileG3 = 'ALVARADO_G_part';
  %}

matName1 = strcat('trained_',File1,'_exp','.mat');
source1_results = load(matName1);
matName2 = strcat('trained_',File2,'_exp','.mat');
source2_results = load(matName2);
matName3 = strcat('trained_',File3,'_exp','.mat');
source3_results = load(matName3);

w_lik = 1e-4;

s1_param1 = source1_results.param1;
s1_param2 = source1_results.param2;
s1_W = source1_results.Wnmf;
s2_param1 = source2_results.param1;
s2_param2 = source2_results.param2;
s2_W = source2_results.Wnmf;
s3_param1 = source3_results.param1;
s3_param2 = source3_results.param2;
s3_W = source3_results.Wnmf;

param1 = {s1_param1,s2_param1,s3_param1};
param2 = {s1_param2,s2_param2,s3_param2};
Wnmf = {s1_W,s2_W,s3_W};

link = @(g) log(1+exp(g-1)); % link function

% likfunc = @likModulatorNMFPowerSq;
likfunc = @likModulatorPreCalcwn;

J = 3; % number of sources
kernel1 = {'exp','exp','exp'};
kernel2 = {'matern52','matern52','matern52'};

inf_type_sep = 'IHGP_EP';

p_cubature = 9;
ep_fraction = 0.75;
ep_damping = 0.025;
ep_itts_test = 10;

% pre-calculate the sigma-point weights
[wn,xn_unscaled,~] = utp_ws(p_cubature,J*source1_results.N);
  
% Moments
% mom = @(hyp,mu,s2,nmfW,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,p_cubature,ep_fraction,'infEP');
mom = @(hyp,mu,s2,nmfW,ep_frac,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,ep_frac,wn,xn_unscaled,'infEP');
  
% State space model
% ss = @(x,p1,p2,kern1,kern2) ss_mixture(p1,p2,kern1,kern2);
ss = @(x,p1,p2,kern1,kern2) ss_modulators_nmf(p1,p2,kern1,kern2);

num_lik_params = length(w_lik);

%%
  fs_ = 16000; % sampling rate of file
  [y_,fs] = audioread([soundPathTest,FileM,'.wav']); % reads in the file
  y = resample(y_, fs_, fs); % downsample the input
  fs = fs_;
  normaliser = sqrt(var(y));
  y_norm = y/normaliser; % rescale the input to unit variance
  
test_ind = [96001:192000];%length(y_norm);

  yTest = y_norm(test_ind);
  T = length(yTest);
  t = linspace(1,T,T)';
  
[y__,fs] = audioread([soundPathTest,FileG1,'.wav']); % reads in the file
  yG1 = resample(y__, fs_, fs); % downsample the input
  fs = fs_;
  y_normG1 = yG1/normaliser; % rescale the input to unit variance
  y_normG1Test = y_normG1(test_ind)./3;
[y__,fs] = audioread([soundPathTest,FileG2,'.wav']); % reads in the file
  yG2 = resample(y__, fs_, fs); % downsample the input
  fs = fs_;
  y_normG2 = yG2/normaliser; % rescale the input to unit variance
  y_normG2Test = y_normG2(test_ind)./3;
[y__,fs] = audioread([soundPathTest,FileG3,'.wav']); % reads in the file
  yG3 = resample(y__, fs_, fs); % downsample the input
  fs = fs_;
  y_normG3 = yG3/normaliser; % rescale the input to unit variance
  y_normG3Test = y_normG3(test_ind)./3;

%%
  w_mix = {log(w_lik), param1, param2, Wnmf};
  tic
    switch inf_type_sep
      case 'IHGP_EP'
        [Eft,Varft,Covft,lb,ub] = ihgp_ep_mods_nmf_mixture(w_mix,t,yTest,ss,mom,t,kernel1,kernel2,J, ...
                                                             ep_fraction,ep_damping,ep_itts_test);
      case 'EP'
        [Eft,Varft,Covft,lb,ub] = gf_ep_mods_nmf_mixture(w_mix,t,yTest,ss,mom,t,kernel1,kernel2,J, ...
                                                             ep_fraction,ep_damping,ep_itts_test);
      otherwise
        error('inf_type not recognised');
    end
  toc



%%

separation_results = struct;
separation_results.instrument = instrument;
separation_results.Eft = Eft;
separation_results.Varft = Varft;
separation_results.lb = lb;
separation_results.ub = ub;
separation_results.yTest = yTest;
separation_results.y_normG1Test = y_normG1Test;
separation_results.y_normG2Test = y_normG2Test;
separation_results.y_normG3Test = y_normG3Test;
matName = strcat('source_sep_',instrument,inf_type_sep,'_',string(ep_itts_test),'_',kernel1{1},'.mat');
save(matName,'-struct','separation_results');


%% plotting
  
  D_ = source1_results.D;
  D = J*source1_results.D;
  N = J*source1_results.N;
  Eft_mod = zeros([N,size(Eft,2)]);
  Varft_mod = zeros([N,size(Varft,2)]);
  s = 100;
  sub_samp = zeros([D,size(Eft,2),s]);
  mod_samp = zeros([N,size(Eft,2),s]);
  
%   figure(2);clf
  for i=1:D
    sub_samp(i,:,:) = bsxfun(@plus,bsxfun(@times,randn(T,s),sqrt(Varft(i,:))'),Eft(i,:)');
    %{
    subplot(ceil(D/2),2,i)
    hold on
    plot(Eft(i,:),'r')
    plot(lb(i,:),'g--')
    plot(ub(i,:),'g--')
    set(gca,'xtick',[])
%     title(sprintf('channel %d - subband',i))
    if i==1, legend('posterior mean','95% confidence','fontsize',10);
        title('Frequency channels'), end
    %}
  end
  
%   figure(3);clf
  for i=1:N
%     subplot(ceil(N/2),2,i)
%     hold on
      mod_samp(i,:,:) = bsxfun(@plus,bsxfun(@times,randn(T,s),sqrt(Varft(D+i,:))'),Eft(D+i,:)');
      Eft_mod(i,:) = mean(link(mod_samp(i,:,:)),3);
      Varft_mod(i,:) = var(link(mod_samp(i,:,:)),0,3);
      %{
      hold on
      plot(Eft_mod(i,:),'r')
      plot(Eft_mod(i,:)+1.96*sqrt(Varft_mod(i,:)),'g--')
      plot(Eft_mod(i,:)-1.96*sqrt(Varft_mod(i,:)),'g--')
      ylim([0, Inf])
%       title(sprintf('NMF component %d',i))
      if i==1, title('NMF components'), end
      %}
  end
  
  W_all = blkdiag(Wnmf{:});
  sig_samp = zeros(T,s);
  sig_samp1 = zeros(T,s);
  sig_samp2 = zeros(T,s);
  sig_samp3 = zeros(T,s);
  envs = zeros(T,D,s);
  channel_samp = zeros(T,D,s);
  for v=1:s
    envs(:,:,v) = sqrt(W_all * link(mod_samp(:,:,v)))';
    channel_samp(:,:,v) = (sqrt(W_all * link(mod_samp(:,:,v))) .* sub_samp(:,:,v))';
    sig_samp(:,v) = sum(sqrt(W_all * link(mod_samp(:,:,v))) .* sub_samp(:,:,v),1);
    sig_samp1(:,v) = sum(channel_samp(:,1:D_,v),2);
    sig_samp2(:,v) = sum(channel_samp(:,D_+1:2*D_,v),2);
    sig_samp3(:,v) = sum(channel_samp(:,2*D_+1:3*D_,v),2);
  end
  envs = mean(envs,3);
  Esig = mean(sig_samp,2);
  Vsig = var(sig_samp,0,2);
  
  Esig1 = mean(sig_samp1,2);
  Vsig1 = var(sig_samp1,0,2);
  Esig2 = mean(sig_samp2,2);
  Vsig2 = var(sig_samp2,0,2);
  Esig3 = mean(sig_samp3,2);
  Vsig3 = var(sig_samp3,0,2);
  
  separation_results.envs = envs;
  separation_results.Esig = Esig;
  separation_results.Vsig = Vsig;
  separation_results.Esig1 = Esig1;
  separation_results.Vsig1 = Vsig1;
  separation_results.Esig2 = Esig2;
  separation_results.Vsig2 = Vsig2;
  separation_results.Esig3 = Esig3;
  separation_results.Vsig3 = Vsig3;
  matName = strcat('source_sep_',instrument,inf_type_sep,'_',string(ep_itts_test),'_',kernel1{1},'.mat');
  save(matName,'-struct','separation_results');
  
  figure(1);clf
  plot(yTest)
  hold on
  plot(real(Esig),'r')
  plot(Esig + 1.96*sqrt(Vsig),'g--')
  plot(Esig - 1.96*sqrt(Vsig),'g--')
  legend('ground truth','posterior mean','95% confidence')
  title('signal (via NMF)')
  
  figure(2);clf
  subplot(3,1,1)
  hold on
  plot(real(Esig1),'r')
  plot(Esig1 + 1.96*sqrt(Vsig1),'g--')
  plot(Esig1 - 1.96*sqrt(Vsig1),'g--')
  legend('posterior mean','95% confidence')
  title('separation - source 1')
  ylim([min(yTest-0.5), max(yTest+0.5)])
  subplot(3,1,2)
  hold on
  plot(real(Esig2),'r')
  plot(Esig2 + 1.96*sqrt(Vsig2),'g--')
  plot(Esig2 - 1.96*sqrt(Vsig2),'g--')
  legend('posterior mean','95% confidence')
  title('separation - source 2')
  ylim([min(yTest-0.5), max(yTest+0.5)])
  subplot(3,1,3)
  hold on
  plot(real(Esig3),'r')
  plot(Esig3 + 1.96*sqrt(Vsig3),'g--')
  plot(Esig3 - 1.96*sqrt(Vsig3),'g--')
  legend('posterior mean','95% confidence')
  title('separation - source 3')
  ylim([min(yTest-0.5), max(yTest+0.5)])
  
  fprintf('RMSE: %d \n',sqrt(mean((yTest-real(Esig)).^2)))  
  
  figure(3);clf
  subplot(3,1,1)
  hold on
  plot(real(Esig1),'r')
  legend('posterior mean')
  title('separation - source 1')
  ylim([min(yTest-0.5), max(yTest+0.5)])
  subplot(3,1,2)
  hold on
  plot(real(Esig2),'r')
  title('separation - source 2')
  ylim([min(yTest-0.5), max(yTest+0.5)])
  subplot(3,1,3)
  hold on
  plot(real(Esig3),'r')
  title('separation - source 3')
  ylim([min(yTest-0.5), max(yTest+0.5)])
                    
                                                            