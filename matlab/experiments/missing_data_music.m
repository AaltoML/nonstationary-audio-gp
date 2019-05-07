clear; close all;
addpath('../');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals

filenames = {'bamboo_flute','cello','clarinet','flute',...
             'guitar','ocarina','piano','piccolo',...
             'sax','toy-accordian'};

%   file_num = 1;
         
%% load trained model

inf_types = {'IHGP_EP'};%{'EP','IHGP_EP','EKF'};

for inf = 1:length(inf_types)
inf_type = inf_types{inf};
for ep_itts = [20]%[1,20]
kernel1 = 'exp';
kernel2 = 'matern52';
% ep_itts = 20;

ep_fraction = 0.75;
ep_damping = 0.1;
if strcmp(inf_type,'IHGP_EP') && ep_itts==20
    ep_damping = 0.0025;
end

w_lik = 1e-4;
p_cubature = 9;
% likfunc = @likModulatorNMFPowerSq;
likfunc = @likModulatorPreCalcwn;
link = @(g) log(1+exp(g-1)); % link function

Esig_store = {};
Vsig_store = {};
RMSE_gaps_store = {};
snr_y_store = {};

for file_num = [1,length(filenames)-1,length(filenames)]%1:length(filenames)
    
    
matName = strcat('trained_',filenames{file_num},'_',kernel1,'.mat');
trained_model = load(matName);




  % number of gap lengths to consider
  numgaps = 1;
  gapLim = [10,320];
%   gapLim = [10,600];

  % set gaps manually to non-silent regions
  switch file_num
      case 1
        gapPos = [1500,5000,7000,9000,13000,18000]; % bamboo_flute
      case 2
        gapPos = [500,1500,3500,5000,8000,10000]; % cello
      case 3
        gapPos = [1000,2500,4000,5500,6000,7500]; % clarinet
      case 4
        gapPos = [800,2200,5000,6500,10000,12500]; % flute
      case 5
        gapPos = [700,1600,2500,6000,8000,11000]; % guitar
      case 6
        gapPos = [700,2000,3000,4000,5000,7000]; % ocarina
      case 7
        gapPos = [800,2000,3000,4000,10000,11000]; % piano
      case 8
        gapPos =[700,2000,5000,6000,9000,10000]; % piccolo
      case 9
        gapPos = [1000,2000,5000,8000,12000,13000]; % sax
      case 10
        gapPos = [1000,3500,7500,8500,10000,15000]; % toy-accordian
  end


D = trained_model.D;
N = trained_model.N;

% link = trained_model.link;
% w_lik_opt = trained_model.w_lik_opt;%trained_model.w_lik_test;
Wnmf = trained_model.Wnmf;
y_norm = trained_model.y_norm;
yTest = y_norm;
T = length(yTest);
t = linspace(1,T,T)';

% pre-calculate the sigma-point weights
[wn,xn_unscaled,~] = utp_ws(p_cubature,N);

% State space model
ss = @(x,p1,p2,kern1,kern2) ss_modulators_nmf(p1,p2,kern1,kern2);

% Moments
%   mom = @(hyp,mu,s2,nmfW,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,p_cubature,ep_fraction,'infEP');
mom = @(hyp,mu,s2,nmfW,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,ep_fraction,wn,xn_unscaled,'infEP');

w_opt = log([w_lik;...
             trained_model.param1(1:D);...
             trained_model.param1(D+1:2*D);...
             trained_model.param1(2*D+1:end);...
             trained_model.param2(1:N);...
             trained_model.param2(N+1:end);...
             Wnmf(:)]);

         
  NGaps = length(gapPos);
  gaps = ceil(linspace(gapLim(1),gapLim(2),numgaps));
  ind = [];
  l = 1; % add a loop here later
  for ng=1:NGaps
    ind = [ind,gapPos(ng)+[-ceil(gaps(l)/2):+ceil(gaps(l)/2)]];
  end
  yTest(ind) = NaN;
  
  
tic
    switch inf_type
      case 'IHGP_EP'
        [Eft,Varft,Covft,lb,ub] = ihgp_ep_modulator_nmf(w_opt,t,yTest,ss,mom,t,kernel1,kernel2,length(w_lik),D,N, ...
                                                             ep_fraction,ep_damping,ep_itts);
      case 'EP'
        [Eft,Varft,Covft,lb,ub] = gf_ep_modulator_nmf(w_opt,t,yTest,ss,mom,t,kernel1,kernel2,length(w_lik),D,N, ...
                                                             ep_fraction,ep_damping,ep_itts);
      case 'EKF'
        [Eft,Varft,Covft,lb,ub] = gf_giekf_modulator_nmf(w_opt,t,yTest,ss,mom,t,kernel1,kernel2,length(w_lik),D,N, ...
                                                             ep_itts,1);
      otherwise
        error('inf_type not recognised');
    end
toc


%% plotting

  Eft_mod = zeros(size(Eft(D+1:end,:)));
  Varft_mod = zeros(size(Varft(D+1:end,:)));
  s = 200;
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
    set(gca,'XTick',[])
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
      set(gca,'XTick',[])
%       title(sprintf('NMF component %d',i))
      if i==1, title('NMF components'), end
    end
  end
  sig_samp = zeros(T,s);
  for v=1:s
%     sig_samp(:,v) = sum((Wnmf * link(mod_samp(:,:,v))) .* sub_samp(:,:,v),1);
    sig_samp(:,v) = sum(sqrt(Wnmf * link(mod_samp(:,:,v))) .* sub_samp(:,:,v),1); % if using sqrt model
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
  
  fprintf('RMSE: %d \n',sqrt(mean((y_norm-Esig).^2)))
  
  
%%
  grey = [0.8 0.8 0.8];
  Esig_ = Esig;
  Esig_(ind) = NaN;
  figure(5); clf
  subplot(2,1,1)
  hold on
  plot(y_norm,'r')
  plot(yTest,'b')
  title('actual signal')
  y_limits = ylim;
  subplot(2,1,2)
  hold on
  plot(Esig,'r')
  plot(Esig_,'b')
  title(sprintf('reconstruction with GTF-NMF model'))
  ylim(y_limits)
  
  
  gaps_to_plot = [1, 2, 3, 4];
  ind1gap = gapPos(gaps_to_plot(1))+[-ceil(gaps(numgaps)/2):+ceil(gaps(numgaps)/2)];
  ind2gap = gapPos(gaps_to_plot(2))+[-ceil(gaps(numgaps)/2):+ceil(gaps(numgaps)/2)];
  ind3gap = gapPos(gaps_to_plot(3))+[-ceil(gaps(numgaps)/2):+ceil(gaps(numgaps)/2)];
  ind4gap = gapPos(gaps_to_plot(4))+[-ceil(gaps(numgaps)/2):+ceil(gaps(numgaps)/2)];
  ind1 = gapPos(gaps_to_plot(1))+[-ceil(2.5*gaps(numgaps)/2):+ceil(2.5*gaps(numgaps)/2)];
  ind2 = gapPos(gaps_to_plot(2))+[-ceil(2.5*gaps(numgaps)/2):+ceil(2.5*gaps(numgaps)/2)];
  ind3 = gapPos(gaps_to_plot(3))+[-ceil(2.5*gaps(numgaps)/2):+ceil(2.5*gaps(numgaps)/2)];
  ind4 = gapPos(gaps_to_plot(4))+[-ceil(2.5*gaps(numgaps)/2):+ceil(2.5*gaps(numgaps)/2)];
  figure(6); clf
  subplot(2,4,1)
  plot(t(ind1),y_norm(ind1),'Color',grey)
  hold on
  plot(t(ind1gap),y_norm(ind1gap), 'k')
  title('actual signal')
  y_limits = ylim;
  subplot(2,4,5)
  plot(t(ind1),Esig(ind1),'Color',grey)
  hold on
  plot(t(ind1gap),Esig(ind1gap), 'b')
  title(sprintf('reconstruction with GTF-NMF model'))
  ylim(y_limits)
  
  subplot(2,4,2)
  plot(t(ind2),y_norm(ind2),'Color',grey)
  hold on
  plot(t(ind2gap),y_norm(ind2gap), 'k')
  title('actual signal')
  y_limits = ylim;
  subplot(2,4,6)
  plot(t(ind2),Esig(ind2),'Color',grey)
  hold on
  plot(t(ind2gap),Esig(ind2gap), 'b')
  title(sprintf('reconstruction with GTF-NMF model'))
  ylim(y_limits)
  subplot(2,4,3)
  plot(t(ind3),y_norm(ind3),'Color',grey)
  hold on
  plot(t(ind3gap),y_norm(ind3gap), 'k')
  title('actual signal')
  y_limits = ylim;
  subplot(2,4,7)
  plot(t(ind3),Esig(ind3),'Color',grey)
  hold on
  plot(t(ind3gap),Esig(ind3gap), 'b')
  title(sprintf('reconstruction with GTF-NMF model'))
  ylim(y_limits)
  subplot(2,4,4)
  plot(t(ind4),y_norm(ind4),'Color',grey)
  hold on
  plot(t(ind4gap),y_norm(ind4gap), 'k')
  title('actual signal')
  y_limits = ylim;
  subplot(2,4,8)
  plot(t(ind4),Esig(ind4),'Color',grey)
  hold on
  plot(t(ind4gap),Esig(ind4gap), 'b')
  title(sprintf('reconstruction with GTF-NMF model'))
  ylim(y_limits)
  
%% SNR
  snr_y = snr(y_norm(ind,:),y_norm(ind,:)-Esig(ind,:));
  RMSE_gaps = sqrt(mean((y_norm(ind,:)-Esig(ind,:)).^2));
  fprintf('Signal SNR: %d \n',snr_y)
  fprintf('Gaps RMSE: %d \n',sqrt(mean((y_norm(ind,:)-Esig(ind,:)).^2)))
%   sound(Esig*normaliser,fs)
  
  Esig_store{file_num} = Esig;
  Vsig_store{file_num} = Vsig;
  RMSE_gaps_store{file_num} = RMSE_gaps;
  snr_y_store{file_num} = snr_y;
  
%% saving
  results = struct;
  results.y_norm = y_norm;
  results.D = D; % number of frequency channels
  results.N = N; % number of NMF components
  results.kernel1 = kernel1; % kernel for subbands
  results.kernel2 = kernel2; % kernel for amplitude envelopes
  results.link = link; % link function
  results.inf_type = inf_type;
  results.ep_itts = ep_itts;
  results.Esig = Esig_store;
  results.Vsig = Vsig_store;
  results.RMSE_gaps = RMSE_gaps_store;
  results.snr_y = snr_y_store;
  results.yTest = yTest;
  matName = strcat('missing_data_music_',inf_type,'_',string(ep_itts),'_',kernel1,'.mat');
  save(matName,'-struct','results');
drawnow;
end

end
end
