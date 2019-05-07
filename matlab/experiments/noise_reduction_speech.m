clear; close all;
addpath('../');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals

filenames = {'speech0_female','speech1_male','speech2_male','speech3_male',...
             'speech4_male','speech5_male','speech6_female','speech7_female',...
             'speech8_female','speech9_female'};

%   file_num = 1;
         
%% load trained model

% inf_types = {'EP','IHGP_EP','EKF'};
% inf_types = {'IHGP_EP'};
inf_types = {'EP'};

for inf = 1:length(inf_types)
inf_type = inf_types{inf};
% for ep_itts = [1,20]
for ep_itts = [30]

% inf_type = 'IHGP_EP';
kernel1 = 'exp';
kernel2 = 'matern52';
% ep_itts = 1;

ep_fraction = 0.75;
ep_damping = 0.1;
if strcmp(inf_type,'IHGP_EP') && ep_itts==20
    ep_damping = 0.01;
end

noise_levels = [0.01 0.05 0.1 0.3 0.5];
%%%%%%%%%%%%%
% w_lik = 3e-1;
%%%%%%%%%%%%%

p_cubature = 9;
% likfunc = @likModulatorNMFPowerSq;
likfunc = @likModulatorPreCalcwn;
link = @(g) log(1+exp(g-1)); % link function

Esig_store = {};
Vsig_store = {};
RMSE_store = cell(length(noise_levels),length(filenames));
snr_y_store = cell(length(noise_levels),length(filenames));

for file_num = 1%:length(filenames)
    
for j = 4%1:length(noise_levels)
    matName = strcat('trained_',filenames{file_num},'.mat');
    trained_model = load(matName);


    D = trained_model.D;
    N = trained_model.N;

    % link = trained_model.link;
    % w_lik_opt = trained_model.w_lik_opt;%trained_model.w_lik_test;
    Wnmf = trained_model.Wnmf;
    y_norm = trained_model.y_norm;
    yTest = y_norm + sqrt(noise_levels(j))*randn(size(y_norm));
    T = length(yTest);
    t = linspace(1,T,T)';

    % pre-calculate the sigma-point weights
    [wn,xn_unscaled,~] = utp_ws(p_cubature,N);

    % State space model
    ss = @(x,p1,p2,kern1,kern2) ss_modulators_nmf(p1,p2,kern1,kern2);

    % Moments
    %   mom = @(hyp,mu,s2,nmfW,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,p_cubature,ep_fraction,'infEP');
    mom = @(hyp,mu,s2,nmfW,yall,k) feval(likfunc,link,hyp,yall(k),mu,s2,nmfW,ep_fraction,wn,xn_unscaled,'infEP');

    w_opt = log([noise_levels(j);...
                 trained_model.param1(1:D);...
                 trained_model.param1(D+1:2*D);...
                 trained_model.param1(2*D+1:end);...
                 trained_model.param2(1:N);...
                 trained_model.param2(N+1:end);...
                 Wnmf(:)]);




    tic
        switch inf_type
          case 'IHGP_EP'
            [Eft,Varft,Covft,lb,ub] = ihgp_ep_modulator_nmf(w_opt,t,yTest,ss,mom,t,kernel1,kernel2,length(noise_levels(j)),D,N, ...
                                                                 ep_fraction,ep_damping,ep_itts);
          case 'EP'
            [Eft,Varft,Covft,lb,ub] = gf_ep_modulator_nmf(w_opt,t,yTest,ss,mom,t,kernel1,kernel2,length(noise_levels(j)),D,N, ...
                                                                 ep_fraction,ep_damping,ep_itts);
          case 'EKF'
            [Eft,Varft,Covft,lb,ub] = gf_giekf_modulator_nmf(w_opt,t,yTest,ss,mom,t,kernel1,kernel2,length(noise_levels(j)),D,N, ...
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
      figure(5); clf
      subplot(3,1,1)
      hold on
      plot(y_norm,'r')
      title('actual signal')
      subplot(3,1,2)
      plot(yTest,'b')
      title('noisey signal')
      y_limits = ylim;
      subplot(3,1,3)
      hold on
      plot(Esig,'k')
      title(sprintf('reconstruction with GTF-NMF model'))
      ylim(y_limits)

    %% SNR
      snr_y = snr(y_norm,y_norm-Esig);
      RMSE = sqrt(mean((y_norm-Esig).^2));
      fprintf('Signal SNR: %d \n',snr_y)
      fprintf('RMSE: %d \n',sqrt(mean((y_norm-Esig).^2)))
    %   sound(Esig*normaliser,fs)

      Esig_store{file_num} = Esig;
      Vsig_store{file_num} = Vsig;
      RMSE_store{j,file_num} = RMSE;
      snr_y_store{j,file_num} = snr_y;

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
      results.RMSE = RMSE_store;
      results.snr_y = snr_y_store;
      results.yTest = yTest;
      matName = strcat('noise_reduction_speech_',inf_type,'_',string(ep_itts),'_',kernel1,'.mat');
      save(matName,'-struct','results');
    drawnow;
end
end

end
end