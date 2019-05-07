clear; close all;
addpath('../');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals

filenames = {'EP_1_exp','EP_20_exp',...
             'IHGP_EP_1_exp','IHGP_EP_20_exp',...
             'EKF_1_exp','EKF_20_exp'};
         
%% load trained model

% Get colors
colorsetup

noise_levels = [0.01 0.05 0.1 0.3 0.5];

grey = [0.5,0.5,0.5];
%red = [1 0 0];
%blue = [0 0 1];
%green = [0 1 0];
% t = (1:10:5000)/16000*1000;

snr_y_mean = zeros(length(filenames),length(noise_levels));
snr_y_std = zeros(length(filenames),length(noise_levels));

snr_base = load('noise_reduction_baseline');
snr_base_mean = mean(cell2mat(snr_base.snr_y_store),2)';
snr_base_std = std(cell2mat(snr_base.snr_y_store),[],2)';

for n = 1:length(filenames)
    matName = strcat('noise_reduction_speech_',filenames{n},'.mat');
    results = load(matName);
    snr_y_mean(n,:) = mean(cell2mat(results.snr_y)');
    snr_y_std(n,:) = std(cell2mat(results.snr_y)');

% fprintf('%s %i snr_y median: %g\n',results.inf_type,results.ep_itts,median(cell2mat(results.snr_y)))
%     end
end

figure(1);clf
hold on
plot(noise_levels,snr_y_mean(1,:),'--','Color',color(1,:),'LineWidth',1.)
plot(noise_levels,snr_y_mean(2,:),'Color',color(1,:),'LineWidth',1.)
plot(noise_levels,snr_y_mean(3,:),'--','Color',color(2,:),'LineWidth',1.)
plot(noise_levels,snr_y_mean(4,:),'Color',color(2,:),'LineWidth',1.)
plot(noise_levels,snr_y_mean(5,:),'--','Color',color(3,:),'LineWidth',1.)
plot(noise_levels,snr_y_mean(6,:),'Color',color(3,:),'LineWidth',1.)
plot(noise_levels,snr_base_mean,'Color',grey,'LineWidth',1.)
f1=patch([noise_levels,fliplr(noise_levels)],[snr_y_mean(1,:)-snr_y_std(1,:),fliplr(snr_y_mean(1,:)+snr_y_std(1,:))],color(1,:));
f2=patch([noise_levels,fliplr(noise_levels)],[snr_y_mean(2,:)-snr_y_std(2,:),fliplr(snr_y_mean(2,:)+snr_y_std(2,:))],color(1,:));
f3=patch([noise_levels,fliplr(noise_levels)],[snr_y_mean(3,:)-snr_y_std(3,:),fliplr(snr_y_mean(3,:)+snr_y_std(3,:))],color(2,:));
f4=patch([noise_levels,fliplr(noise_levels)],[snr_y_mean(4,:)-snr_y_std(4,:),fliplr(snr_y_mean(4,:)+snr_y_std(4,:))],color(2,:));
f5=patch([noise_levels,fliplr(noise_levels)],[snr_y_mean(5,:)-snr_y_std(5,:),fliplr(snr_y_mean(5,:)+snr_y_std(5,:))],color(3,:));
f6=patch([noise_levels,fliplr(noise_levels)],[snr_y_mean(6,:)-snr_y_std(6,:),fliplr(snr_y_mean(6,:)+snr_y_std(6,:))],color(3,:));
f7=patch([noise_levels,fliplr(noise_levels)],[snr_base_mean-snr_base_std,fliplr(snr_base_mean+snr_base_std)],grey);
f1.EdgeAlpha=0;f2.EdgeAlpha=0;f3.EdgeAlpha=0;f4.EdgeAlpha=0;f5.EdgeAlpha=0;f6.EdgeAlpha=0;f7.EdgeAlpha=0;
alpha(f1,0.15);alpha(f2,0.15);alpha(f3,0.15);alpha(f4,0.15);alpha(f5,0.15);alpha(f6,0.15);alpha(f7,0.15);
legend('EP 1','EP 20','IHGP 1','IHGP 20','EKF 1','EKF 20','SpecSub')
set(gca(),'XTick',noise_levels([1, 3:5]))
xlim([noise_levels(1),noise_levels(end)])
ylim([0, 15])
xlabel('Corrupting noise variance')
ylabel('SNR [dB]')
box on

% Save figure
  if true
    matlab2tikz('../../paper/figs/noise_reduction_comparison.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./figs/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end
