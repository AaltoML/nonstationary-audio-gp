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


grey = [0.6,0.6,0.6];
%red = [1 0 0];
%blue = [0 0 1];
%green = [0 1 0];
% t = (1:10:5000)/16000*1000;

for n = 1:length(filenames)
matName = strcat('missing_data_music_',filenames{n},'.mat');
results = load(matName);

fprintf('%s %i snr_y mean: %g\n',results.inf_type,results.ep_itts,mean(max(cell2mat(results.snr_y),0)))
% if n == 4
%     keyboard
% end
end

for n = 1:length(filenames)
matName = strcat('missing_data_music_',filenames{n},'.mat');
results = load(matName);

fprintf('%s %i RMSE mean: %g\n',results.inf_type,results.ep_itts,mean(max(cell2mat(results.RMSE_gaps),0)))
end

matName = strcat('missing_data_music_',filenames{2},'.mat');
results = load(matName);
matName = strcat('missing_data_music_',filenames{4},'.mat');
results2 = load(matName);
matName = strcat('missing_data_music_',filenames{6},'.mat');
results3 = load(matName);
T = length(results.yTest);
t = linspace(1,T,T)';

sig_num = 1;
EsigEP1 = results.Esig{sig_num};
VsigEP1 = results.Vsig{sig_num};
y_norm = results.y_norm;
EsigIHGP1 = results2.Esig{sig_num};
VsigIHGP1 = results2.Vsig{sig_num};
EsigEKF1 = results3.Esig{sig_num};
VsigEKF1 = results3.Vsig{sig_num};

numgaps = 1;
  gapLim = [10,320];
  
switch sig_num
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
NGaps = length(gapPos);
  gaps = ceil(linspace(gapLim(1),gapLim(2),numgaps));
  ind = [];
  l = 1; % add a loop here later
  for ng=1:NGaps
    ind = [ind,gapPos(ng)+[-ceil(gaps(l)/2):+ceil(gaps(l)/2)]];
  end
  yTest(ind) = NaN;

signames = {'bamboo_flute','cello','clarinet','flute',...
             'guitar','ocarina','piano','piccolo',...
             'sax','toy-accordian'};
s0 = load(strcat('trained_',signames{sig_num},'_exp.mat'));
  y_norm = s0.y_norm;
  
gaps_to_plot = [3, 2, 3, 4];
  ind1gap = gapPos(gaps_to_plot(1))+[-ceil(gaps(numgaps)/2):+ceil(gaps(numgaps)/2)];
  ind2gap = gapPos(gaps_to_plot(2))+[-ceil(gaps(numgaps)/2):+ceil(gaps(numgaps)/2)];
  ind3gap = gapPos(gaps_to_plot(3))+[-ceil(gaps(numgaps)/2):+ceil(gaps(numgaps)/2)];
  ind4gap = gapPos(gaps_to_plot(4))+[-ceil(gaps(numgaps)/2):+ceil(gaps(numgaps)/2)];
  ind1 = gapPos(gaps_to_plot(1))+[-ceil(2.5*gaps(numgaps)/2):+ceil(2.5*gaps(numgaps)/2)];
  ind2 = gapPos(gaps_to_plot(2))+[-ceil(2.5*gaps(numgaps)/2):+ceil(2.5*gaps(numgaps)/2)];
  ind3 = gapPos(gaps_to_plot(3))+[-ceil(2.5*gaps(numgaps)/2):+ceil(2.5*gaps(numgaps)/2)];
  ind4 = gapPos(gaps_to_plot(4))+[-ceil(2.5*gaps(numgaps)/2):+ceil(2.5*gaps(numgaps)/2)];
%{
  figure(6); clf
  subplot(2,4,1)
  plot(t(ind1),y_norm(ind1),'Color',grey)
  hold on
  plot(t(ind1gap),y_norm(ind1gap), 'k')
  title('actual signal')
  y_limits = ylim;
  subplot(2,4,5)
  plot(t(ind1),EsigEP1(ind1),'Color',grey)
  hold on
  plot(t(ind1gap),EsigEP1(ind1gap), 'b')
  title(sprintf('reconstruction with GTF-NMF model'))
  ylim(y_limits)
  
  subplot(2,4,2)
  plot(t(ind2),y_norm(ind2),'Color',grey)
  hold on
  plot(t(ind2gap),y_norm(ind2gap), 'k')
  title('actual signal')
  y_limits = ylim;
  subplot(2,4,6)
  plot(t(ind2),EsigEP1(ind2),'Color',grey)
  hold on
  plot(t(ind2gap),EsigEP1(ind2gap), 'b')
  title(sprintf('reconstruction with GTF-NMF model'))
  ylim(y_limits)
  subplot(2,4,3)
  plot(t(ind3),y_norm(ind3),'Color',grey)
  hold on
  plot(t(ind3gap),y_norm(ind3gap), 'k')
  title('actual signal')
  y_limits = ylim;
  subplot(2,4,7)
  plot(t(ind3),EsigEP1(ind3),'Color',grey)
  hold on
  plot(t(ind3gap),EsigEP1(ind3gap), 'b')
  title(sprintf('reconstruction with GTF-NMF model'))
  ylim(y_limits)
  subplot(2,4,4)
  plot(t(ind4),y_norm(ind4),'Color',grey)
  hold on
  plot(t(ind4gap),y_norm(ind4gap), 'k')
  title('actual signal')
  y_limits = ylim;
  subplot(2,4,8)
  plot(t(ind4),EsigEP1(ind4),'Color',grey)
  hold on
  plot(t(ind4gap),EsigEP1(ind4gap), 'b')
  title(sprintf('reconstruction with GTF-NMF model'))
  ylim(y_limits)
%}
%%
ind1_ = ind1(61:end-60);
t0=t(ind1_(1));
tt = (t(ind1_)-t0)/16;
figure(1);clf
hold on
% plot(tt,EsigEP1(ind1),'Color',grey)
plot(tt,y_norm(ind1_),'Color',grey,'LineWidth',0.5)
 plot((t(ind1gap)-t0)/16,EsigEP1(ind1gap), 'Color', color(1,:),'LineWidth',0.5)
 plot((t(ind1gap)-t0)/16,EsigIHGP1(ind1gap), 'Color', color(2,:),'LineWidth',0.5)
 plot((t(ind1gap)-t0)/16,EsigEKF1(ind1gap), 'Color', color(3,:),'LineWidth',0.5)
f1=patch([tt',fliplr(tt')],[EsigEP1(ind1_)'-1.96.*sqrt(VsigEP1(ind1_))',fliplr(EsigEP1(ind1_)'+1.96.*sqrt(VsigEP1(ind1_))')],color(1,:));
f1.EdgeAlpha=0;
% ylim([-3.1,3.1])
alpha(f1,0.15);
 xlabel('Time[ms]')
 box on
  
%   title(sprintf('reconstruction with GTF-NMF model'))
%   ylim(y_limits)
xlim([(t(ind1_(1))-t0)/16 (t(ind1_(end))-t0)/16])
set(gca,'layer','top')
legend('y','EP','IHGP','EKF')
box on
