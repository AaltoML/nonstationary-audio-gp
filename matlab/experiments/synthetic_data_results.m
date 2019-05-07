clear; close all;
addpath('../');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals

filenames = {'EP_1','EP_20',...
             'IHGP_EP_1','IHGP_EP_20',...
             'EKF_1','EKF_20'};
         
%% load trained model

% Get colors
colorsetup


grey = [0.2,0.2,0.2];
%red = [1 0 0];
%blue = [0 0 1];
%green = [0 1 0];
t = (1:10:5000)/16000*1000;
matName = strcat('synthetic_results_',filenames{1},'.mat');
results = load(matName);
link = results.link;

figure(1);clf
xlabel('Time [ms]')
ylabel('$g_1(t)$')
% title('NMF component 1')
hold on
ylim([0,3.1])
xlim([t(1) t(end)])
box on
set(gca,'layer','top')

% Draw shade
fill([2001 4000 4000 2001 2001]/16,[0 0 3.1 3.1 0],1, ...
    'FaceColor',lcolor(5,:),'EdgeColor',lcolor(5,:))

% Draw gt
plot(t,link(results.y_s(1:10:end,1)),'k','LineWidth',1,'Color',grey)

figure(2);clf
xlabel('Time [ms]')
ylabel('$g_2(t)$')
% title('NMF component 2')
hold on
ylim([0,1.5])
xlim([t(1) t(end)])
box on
set(gca,'layer','top')

% Draw shade
fill([2001 4000 4000 2001 2001]/16,[0 0 3.1 3.1 0],1, ...
    'FaceColor',lcolor(5,:),'EdgeColor',lcolor(5,:))

% Draw gt
h = nan(7,1);
h(1) = plot(t,link(results.y_s(1:10:end,2)),'k','LineWidth',1,'Color',grey)


%colours = {blue,red,green};
c = 1;
for n = 1:length(filenames)
    matName = strcat('synthetic_results_',filenames{n},'.mat');
    results = load(matName);
    if mod(n,2)
        figure(1)
        plot(t,link(results.Emod(1:10:end,1)),'--','Color',color(c,:),'LineWidth',.75)
        figure(2)
        h(n+1) = plot(t,link(results.Emod(1:10:end,2)),'--','Color',color(c,:),'LineWidth',.75)
    else
        figure(1)
        plot(t,link(results.Emod(1:10:end,1)),'Color',color(c,:),'LineWidth',.75)
        figure(2)
        h(n+1) = plot(t,link(results.Emod(1:10:end,2)),'Color',color(c,:),'LineWidth',.75)
        c = c + 1;
    end
end

figure(2)
legend(h,'Ground truth','EP 1 (ADF)','EP 20','IHGP 1','IHGP 20','EKF 1','EKF 20','Location','NorthWest')

figure(1)
% Save figure
  if false
    matlab2tikz('../../paper/figs/synthetic_data_mod1_2.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./figs/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end
figure(2)
% Save figure
  if false
    matlab2tikz('../../paper/figs/synthetic_data_mod2_2.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./figs/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end
  
%%
dind = 4;
t = 2001:5:4000;

figure(3);clf
hold on

% Shade
fill([2001 4000 4000 2001 2001]/16,[-1 -1 1 1 -1]*0.75,1, ...
    'FaceColor',lcolor(5,:),'EdgeColor',lcolor(5,:))

plot(t/16,results.y_f(t,dind),'k-','Color',grey,'LineWidth',0.75)
xlabel('Time [ms]')
ylabel('$z_4(t)$')
c=1;
ylim([-.75 .75])

for n = [2,4,6]%1:length(filenames)
    matName = strcat('synthetic_results_',filenames{n},'.mat');
    results = load(matName);
    if mod(n,2)
        figure(3)
        plot(t/16,results.Esub(t,dind),'Color',color(c,:),'LineWidth',0.75)
%         subplot(2,1,2)
%         plot(t,results.Esub(t,2),'Color',colours{c},'LineWidth',0.5)
    else
        figure(3)
        plot(t/16,results.Esub(t,dind),'-','Color',color(c,:),'LineWidth',0.75)
%         subplot(2,1,2)
%         plot(t,results.Esub(t,2),'--','Color',colours{c},'LineWidth',0.5)
        c = c + 1;
    end
end
xlim([t(1)/16 t(end)/16])
box on
set(gca,'layer','top')

% Save figure
  if false
    matlab2tikz('../../paper/figs/synthetic_data_sub_2.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./figs/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end

  
%%
for n = 1:length(filenames)
    matName = strcat('synthetic_results_',filenames{n},'.mat');
    results = load(matName);
    fprintf('%s RMSE_sig: %g\n',results.inf_type,results.RMSE_sig)
end
for n = 1:length(filenames)
    matName = strcat('synthetic_results_',filenames{n},'.mat');
    results = load(matName);
    fprintf('%s RMSE_sub: %g\n',results.inf_type,results.RMSE_sub)
end
for n = 1:length(filenames)
    matName = strcat('synthetic_results_',filenames{n},'.mat');
    results = load(matName);
    fprintf('%s RMSE_mod: %g\n',results.inf_type,results.RMSE_mod)
end



%%
matName = strcat('synthetic_results_',filenames{2},'.mat');
results = load(matName);
Emod = results.Emod;
y_s = results.y_s;
Wnmf = exp(reshape(results.w_true(end-results.D*results.N+1:end),[results.D,results.N]));

[EU,ES,EV] = svd(Emod);
[yU,yS,yV] = svd(y_s);
[WU,WS,WV] = svd(Wnmf);

figure(4);clf
subplot(2,3,1)
plot(Emod)
subplot(2,3,4)
plot(y_s)
subplot(2,3,2)
plot(EU*ES)
subplot(2,3,5)
plot(yU*yS)
subplot(2,3,3)
plot((WV*Emod')')
subplot(2,3,6)
plot((WV*y_s')')

figure(5);clf
subplot(2,3,1)
plot(link(Emod))
subplot(2,3,4)
plot(link(y_s))
subplot(2,3,2)
plot(link(EU*ES))
subplot(2,3,5)
plot(link(yU*yS))
subplot(2,3,3)
plot((WV*link(Emod'))')
subplot(2,3,6)
plot((WV*link(y_s'))')

