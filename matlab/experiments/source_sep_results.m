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

% ss_results_c = load('../../../inf-hor-mat-files/source_sep_clarinetIHGP_EP_10_exp.mat');
% ss_results_g = load('../../../inf-hor-mat-files/source_sep_guitarIHGP_EP_10_exp.mat');
% ss_results_p = load('../../../inf-hor-mat-files/source_sep_pianoIHGP_EP_10_exp.mat');
ss_results_p = load('source_sep_piano_small.mat');
%%
%{
ss_results_p2 = struct;
ss_results_p2.intrument = ss_results_p.instrument;
ss_results_p2.Esig = ss_results_p.Esig;
ss_results_p2.Vsig = ss_results_p.Vsig;
ss_results_p2.Esig1 = ss_results_p.Esig1;
ss_results_p2.Vsig1 = ss_results_p.Vsig1;
ss_results_p2.Esig2 = ss_results_p.Esig2;
ss_results_p2.Vsig2 = ss_results_p.Vsig2;
ss_results_p2.Esig3 = ss_results_p.Esig3;
ss_results_p2.Vsig3 = ss_results_p.Vsig3;
ss_results_p2.y_normG1Test = ss_results_p.y_normG1Test;
ss_results_p2.y_normG2Test = ss_results_p.y_normG2Test;
ss_results_p2.y_normG3Test = ss_results_p.y_normG3Test;
ss_results_p2.yTest = ss_results_p.yTest;

% matName = strcat('source_sep_piano_small.mat');
save('source_sep_piano_small.mat','-struct','ss_results_p2');
%}

%%%%%%%% THIS IS THE MIXTURE %%%%%%%%

mixture = ss_results_p.yTest;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
colorsetup

grey = [0.2,0.2,0.2];

T = size(ss_results_p.Esig1,1);
t = linspace(1,T/16000,T)';
t_ind = 1:10:T;
figure(1);clf
plot(t(t_ind),ss_results_p.Esig1(t_ind),'Color',color(1,:))
xlabel('Time [secs]')
ylim([-5.5, 4.5])
% Save figure
  if true
    matlab2tikz('../../paper/figs/source_sep1.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./figs/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false) 
  end
  
figure(2);clf
plot(t(t_ind),ss_results_p.Esig2(t_ind),'Color',color(2,:))
ylim([-5.5, 4.5])
xlabel('Time [secs]')
% Save figure
  if true
    matlab2tikz('../../paper/figs/source_sep2.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./figs/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end
  
figure(3);clf
plot(t(t_ind),ss_results_p.Esig3(t_ind),'Color',color(3,:))
ylim([-6.5, 4.5])
xlabel('Time [secs]')
% Save figure
  if true
    matlab2tikz('../../paper/figs/source_sep3.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./figs/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end
  
  
%% Joint figure  
  
colorsetup

grey = [0.2,0.2,0.2];

T = size(ss_results_p.Esig1,1);
t = linspace(1,T/16000,T)';
t_ind = 1:10:T;

a = 8;
b = 6;

figure(1); clf; hold on
  plot(t(t_ind), 1.5*a+min(max(.8*mixture(t_ind),-b-2),b+2),'-','color',grey)
  plot(t(t_ind), 0*a+min(max(ss_results_p.Esig1(t_ind),-b),b),'Color',color(1,:))
  plot(t(t_ind),-1*a+min(max(ss_results_p.Esig2(t_ind),-b),b),'Color',color(2,:))
  plot(t(t_ind),-2*a+min(max(ss_results_p.Esig3(t_ind),-b),b),'Color',color(3,:))
  
  ylim([-23, 16.5])
  set(gcf,'color','white')
  lims = axis;
  
  % Convert to bitmap
  axis off
  A = getframe(gca); cla; A = A.cdata; %A = A.cdata(end:-1:1,:,:); % Flip to correct
  figure(1); clf
  imagesc(linspace(lims(1),lims(2),size(A,2)), ...
          linspace(lims(3),lims(4),size(A,1)),A)  
  axis(lims)
  set(gca,'Ytick',[])
  xlabel('Time [secs]')
  
  labels = {'Source one: piano note C', ...
            'Source two: piano note E', ...
            'Source three: piano note G'};
  
  text(6,-22,sprintf('\\tiny\\bf Input audio, $y$'),'HorizontalAlign','right')  
  for j=0:2
    text(6,-18+a+a*j, ...
        sprintf('\\tiny\\textcolor{mycolor%i}{\\bf %s}',j,labels{j+1}), ...
        'HorizontalAlign','right')  
  end
        
  
% Save figure
  if true
    matlab2tikz('../../paper/figs/source_sep.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./figs/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end

