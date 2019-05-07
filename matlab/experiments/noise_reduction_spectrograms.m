clear; close all;
addpath('../');
addpath('../../../inf-hor-mat-files/');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals

filenames = {'EP_30_exp'};%,'EP_1_exp','EP_20_exp',...
            % 'IHGP_EP_1_exp','IHGP_EP_20_exp',...
            % 'EKF_1_exp','EKF_20_exp'};
         
%% load trained model

% Get colors
colorsetup


grey = [0.4,0.4,0.4];
%red = [1 0 0];
%blue = [0 0 1];
%green = [0 1 0];
% t = (1:10:5000)/16000*1000;

for n = 1:length(filenames)
    if ~(n == 4)
matName = strcat('noise_reduction_speech_',filenames{n},'.mat');
results = load(matName);

fprintf('%s %i snr_y median: %g\n',results.inf_type,results.ep_itts,median(cell2mat(results.snr_y)))
    end
end

for n = 1:length(filenames)
    if ~(n == 4)
matName = strcat('noise_reduction_speech_',filenames{n},'.mat');
results = load(matName);

fprintf('%s %i RMSE median: %g\n',results.inf_type,results.ep_itts,median(cell2mat(results.RMSE)))
    end
end

matName = strcat('noise_reduction_speech_',filenames{1},'.mat');
results = load(matName);
yClean = results.y_norm;
yNoise = results.yTest;
T = length(results.yTest);
t = linspace(1,T,T)';

EsigEP1 = results.Esig{1};
VsigEP1 = results.Vsig{1};
  
s0 = load('trained_speech0_female.mat');
  param1 = s0.param1;

D = results.D;
[A1,Q1,H1,Pinf1,K1,tau1] = get_disc_model(1./param1(D+1:2*D),param1(1:D),param1(2*D+1:end),D,'exp',6);
  
Esig = results.Esig{1};
  
%% Run GP prob amp demodulation
  ZClean = kernel_ss_probFB(yClean,A1,Q1,H1,Pinf1,K1,0,tau1);
  AClean = abs(ZClean).^0.5;
  ZNoise = kernel_ss_probFB(yNoise,A1,Q1,H1,Pinf1,K1,0,tau1);
  ANoise = abs(ZNoise).^0.5;
  ZE = kernel_ss_probFB(Esig,A1,Q1,H1,Pinf1,K1,0,tau1);
  AE = abs(ZE).^0.5;
  
%
figure(1);clf
imagesc(AClean(:,1:16:end))
set(gca,'YDir','normal')
xlabel('Time[ms]')
ylabel('Frequency channel')
box off


figure(2);clf
imagesc(ANoise(:,1:16:end))
set(gca,'YDir','normal')
xlabel('Time[ms]')
ylabel('Frequency channel')
box off


figure(3);clf
imagesc(AE(:,1:16:end))
set(gca,'YDir','normal')
xlabel('Time[ms]')
ylabel('Frequency channel')
box off

