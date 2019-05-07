
% clear; close all;
addpath('../');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals
addpath(genpath('gppad/'));
addpath(genpath('nmf/'));
addpath(genpath('toolsGP/'));
soundPath = '../../audio/source_sep/training_data/';

  opts.fs = 16000; % sampling rate of file
  opts.D = 16; % number of frequency channels
  opts.N = 3; % number of NMF components
  opts.kernel1 = 'exp'; % kernel for subbands
%   opts.kernel2 = 'matern52'; % kernel for amplitude envelopes
  opts.link = @(g) log(1+exp(g-1)); % link function
  opts.w_lik = 1e-5; % observation noise (testing)
%   opts.inf_type = 'EP'; % 'EP', 'IHGP_EP', 'EKF'
%   opts.p_cubature = 9; % order of cubature for Gaussian integral
%   opts.max_iters = 15; % maximum number of iterations  #######?########
%   likfunc = @likModulatorNMFPowerSq; % models the square of the amplitudes (i.e the spectrogram)
%   opts.likfunc = @likModulatorPreCalcwn; % models the square, and takes pre-calculated sigma pts
%   opts.ep_fraction = 0.75; % Power EP fraction (optimal value varies by sound class)
  %%% training params %%%
%   opts.ep_damping_train = 0.1; % EP damping during training
%   opts.ep_itts_train = 2; % EP iterations during training #######?########
  %%% testing params %%%
%   opts.ep_damping_test = 0.1; % EP damping during testing (0.1 seems best)
%   opts.ep_itts_test = 20; % EP iterations during testing
%   opts.l_iter = 1; % local EKF (inner loop) iterations
%   opts.train_duration = 16000; %  #######?########
  
%   opts.suffix = '_MAT32_MAT52'; % suffix for file - relates to PEP fraction
  
  
%%

filenames = {'011PFNOM_M60_train','011PFNOM_M64_train','011PFNOM_M67_train','131EGLPM_M60_train',...
             '131EGLPM_M64_train','131EGLPM_M67_train','311CLNOM_M60_train','311CLNOM_M64_train',...
             '311CLNOM_M67_train'};
    
%%
for n=1:length(filenames)
    fprintf('training source sep sound %i',n);
%     train_GTFNMF(soundPath,filenames{n},opts);
    train_GTFNMF_fixed_params(soundPath,filenames{n},opts);
end

opts_mat32 = opts;
opts_mat32.kernel1 = 'matern32';

for n=1:length(filenames)
    fprintf('training matern32 kernel on source sep sound %i',n);
    % train_GTFNMF(soundPath,filenames{num},opts);
    train_GTFNMF_fixed_params(soundPath,filenames{n},opts_mat32);
end


%{
%%
disp('training source sep sound 1')
train_GTFNMF(soundPath,'011PFNOM_M60_train',opts);

%%
disp('training source sep sound 2')
train_GTFNMF(soundPath,'011PFNOM_M64_train',opts);

%%
disp('training source sep sound 3')
train_GTFNMF(soundPath,'011PFNOM_M67_train',opts);

%%
disp('training source sep sound 4')
train_GTFNMF(soundPath,'131EGLPM_M60_train',opts);

%%
disp('training source sep sound 5')
train_GTFNMF(soundPath,'131EGLPM_M64_train',opts);

%%
disp('training source sep sound 6')
train_GTFNMF(soundPath,'131EGLPM_M67_train',opts);

%%
disp('training source sep sound 7')
train_GTFNMF(soundPath,'311CLNOM_M60_train',opts);

%%
disp('training source sep sound 8')
train_GTFNMF(soundPath,'311CLNOM_M64_train',opts);

%%
disp('training source sep sound 9')
train_GTFNMF(soundPath,'311CLNOM_M67_train',opts);

%}

