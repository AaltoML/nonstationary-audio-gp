
% clear; close all;
addpath('../');
addpath('../prob_filterbank/'); % filterbank code
addpath('../symmetric-cubature-rules/'); % code for approximate Gaussian integrals
addpath(genpath('gppad/'));
addpath(genpath('nmf/'));
addpath(genpath('toolsGP/'));
soundPath = '../../audio/textures/';

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

filenames = {'stim23_bees_buzzing','stim35_boiling_water','stim41_breathing','stim50_camera_snapping_photos',...
             'stim78_chimes_in_the_wind','stim114_dishes_clanking','stim211_keys_jingling','stim224_walking_on_leaves',...
             'stim312_wind','stim398_walking_on_gravel'};
    
%%
for n=1:length(filenames)
    fprintf('training texture sound %i',n);
%     train_GTFNMF(soundPath,filenames{n},opts);
    train_GTFNMF_fixed_params(soundPath,filenames{n},opts);
end

opts_mat32 = opts;
opts_mat32.kernel1 = 'matern32';

for n=1:length(filenames)
    fprintf('training matern32 kernel on texture sound %i',n);
    % train_GTFNMF(soundPath,filenames{num},opts);
    train_GTFNMF_fixed_params(soundPath,filenames{n},opts_mat32);
end

%{
%%
disp('training texture sound 1')
train_GTFNMF(soundPath,'stim23_bees_buzzing',opts);

%%
disp('training texture sound 2')
train_GTFNMF(soundPath,'stim35_boiling_water',opts);

%%
disp('training texture sound 3')
train_GTFNMF(soundPath,'stim41_breathing',opts);

%%
disp('training texture sound 4')
train_GTFNMF(soundPath,'stim50_camera_snapping_photos',opts);

%%
disp('training texture sound 5')
train_GTFNMF(soundPath,'stim78_chimes_in_the_wind',opts);

%%
disp('training texture sound 6')
train_GTFNMF(soundPath,'stim114_dishes_clanking',opts);

%%
disp('training texture sound 7')
train_GTFNMF(soundPath,'stim211_keys_jingling',opts);

%%
disp('training texture sound 8')
train_GTFNMF(soundPath,'stim224_walking_on_leaves',opts);

%%
disp('training texture sound 9')
train_GTFNMF(soundPath,'stim312_wind',opts);

%%
disp('training texture sound 10')
train_GTFNMF(soundPath,'stim398_walking_on_gravel',opts);
             %}


