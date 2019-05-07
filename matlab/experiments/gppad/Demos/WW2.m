% This demo demonstrates how to use GPPAD to do sub-band
% demodulation. It takes a spoken sentence sound, filters it
% through a Gammatone Filter bank, and then uses GPMAP.m
% to demodulate the filter bank. By placing your sound
% file in the director called '/Data' and replacing the
% current file name below ('74 - Sentences.wav') you will
% demodulate your sound instead. 
%
% This example took about 8 minutes on my laptop.

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Paths
addpath ../
LoadLocalPaths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading the sentence sound

File = 'juice.wav'; % Name of file to load
[y,FSamp] = audioread([SoundPath,File]); % reads in the file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering the sound

% Parameters of the filter bank
FLim = [200,3000]; % upper and lower frequency limits on the filters
D = 10; % Number of channels in the filter bank

% Construct the filter using Malcolm Slaney's code
[forward,feedback,CF,ERB,B] = GammaToneMake(FSamp,D, ...
                                            FLim(1),FLim(2),'Moore'); 
DS=1;
% Apply the filter using Malcolm Slaney's code
Y = GammaToneApply(y,forward,feedback)'; %filter into individual channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick the time-scale of the demodulation

len = 200/DS; % Time-scale in samples (in this version of
              % GPPAD this is the one free parameter which
              % is not learned from the data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the demodulation algorithm

[A,C] = GPPAD(Y,len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Results

PlotSpectrogram(A,FSamp,CF,[]);

d = 2;

PlotSingleChannel(Y(:,d),A(:,d),FSamp);
