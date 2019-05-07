% This demo demonstrates how to use GPPAD to do basic
% demodulation. It takes a spoken sentence sound and then
% demodulates it using GPPAD.m using a user specified
% time-scale. By placing your sound file in the director
% called '/Data' and replacing the current file name below
% ('74 - Sentences.wav') you will demodulate your sound
% instead. You also need to set the characteristic
% time-scale of the demodulation.

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Paths
addpath ../
LoadLocalPaths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the sentence sound

File = '74 - Sentences.wav'; % Name of file to load
[y,FSamp] = audioread([SoundPath,File]); % reads in the file
RngLim = round([FSamp*1/2+1,2.1*FSamp]);  % Picks a small range
                                   % over which to
                                   % demodulate
DS = 4; % Down sample the original if desired
y = y(RngLim(1):DS:RngLim(2)); 
FSamp = FSamp/DS;

y = y/sqrt(var(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick the time-scale of the demodulation

len = 200/DS; % Time-scale in samples (in this version of
              % GPPAD this is the one free parameter which
              % is not learned from the data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the demodulation algorithm

[a,c] = GPPAD(y,len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Results

PlotSingleChannel(y,a,FSamp);
