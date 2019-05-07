% This demo demonstrates how to use GPPAD to demodulate a
% signal and learn an appropriate time-scale. The method
% also returns error-bars on the modulator. It takes a
% spoken sentence sound and then demodulates it using
% GPPAD.m.  The user must specify an initial guess at the
% time-scale. By placing your sound file in the director
% called '/Data' and replacing the current file name below
% ('74 - Sentences.wav') you will demodulate your sound
% instead. You will probably also need to alter the initial
% characteristic time-scale of the demodulation, len0.
%
% This demo took about 5-10mins on my machine

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
DS = 8; % Down sample the original if desired
y = y(RngLim(1):DS:RngLim(2)); 
FSamp = FSamp/DS;

y = y/sqrt(var(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick the time-scale of the demodulation

len0 = 40; % Time-scale in samples (in this version of GPPAD
           % this is an intial guess which is then altered to
           % fit the data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the demodulation algorithm

[a,c,Params,Info,x,varX,aupper,alower] = GPPAD(y,len0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Results

PlotSingleChannelErrorBars(y,x,varX,Params,FSamp);

