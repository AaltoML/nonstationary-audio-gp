% This demo demonstrates how to use GPPAD to demodulate a
% signal and learn an appropriate time-scale using a grid
% search. In contrast to GPPADDemo2.m, which also learns
% the time-scale, this method can detect whether there
% are multiple local optima. However it is *really*
% slow. This example took 30mins on my laptop.
%
% By placing your sound file in the director called '/Data'
% and replacing the current file name below ('74 -
% Sentences.wav') you will demodulate your sound
% instead. You will probably also need to alter the range
% of the search.

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Paths
addpath ../
LoadLocalPaths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the sentence sound

File = '74 - Sentences.wav'; % Name of file to load
[y,FSamp] = wavread([SoundPath,File]); % reads in the file
RngLim = round([FSamp*1/2+1,2.1*FSamp]);  % Picks a small range
                                   % over which to
                                   % demodulate
DS = 2; % Down sample the original if desired
y = y(RngLim(1):DS:RngLim(2)); 
FSamp = FSamp/DS;

y = y/sqrt(var(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick the time-scale of the demodulation

Opts.LenRng = [8,200];
Opts.NumLens = 25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the demodulation algorithm

[a,c,Params,Info,x,varX,aupper,alower] = GPPAD(y,[],2,Opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Results

PlotSingleChannelErrorBars(y,x,varX,Params,FSamp);

PlotLikelihood(Info{2}.Lens,Info{2}.Objs)

