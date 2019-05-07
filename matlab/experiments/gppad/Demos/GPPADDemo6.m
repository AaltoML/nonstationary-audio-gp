% This demo demonstrates how to use GPPAD to do cascade
% demodulation. It takes a spoken sentence sound and then
% recursively demodulates it using GPPAD.m using two user
% specified time-scales. By placing your sound file in the
% director called '/Data' and replacing the current file
% name below ('74 - Sentences.wav') you will recursively
% demodulate your sound instead. You will probably need to
% alter the characteristic time-scales of the
% demodulation cascade (len).
%
% This demo took about a 1.5mins on my machine


clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Paths
addpath ../
LoadLocalPaths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the sentence sound

File = '74 - Sentences.wav'; % Name of file to load
[y,FSamp] = wavread([SoundPath,File]); % reads in the file
RngLim = round([FSamp*1/2+1,5.1*FSamp]);  % Picks a small range
                                   % over which to
                                   % demodulate
DS = 8; % Down sample the original if desired
y = y(RngLim(1):DS:RngLim(2)); 
FSamp = FSamp/DS;

y = y/sqrt(var(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick the time-scale of the demodulation

len = [200/DS;4000/DS]; % Time-scales in samples (in this version of
                      % GPPAD this is the one free parameter which
                      % is not learned from the data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the demodulation algorithm

[A,c] = GPPAD(y,len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Results

PlotCas(A,y,FSamp,[]);

