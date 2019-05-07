% This file produces a three level cascade on the sentence
% sound. By placing your sound file in the director called
% '/Data' and replacing the current file name below ('74 -
% Sentences.wav') you will demodulate you sound instead.
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Paths

addpath ..
LoadLocalPaths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading the sentence sound

File = '74 - Sentences.wav'; % Name of file to load
[y,FSamp] = wavread([SoundPath,File]); % reads in the file
RngLim = [FSamp*1/2+1,6.5*FSamp];

DS = 2; % Down sample the original if desired
y = y(RngLim(1):DS:RngLim(2)); 
FSamp = FSamp/DS;

y = y/sqrt(var(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the options for the demodulation and run it

CASGPMAPOpts; % Load the standard options

% Can modify the above options here:
% e.g. to reduce the number of iterations from 2000 to
% 1000, set Opts.NumIts = 1000;
% see FBGPMAPOpts.m for all of the settings which can be modified
%Opts.NumItsCas = 500;

lenx = [15,250,5000]/DS; % Time-scale in samples (this is the one
                         % free parameter which is not learned from
                         % the data)

% Run the demodulation algorithm
[A,X,c,Params,Info] = CasGPMAP(y,lenx,Opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Results
PlotCas(A,y,FSamp,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you want to save the results

y = single(y);
X = single(X);
c = single(c);
A = single(A);

%save(['YourDirectory',YourSaveName],'y','A', ...
%     'c','X','Params','Opts','FSamp','Info')