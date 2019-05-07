% This demo demonstrates how to use GPPAD to estimate a
% modulator in a noisy signal that also has missing data.
% It takes a clean spoken sentence sound, adds some noise to
% it, removes the signal from the missing sections, and then
% estimates the modulator using GPPAD. For comparison, the
% modulator is also estimated from the clean signal. The
% modulators estimated on the noisy signal are reasonably
% close to those estimated on the clean signal, and the
% clean modulators lie with in the error-bars of the noisy
% modulators.
%
% By placing your sound file in the director called '/Data'
% and replacing the current file name below ('74 -
% Sentences.wav') you will demodulate your own
% noise-corrupted sound with missing sections instead. You
% will probably also need to alter the initial
% characteristic time-scale of the demodulation, len0.  
%
% This demo takes about a minute to run on my machine.

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
% Making a noisy version of the clean signal
T = length(y);
vary = 0.2*ones(T,1); % Baseline noise variance
varyMissing = 100; % variance of noise in 'missing' sections

MissCent =  [2400,14800,19800]/DS;
MissLength = 200;

yN = y+randn(T,1).*sqrt(vary);

for c=1:length(MissCent)
  ind = [round(MissCent(c)-MissLength/2):1: round(MissCent(c)+MissLength/2)];
  vary(ind) = varyMissing; 
  yN(ind) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters learned from a training sample 
Params.len = 350/DS;
Params.varx = 3.6;
Params.varc = 0.5;
Params.mux =  0.16;
Params.vary = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for noisy sound
ParamsN = Params;
ParamsN.vary = vary;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the demodulation algorithm on the clean and noisy
% signals

[aN,cN,ParamsN,InfoN,xN,varXN,aupperN,alowerN] = GPPAD(yN,ParamsN);
[a,c,Params,Info,x,varX,aupper,alower] = GPPAD(y,Params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Results

PlotSingleChannelErrorBarsCompare(yN,y,xN,varXN,x,varX,ParamsN,Params,FSamp);

