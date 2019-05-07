function freqSq = getGPSEVarSpec(len);
  
% function freqSq = getGPSEVarSpec(len);
%
% Returns the variance of the spectrum of a GP with a squared
% exponential kernel.
%
% powSpec =sqrt(2*pi*len^2)*exp(-omega.^2*len.^2/2);
%
% INPUTS
% len = length scale of the GP
%
% OUTPUTS
% freqSq = average square frequency

  freqSq = 1/(2*pi*len)^2;

%    omega = 2*pi*[[0:floor(T/2)],[-floor(T/2)+1:1:-1]]'/T;
%  fftCov =sqrt(2*pi*len^2)*exp(-omega.^2*len.^2/2);


  