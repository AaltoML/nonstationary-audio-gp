function [y, a, x] = FMGP(Params,T)
  
% function [y, x, a] = FMGP(Params,T)
%  
% Forward model for GPPAD when the covariance function is
% a squared exponential
%
% INPUTS
% Params = Structure of parameters containing
%   len = length scale of the GP
%   mux = shift parameter of transformed envelopes 
%   varc = variance of carrier noise 
%   varx = variance of the amplitude noise
%   vary = observations noise (usually set to zero)
% T = number of observations to generate  
%  
% OUTPUTS  
% y = observations [T,1]
% x = transformed envelopes [T,1]
% a = amplitudes [T,1]
% c = carriers  

[varx,len,mux,varc,vary] = UnpackParamsGP(Params);

% Generate the transformed amplitudes
%Col = exp(-1/(2*len^2)*([0:T-1]').^2);
%Col = [Col;Col(T-1:-1:2)];
%fftCol = fft(Col);
fftCov = GetFFTCov(len,T);

wn = randn(2*T-2,1);

ifftwn = ifft(wn);

x = mux+ifft(sqrt(abs(fftCov)).*fft(wn))*sqrt(varx);
x = x(1:T);

% Generate the amplitudes
a = log(1+exp(x));

% Generate the carriers
c = randn(T,1)*sqrt(varc);

% Generate the observations
y = a.*c+sqrt(vary)*randn(T,1);

