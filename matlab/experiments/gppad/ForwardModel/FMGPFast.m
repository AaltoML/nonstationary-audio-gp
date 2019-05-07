function [y, a, x] = FMGPFast(Params,Dims)
  
% function [y, x, a] = FMGPFast(Params,Dims)
%  
% Forward model for GPPAD when the covariance function is
% a squared exponential using the FAST version where the
% duration of the latents is 2^(ceil(log2(T+tol*Params.len)))
%
% INPUTS
% Params = Structure of parameters containing
%   len = length scale of the GP
%   mux = shift parameter of transformed envelopes 
%   varc = variance of carrier noise 
%   varx = variance of the amplitude noise
%   vary = observations noise (usually set to zero)
% Dims = structure containing
% T = number of observations to generate  
% Tx = number of latents (includes missing data)  
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

[T,Tx] = UnpackDimsGPFast(Dims); 

fftCov = GetFFTCovFast(len,Tx);

wn = randn(Tx,1);

ifftwn = ifft(wn);

x = mux+real(ifft(sqrt(abs(fftCov)).*fft(wn))*sqrt(varx));

% Generate the amplitudes
a = log(1+exp(x(1:T)));

% Generate the carriers
c = randn(T,1)*sqrt(varc);

% Generate the observations
y = a.*c+sqrt(vary).*randn(T,1);

