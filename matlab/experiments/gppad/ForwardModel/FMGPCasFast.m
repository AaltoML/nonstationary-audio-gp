function [y, A, X] = FMGPCasFast(Params,Dims)
  
% function [y, A, X] = FMGPCasFast(Params,Dims)
%  
% Forward model for Cascaded version of GPPAD when the
% covariance function is a squared exponential using the
% FAST version where the duration of the latents is
% 2^(ceil(log2(T+tol*Params.len)))
%
% INPUTS
% Params = Structure of parameters containing
%   len = length scales of the GP [M,1]
%   mux = shift parameters of transformed envelopes [M,1]
%   varc = variance of carrier noise 
%   varx = variances of the amplitude noise [M,1]
%   vary = observations noise (usually set to zero)
% T = number of observations to generate  
%  
% OUTPUTS  
% y = observations [T,1]
% X = transformed envelopes [T,M]
% A = amplitudes [T,M]
% c = carriers  
   
[varx,len,mux,varc,vary] = UnpackParamsGP(Params);

[T,Tx,M] = UnpackDimsCas(Dims); 

fftCov = zeros(Tx,M);
X = zeros(Tx,M);

for m=1:M
  fftCov(:,m) = GetFFTCovFast(len(m),Tx);

  wn = randn(Tx,1);

  ifftwn = ifft(wn);

  X(:,m) = mux(m)+ifft(sqrt(abs(fftCov(:,m))).*fft(wn))*sqrt(varx(m));

end

% Generate the amplitudes
A = log(1+exp(X(1:T,:)));

% Generate the carriers
c = randn(T,1)*sqrt(varc);

% Generate the observations
y = prod(A,2).*c+sqrt(vary)*randn(T,1);

