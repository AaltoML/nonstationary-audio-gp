function x = sampleGPSE(varx,lenx,T)
  
% function x = sampleGPSE(varx,lenx,T)
%  
% Returns a sample from a GP with a squared exponential kernel:
%
% <x(t)x(t')> = varx*exp(-1/(2*lenx^2)*(t-t')^2)
%
% INPUTS
% lenx = length scale of the GP
% varx = variance of the process
% T = duration of the signal to sample
%
% OUTPUTS  
% x = realisation from a GP with a squared exponential kernel, size [T,1]


tau = 5*lenx; % offset added to avoid wrap-around effects 

Tx = 2^ceil(log2(T+tau)); % make the duration a power of 2 so the
                          % fft is faster

fftCov = getGPSESpec(lenx,Tx);

wn = randn(Tx,1);

x = real(ifft(sqrt(fftCov).*fft(wn)))*sqrt(varx);

x = x(1:T);

