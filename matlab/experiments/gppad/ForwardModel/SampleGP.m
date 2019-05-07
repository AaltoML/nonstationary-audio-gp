function z = SampleGP(fftCov)

  % function z = SampleGP(fftCov)
  %
  % Samples from a stationary GP  
    
T = length(fftCov);
wn = randn(T,1);
ifftwn = ifft(wn);
z = ifft(sqrt(abs(fftCov)).*fft(wn));

