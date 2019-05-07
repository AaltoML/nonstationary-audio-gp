function fftCov = GetFFTCovFast(len,T);
  
  
  % Col = exp(-1/(2*len^2)*([0:T-1]').^2);
  % Col(1) = Col(1)+1e-6;
%   Col = [Col;Col(T-1:-1:2)];
  
%   fftCov = fft(Col);
  
  tiny = 1e-6;
    
  omega = 2*pi*[[0:floor(T/2)],[-floor(T/2)+1:1:-1]]'/T;

  fftCov = sqrt(2*pi*len^2)*exp(-omega.^2*len.^2/2);

  fftCov = fftCov+tiny;
  