function [fftCov,varargout] = getGPSESpec(len,T);
  

% function [fftCov,varargout] = getGPSESpec(len,T);
%
% Returns the spectrum of a GP with a squared exponential kernel.
%
% <x(t)x(t')> = exp(-1/(2*len^2)*(t-t')^2)
%
% INPUTS
% len = length scale of the GP, [D,1]
% T = number of frequencies over which to evaluate the spectrum
%
% OUTPUTS
% fftCov = power spectrum of the GP, size [T,D]
% dfftCov = derivative wrt length-scales [T,D]
   
  omega = 2*pi*[[0:ceil(T/2)],[-floor(T/2)+1:1:-1]]'/T;
  
  D = length(len);
  fftCov = zeros(T,D);

  for d=1:D
    fftCov(:,d) = sqrt(2*pi*len(d)^2)*exp(-omega.^2*len(d).^2/2);
  end

  if nargout>1
      dfftCov = zeros(T,D);

    for d=1:D
      dlogfftCovdlen = 1./len-omega.^2*len;
      dfftCov(:,d) = fftCov(:,d).*dlogfftCovdlen;
    end
    varargout{1} = dfftCov;
  end
  
  %  tiny = 1e-6;
%  fftCov = fftCov+tiny;
  