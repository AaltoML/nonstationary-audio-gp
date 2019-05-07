function [obj,dobj] = getObjSEGP(param,specy,varargin)

% function [obj,dob] = getObjSEGP(param,specy)
%
% negative likelihood of a squared exponential GP evaluated on
% regularly sampled data using the spectral representation of the GP
%
% specy = abs(fft(y)).^2
% <y(t)y(t')> = varx*exp(-1/(2*len^2)*(t-t')^2)+vary
%
% evaluates negative likelihood and the derivatives wrt log(len),log(varx),log(vary)
%
% function [obj,dob] = getObjSEGP(param,specy,vary)
%
% user can specify the obervationtion noise, param = 2*1 vector
%
%
% INPUTS 
% param = either a 
%   length-2 vector containing log(len) and log(varx), or
%   length-3 vector containing the above and log(vary)
% specy = spectrum of the data (abs(fft(y)).^2) [T,1]
% optional input if length(param)==2
% vary = observation noise
%
% OUPUTS
% obj = negative log-likelihood
% dobj = derivatives of the log-likelihood wrt log(len) log(varx)
%        and optionally log(vary) [2,1] or [3,1]


T = length(specy);

if nargin==2
  lenx = exp(param(1));
  varx = exp(param(2));
  vary = exp(param(3));
  
  
  [fftCov,dfftCov_dlen] = getGPSESpec(lenx,T);

  obj = 1/2*sum(log(varx*fftCov+vary))+1/(2*T)*sum(specy./(varx*fftCov+vary));

  dobdjdspec = 1./(2*(varx*fftCov+vary))-1/(2*T)*specy./((varx*fftCov+vary).^2);

  dobj = zeros(3,1);
  dobj(1) = sum(dobdjdspec.*dfftCov_dlen)*lenx*varx;
  dobj(2) = sum(dobdjdspec.*fftCov)*varx;
  dobj(3) = sum(dobdjdspec)*vary;

elseif nargin==3
  
  lenx = exp(param(1));
  varx = exp(param(2));
  vary = varargin{1};
  
  [fftCov,dfftCov_dlen] = getGPSESpec(lenx,T);
  
  obj = 1/2*sum(log(varx*fftCov+vary))+1/(2*T)*sum(specy./(varx*fftCov+vary));

  dobdjdspec = 1./(2*(varx*fftCov+vary))-1/(2*T)*specy./((varx*fftCov+vary).^2);

  dobj = zeros(2,1);
  dobj(1) = sum(dobdjdspec.*dfftCov_dlen)*lenx*varx;
  dobj(2) = sum(dobdjdspec.*fftCov)*varx;

else
  lenx = exp(param(1));
  varx = varargin{2};
  vary = varargin{1};
  
  [fftCov,dfftCov_dlen] = getGPSESpec(lenx,T);
  
  obj = 1/2*sum(log(varx*fftCov+vary))+1/(2*T)*sum(specy./(varx*fftCov+vary));

  dobdjdspec = 1./(2*(varx*fftCov+vary))-1/(2*T)*specy./((varx*fftCov+vary).^2);

  dobj = sum(dobdjdspec.*dfftCov_dlen)*lenx*varx;

end

obj = obj/T;
dobj = dobj/T;