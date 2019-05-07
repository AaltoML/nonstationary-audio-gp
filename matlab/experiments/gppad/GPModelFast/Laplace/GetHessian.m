function [d,h] = GetHessian(x,y,Params,Dims);
  
  % function [d,h] = GetHessian(x,y,Params,Dims);
  %
  % Returns the second derivatives of the prior (in
  % frequency space) and likelihood for the cascade model 
  %
  % INPUTS
  %
  % x = MAP transformed envelopes [TT,1]
  % y = data [T,1]
  % Params = structure of parameters
  % Dims = structure of Dims
  % 
  % OUTPUTS
  % d = second derivative of the likelihood [Tx,1]  
  % h = 1/second derivative of the prior [Tx,1]
    
    [T,Tx] = UnpackDimsGPFast(Dims); 
    [varx,len,mux,varc,vary] = UnpackParamsGP(Params); 

    tiny = 1e-5;
    vary = vary+tiny;

%%%%%%%%%%%%%%%%%%%%
% Likelihood

a = GetAmp(x(1:T));
dadx = 1./(1+exp(-x(1:T)));
d2adx2 = 1/4*1./cosh(x(1:T,:)/2).^2;

vars = a.^2*varc+vary;
dvarsdx = 2*varc*a.*dadx;

d = (d2adx2.*a+dadx.^2)./vars.*(1-y.^2./vars)*varc ...
    - dadx.*a./(vars.^2).*(1-y.^2./vars)*varc.*dvarsdx ...
    + dadx.*a./vars.*y.^2./(vars.^2)*varc.*dvarsdx;


%%%%%%%%%%%%%%%%%%%%%
% Prior

h = GetFFTCovFast(len,Tx)*varx;
