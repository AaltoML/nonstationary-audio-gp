function [x,Obj] = MAPGPFast(x,y,Params,Dims,NumIts);

% [x,Obj] = MAPGPFast(x,y,Params,Dims,NumIts);
%  
% Updates the envelopes using MAP estimation in the FAST
% version of GPPAD
%
% INPUTS
% x = transformed envelopes [Tx,1]
% y = data [T,1]
% Params = structure of parameters
% Dims = structure of Dims (containing Tx and T)  
% NumIts = number of iterations of conjugate gradients
%
% OUTPUTS
% x = updated transformed envelopes [Tx,1]
% Obj = objective function history [N,1], N>=NumIts  
  

  fftCov = GetFFTCovFast(Params.len,Dims.Tx);

  xIn=x;
  [x,Obj,in] = minimize(x,'GetGPObjFast',NumIts,y,fftCov,Params,Dims);
