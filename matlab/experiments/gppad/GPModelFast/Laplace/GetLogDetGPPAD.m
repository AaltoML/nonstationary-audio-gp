function LogDet = GetLogDetGPPAD(x,y,Params,Dims,Opts)
  
  % function LogDet = GetLogDetGPPAD(x,y,Params,Dims,Opts)
  %
  % Approximate log det using Laplace's approximation for
  % the GPPAD model   
  
[d,h] = GetHessian(x,y,Params,Dims);
[xv,lmb] = GetEigSigDSig(h,d,Opts);

lmb(lmb<0)=0;

LogDet = -sum(log(h(:)))+sum(log(abs(1+lmb)));

