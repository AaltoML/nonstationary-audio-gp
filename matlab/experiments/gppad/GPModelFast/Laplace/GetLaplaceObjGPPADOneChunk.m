function Obj = GetLaplaceObjGPPADOneChunk(x,y,Params,Dims,Opts)
  
  % function Obj = GetLaplaceObjGPPADOneChunk(x,y,Params,Dims,Opts)
  %
  % Laplace's approximation for the marginal likelihood of
  % the GPPAD model for one chunk of data only
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contribution from the determinant of the Hessian
  
[d,h] = GetHessian(x,y,Params,Dims);
[xv,lmb] = GetEigSigDSig(h,d,Opts);

lmb(lmb<0)=0;

Obj0 = -1/2*sum(log(abs(1+lmb)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff needed for the next two components

[varx,len,mux,varc,vary] = UnpackParamsGP(Params); 

tiny = 1e-5;
vary = vary+tiny;
    
[T,Tx] = UnpackDimsGPFast(Dims); 

fftCov = GetFFTCovFast(len,Tx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contribution from the likelihood

a = GetAmp(x(1:T));

vars = a.^2*varc+vary;
    
Obj1 = -1/2*sum(log(vars) + y.^2./vars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contribution from the Prior

xExt = x-mux;
   
fftxExt = fft(xExt);

InvCovX = ifft(fftxExt./fftCov);
    
Obj2 = -1/(2*varx)*InvCovX'*xExt;

Obj2 = real(Obj2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total contribution
Obj = Obj0+Obj1+Obj2;
