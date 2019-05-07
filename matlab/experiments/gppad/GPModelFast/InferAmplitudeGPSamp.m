function [aMn,xMn,xVar,Info] = InferAmplitudeGPSamp(y,Params,Opts)  
  
  % function [aMn,xMn,xVar,Info] = InferAmplitudeGPSamp(y,Params,Opts)  
  % 
  % FAST MAP Inference for the amplitudes  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Figure out the length of the missing latent variables
% so that they are a power of two in length

[varx,len,mux,varc,vary] = UnpackParamsGP(Params);   
T = length(y);
Tx = GetTx(T,len*Opts.tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation

WinSz = floor(len*2);
[x,a,c] = InitPADSample(y,Params,WinSz);

x = [x;Params.mux*ones(Tx-T,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demodulation via sampling

[XSamp,Params,Info] = ESS_HMCMC_GPFAST(x,y,Params,Opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manage Output

if isfield(Opts,'Long')
  xMn = mean(XSamp,2);
  aMn = mean(log(1+exp(XSamp)),2);
  xVar = var(XSamp,0,2);
else
  xMn = mean(XSamp(1:T,:),2);
  aMn = mean(log(1+exp(XSamp(1:T,:))),2);
  xVar = var(XSamp(1:T,:),0,2);
end