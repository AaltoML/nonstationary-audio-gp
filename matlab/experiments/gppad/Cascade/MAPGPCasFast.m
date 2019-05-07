function [X,Obj] = MAPGPCasFast(X,y,Params,Dims,NumIts);
  

% function [X,Obj] = MAPGPCasFast(X,y,Params,Dims,NumIts);
%
% MAP inference of the cascade
  
fftCov = zeros(Dims.Tx,Dims.M);

for m=1:Dims.M
  fftCov(:,m) = GetFFTCovFast(Params.len(m),Dims.Tx);
end

tic;
[x,Obj,in] = minimize(X(:),'GetObjCasGPFAST',NumIts,y,fftCov,Params,Dims);
time = toc;

disp(['Duration of fine tuning ',num2str(time),' seconds'])

X = reshape(X,[Dims.Tx,Dims.M]);