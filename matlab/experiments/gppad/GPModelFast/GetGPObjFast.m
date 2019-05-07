function [Obj,dObjdx] = GetGPObjFast(x,y,fftCov,Params,Dims);
  
  % function [Obj,dObjdx] = GetGPObjFast(x,y,fftCov,Params,Dims);
  %
  %
  % Computes the transformed envelope derivatives for GPPAD
  % QUICKLY
  
    [varx,len,mux,varc,vary] = UnpackParamsGP(Params); 

    tiny = 1e-5;
    vary = vary+tiny;

    [T,Tx] = UnpackDimsGPFast(Dims); 

    %%%%%%%%%%%%%%%%%%%%
    % Likelihood

    a = GetAmp(x(1:T));
    dadx = 1./(1+exp(-x(1:T)));

    vars = a.^2*varc+vary;

    ObjA = sum(1/2*log(vars) + 1/2*y.^2./vars);
    dObjAdx = dadx.*a./vars.*(1-y.^2./vars)*varc;

    %%%%%%%%%%%%%%%%%%%%%
    % Prior
    xExt = x-mux;
    
    fftxExt = fft(xExt);
    
    InvCovX = ifft(fftxExt./fftCov);
    
    ObjB = 1/(2*varx)*InvCovX'*xExt;
    
    ObjB = real(ObjB);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    % InvCovX2 = ifft(fftxExt./fft(x));
    % A2 = 1/2*InvCovX'*x;
    % A3 = sum(x)/fftCov(1); 
    % A4 = (T-1)/fftCov(1);
    % ObjB2 = prea*A2-sh*prea* A3+prea*sh^2*A4;
    
    % [ObjB-ObjB2]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    dObjBdx = real(InvCovX)/varx;

    %%%%%%%%%%%%%%%%%%%%%
    % Output
    
    Obj = (ObjA + ObjB);
    dObjdx = dObjBdx;
    dObjdx(1:T) = dObjdx(1:T) + dObjAdx;

    %Obj = (ObjA)/T;
    %dObjdx = [dObjAdx;zeros(T-2,1)]/T;
    
    % Obj = ObjB/T;
    % dObjdx = dObjBdx/T;