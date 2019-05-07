function [Obj,dObj] = GetGPMargObj(x,y,Params)
  
  % function [Obj,dObj] = GetGPMargObj(x,y,Params)
  %
  % Objective function and gradients for the 'marginal' model:
  % p(x|\theta) = Norm(x;mux,varx)
  % p(y|x,\theta) = Norm(y;0,varc a^2(x)+vary)  
  

    [varx,len,mux,varc,vary,mumu,varmu,alphac,betac, ...
          alphax,betax] = UnpackParamsGPPriors(Params);
    
    tiny = 1e-5;
    vary = vary+tiny;
    
    T = length(y);
    
    %%%%%%%%%%%%%%%%%%%%
    % Likelihood

    a = GetAmp(x(1:T));
    dadx = 1./(1+exp(-x(1:T)));
    vars = a.^2*varc+vary;

    ObjA = sum(1/2*log(vars) + 1/2*y.^2./vars);
    
    dObjA = dadx.*a./vars.*(1-y.^2./vars)*varc;
  
    %%%%%%%%%%%%%%%%%%%%%
    % Prior on x
    
    dx = (x-mux);
    ObjB = 1/2*sum(dx.^2)/varx-T/2*log(varx);
    
    dObjB = dx/varx;
    
    %%%%%%%%%%%%%%%%%%55
    
%     Obj = ObjA;
%     dObj = dObjA;
    
%     Obj = ObjB;
%     dObj = dObjB;
    
    Obj = ObjA+ObjB;
    dObj = dObjA+dObjB;
    