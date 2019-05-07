function [Obj,dObj] = GetGPMargObj2(theta,y,Params)
  
  % function [Obj,dObj] = GetGPMargObj2(theta,x,y,Params)
  %
  % Objective function and gradients for the 'marginal' model:
  % p(x|\theta) = Norm(x;mux,varx)
  % p(y|x,\theta) = Norm(y;0,varc a^2(x)+vary)  
  
    [varx,len,mux,varc,vary,mumu,varmu,alphac,betac, ...
          alphax,betax] = UnpackParamsGPPriors(Params);
    
    tiny = 1e-5;
    vary = vary+tiny;
    
    T = length(y);
    
    x = theta(1:T);
    varc = exp(theta(T+1));
    varx = exp(theta(T+2)); 
    mux = theta(T+3);
        
    %%%%%%%%%%%%%%%%%%%%
    % Likelihood

    a = GetAmp(x(1:T));
    dadx = 1./(1+exp(-x(1:T)));
    vars = a.^2*varc+vary;

    ObjA = sum(1/2*log(vars) + 1/2*y.^2./vars);
    
    dObjA = zeros(T+3,1);
    dObjA(1:T) = dadx.*a./vars.*(1-y.^2./vars)*varc;
  
    dObjA(T+1) = -1/2*sum(a.^2*varc./vars.*(y.^2./vars-1));

    
    %%%%%%%%%%%%%%%%%%%%%
    % Prior on x
    
    dx = (x-mux);
    ObjB = 1/2*sum(dx.^2)/varx-T/2*log(varx);
    
    dObjB = zeros(T+3,1);
    
    dObjB(1:T,1) = dx/varx;
    
    dObjB(T+3,1) = -sum(dx)/varx;
    
    dObjB(T+2,1) = -1/2*sum(dx.^2)/varx-T/2;
    
    %%%%%%%%%%%%%%%%%%%%%
    % Prior on parameters
    
    ObjC = 1/2*(mux-mumu)^2/varmu+(alphac+1)*log(varc)- ...
           betac/varc+(alphax+1)*log(varx)-betax/varx;
    
    dObjC = zeros(T+3,1);
    dObjC(T+1,1) = (alphac+1)+betac/varc;
    
    dObjC(T+3,1) = (mux-mumu)/varmu;
    
    dObjC(T+2,1) = (alphax+1)+betax/varx; 
    
    %%%%%%%%%%%%%%%%%%55
    
%     Obj = ObjA;
%     dObj = dObjA;
    
%     Obj = ObjB;
%     dObj = dObjB;
    
%     Obj = ObjC;
%     dObj = dObjC;
    
    Obj = ObjA+ObjB+ObjC;
    dObj = dObjA+dObjB+dObjC;
    