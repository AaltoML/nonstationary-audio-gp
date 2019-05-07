function [Obj,dObj] = GetObjNumericalInt2(theta,y,Opts) 
  
  % function [Obj,dObj] = GetObjNumericalInt2(theta,y,Opts) 
  %  
  % Computes the objective function for initialising the
  % parameters - i.e. an approximation to the marginal
  % p(y|\theta)=\int p(y|\theta,x)p(x|\theta) from
  % numerical approximation WITH THE GRADIENTS
  %
  % INPUTS
  % theta = parameters [3,1]
  % y = data [T,1]
  % Opts = structure of options including
  %     NumXs = Number of xs to integrate over (2000)
  %     zUpper = Upper limit of x (6)
  %     zLower = Lower limit of x (-12)  
  %     ChSz = Size of data-chunks (avoids memory problems: we produce arrays of  %            size [NumXs * ChSz]) (500) 
  %
  % OUTPUTS
  % Obj = Objective function:   

  varc = exp(theta(1));
  varx = exp(theta(2));
  mux = theta(3);
    
  NumXs = Opts.NumXs; %500; % Number of xs to integrate over
  zUpper = Opts.zUpper; %7;   % Upper limit of x
  zLower = Opts.zLower; %-7; %-12  % Lower limit of x
  
  ChSz = Opts.ChSz; %1000;   % Size of data-chunks (avoids memory
                    % problems: we produce arrays of size [NumXs * ChSz])
  
  z = linspace(zLower,zUpper,NumXs);
  
  tiny = 1e-9;  % A small value to avoid vary=0

  a = GetAmp(z*sqrt(varx)+mux);
  logpz = -1/2*z.^2;

  vary = varc*a.^2+tiny;
    
  T = length(y);  
  NumCh = ceil(T/ChSz);
  
  Obj = 0;
  dObj = zeros(3,1);
  
  dadmux = 1./(1+exp(-(z*sqrt(varx)+mux)));
  dadlogvarx = 1/2*sqrt(varx)*z.*dadmux;
  
  
  for Ch = 1:NumCh
    TStart = (Ch-1)*ChSz+1;
    TEnd = min(T,Ch*ChSz);
    TCh = TEnd-TStart+1;
    
    logpxys =  -y(TStart:TEnd).^2*(1./(2*vary))+...
        repmat(logpz-1/2*log(vary),[TCh,1]);
    
    maxlogpxys = max(logpxys')';
    logpxys = logpxys-maxlogpxys*ones(1,NumXs);
    
    logpys = log(sum(exp(logpxys),2))+maxlogpxys;
    Obj = Obj+sum(logpys);

    A1 = (y(TStart:TEnd).^2)*(a.^2./vary.^2)-repmat(a.^2./vary,[TCh,1]);

    dObj(1) = dObj(1) + 1/2*varc*sum(sum(exp(logpxys).* ...
                                         A1,2)./sum(exp(logpxys),2));
 
    A2 = (y(TStart:TEnd).^2)*(a./vary.^2.*dadlogvarx)- ...
         repmat(a.*dadlogvarx./vary,[TCh,1]);

    dObj(2) = dObj(2) + varc*sum(sum(exp(logpxys).* ...
                                     A2,2)./sum(exp(logpxys),2));

    A3 = (y(TStart:TEnd).^2)*(a./vary.^2.*dadmux)- ...
         repmat(a.*dadmux./vary,[TCh,1]);

    dObj(3) = dObj(3)  + varc*sum(sum(exp(logpxys).* ...
                                      A3,2)./sum(exp(logpxys),2)) ;
       
  end

%  Prior on varc

  ObjPrior = 1/2*(theta(1)-Opts.logvarcMn)^2/Opts.logvarcVar;
  dObjPrior = (theta(1)-Opts.logvarcMn)/Opts.logvarcVar;

  Obj = Obj+T*(log(z(2)-z(1))-log(2*pi));
  
  Obj = -Obj;
  dObj = -dObj;
  
  Obj(1) = Obj(1)+ObjPrior;
  dObj(1) = dObj(1)+dObjPrior;
  