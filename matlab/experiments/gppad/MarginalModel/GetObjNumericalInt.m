function Obj = GetObjNumericalInt(y,varx,varc,mux,Opts) 
  
  % function Obj = GetObjNumericalInt(y,vara,varc,mux,Opts) 
  %  
  % Computes the objective function for initialising the
  % parameters - i.e. an approximation to the marginal
  % p(y|\theta)=\int p(y|\theta,x)p(x|\theta) from
  % numerical approximation
  %
  % INPUTS
  % y = data [T,1]
  % varx = transformed envelope amplitude marginal variance
  % varc = carrier variance
  % mux = mean of the transformed envelope
  % Opts = structure of options including
  %     NumXs = Number of xs to integrate over (2000)
  %     zUpper = Upper limit of x (6)
  %     zLower = Lower limit of x (-12)  
  %     ChSz = Size of data-chunks (avoids memory problems: we produce arrays of
  %            size [NumXs * ChSz]) (500) 
  %
  % OUTPUTS
  % Obj = Objective function:   

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

  end
  
   
  Obj = Obj+T*(log(z(2)-z(1))-log(2*pi));