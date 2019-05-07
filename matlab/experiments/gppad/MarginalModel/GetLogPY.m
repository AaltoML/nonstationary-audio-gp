function logPY = GetLogPY(y,varx,varc,mux,Opts) 
  
  % function logPY = GetLogPY(y,varx,varc,mux,Opts) 
  %  
  % Returns values for the marginal distribution over the data

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
  
  logPY = zeros(T,1);
  
  for Ch = 1:NumCh
    TStart = (Ch-1)*ChSz+1;
    TEnd = min(T,Ch*ChSz);
    TCh = TEnd-TStart+1;
    
    logpxys =  -y(TStart:TEnd).^2*(1./(2*vary))+...
        repmat(logpz-1/2*log(vary),[TCh,1]);
    
    maxlogpxys = max(logpxys')';
    logpxys = logpxys-maxlogpxys*ones(1,NumXs);
    
    logPY(TStart:TEnd) = log(sum(exp(logpxys),2))+maxlogpxys;

  end
  
   
