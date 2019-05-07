function KL = GetKLPY2(varx1,varc1,mux1,varx2,varc2,mux2,Opts) 
  
  % function KL = GetKLPY(y,varx1,varc1,mux1,varx2,varc2,mux2,Opts) 
  %  
  % Estimate the KL of two marginal distributions

  NumYs = Opts.NumYs; 
  yUpper = Opts.yUpper; 
    
  NumXs = Opts.NumXs; %500; % Number of xs to integrate over
  zUpper = Opts.zUpper; %7;   % Upper limit of x
  zLower = Opts.zLower; %-7; %-12  % Lower limit of x
  
  ChSz = Opts.ChSz; %1000;   % Size of data-chunks (avoids memory
                    % problems: we produce arrays of size [NumXs * ChSz])
  
  z = linspace(zLower,zUpper,NumXs);
  y = linspace(0,yUpper,NumYs)';
  
  tiny = 1e-9;  % A small value to avoid vary=0

  a1 = GetAmp(z*sqrt(varx1)+mux1);
  logpz1 = -1/2*z.^2-1/2*log(2*pi);
  vary1 = varc1*a1.^2+tiny;

  a2 = GetAmp(z*sqrt(varx2)+mux2);
  logpz2 = -1/2*z.^2-1/2*log(2*pi);
  vary2 = varc2*a2.^2+tiny;

  T = length(y);  
  NumCh = ceil(T/ChSz);
  
  Obj = 0;
  
  logPY1 = zeros(T,1);
  logPY2 = zeros(T,1);
  
  for Ch = 1:NumCh
    TStart = (Ch-1)*ChSz+1;
    TEnd = min(T,Ch*ChSz);
    TCh = TEnd-TStart+1;
  
    logpxys1 =  -y(TStart:TEnd).^2*(1./(2*vary1))+...
        repmat(logpz1-1/2*log(2*pi*vary1),[TCh,1]);
    
    maxlogpxys1 = max(logpxys1')';
    logpxys1 = logpxys1-maxlogpxys1*ones(1,NumXs);
    
    logPY1(TStart:TEnd) = log(sum(exp(logpxys1),2))+maxlogpxys1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    logpxys2 =  -y(TStart:TEnd).^2*(1./(2*vary2))+...
        repmat(logpz2-1/2*log(2*pi*vary2),[TCh,1]);
    
    maxlogpxys2 = max(logpxys2')';
    logpxys2 = logpxys2-maxlogpxys2*ones(1,NumXs);
    
    logPY2(TStart:TEnd) = log(sum(exp(logpxys2),2))+maxlogpxys2;
    
  end

  dy = y(2)-y(1);
  dz = z(2)-z(1);

  % Factor of two comes from the fact that we only
  % integrated from zero upward and, by symmetry, the
  % result is twice this value
  
  temp = exp(logPY1).*(logPY1-logPY2);
  KL = (2*sum(temp)-temp(1))*dy*dz;

  if 2*sum(exp(logPY1))*dy*dz<0.9
    disp('Warning: KL limits MAY be causing inaccurate values')
  end

  if KL<0
    disp('Error: KL limits ARE causing negative KLs')
    KL = 0;
  end


  