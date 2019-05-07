function [len,Info] = GetLenGradAscent(len,y,Params,Opts)


  % function [Len,Obj,Info] = GetLenGradAscent(Len,y,Params,Opts)
  %
  % Finds the length-scale by stochastic gradient ascent  
    
  
  T = length(y);
  
  NumItsLL = Opts.NumItsLL;
  
  loglenLR0 = Opts.loglenLR0; % initial step size

  Nhalf = NumItsLL/2; % transition point from 'search' to 'fine-tuning'
  
  % Learning rate decreases over time
  loglenLR = loglenLR0./(1+[0:NumItsLL]/Nhalf); 

  m = Opts.m;
  
  LastUp = 0;
  
  Info.LenHist = zeros(NumItsLL,1);
  Info.ObjHist = zeros(NumItsLL,1);
  Info.dObjHist = zeros(NumItsLL,1);
  Info.Ts = zeros(NumItsLL,2);
  Info.loglenLR = loglenLR;

  disp('Learning the time-scale using stochastic gradient ascent')
  
  for it=1:NumItsLL  
    
    disp(['Progress: ',num2str(it),'/',num2str(NumItsLL)])
    disp(['Time-scale: ',num2str(len)])
    disp(['Step-size: ',num2str(loglenLR(it))])
%    disp(['Objective: ',num2str(Obj)])
    
    TxMax = 2*pi*len*Opts.nev/Opts.TruncPoint;  % Opts.TruncPoint>8  
    TxMax = 2^ceil(log2(TxMax));
    TMax = floor(TxMax-Opts.tol*len);
    
    if TMax<T
      tstart = 1+floor(rand*(T-TMax)); %uniform(1,T-Tmax)
      tend = tstart-1+TMax;
    else
      tstart = 1;
      tend = T;
    end
    
    yCur = y(tstart:tend);
  
    % Evaluate the objective function twice to figure out
    % an empirical gradient
    
    lenplusdelta = max([len+1,len*1.01]);
    Len = [len,lenplusdelta]; 
    
    [Len,Objs,InfoCur] = GetLaplaceObjGPPADMulti(Len,yCur,Params,Opts);

    % If there is a prior over the length-scale
    if isfield(Params,'mulen')
      Objs = Objs-(Len-Params.mulen).^2/(2*Params.varlen);
    end
    
    dObj = diff(Objs);
    
    % Update the length-scale based on the gradient  
    %dloglen = sign(dObj)*log(dlen(it))+m*LastUp;
    dloglen = sign(dObj)*loglenLR(it)+m*LastUp;
    
    loglen = log(len)+dloglen;
    
    LastUp = dloglen;
    
    len = exp(loglen);

    Info.LenHist(it) = len;
    Info.ObjHist(it) = Objs(1);
    Info.dObjHist(it) = dObj;
    Info.Ts(it,:) = [tstart,tend];
  end
    

  