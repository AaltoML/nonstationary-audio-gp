function [Len,Objs,Info] = GetLaplaceObjGPPADMulti(Len,y,Params,Opts)
  
  % function [Len,Obj,Info] = GetLaplaceObjGPPADMulti(Len,y,Params,Opts)
  %
  % Evaluates the Laplace Objective at a selection of
  % length-scales specified in Len. The chunk length is
  % set using the smallest of the specified length-scales
  % and is common for all the objective evaluations.
    
  % Figuring out the chunk length  
       
  % The maximum size of Tx is set by the smallest length-scale
  TxMax = 2*pi*min(Len)*Opts.nev/Opts.TruncPoint;  % Opts.TruncPoint>8  

  % Coverting Tx to the nearest smaller power of 2
  TxMax = 2^ceil(log2(TxMax));
  
  % Computing the maximum size of the observed data by
  % removing the extra bit
  TMax = floor(TxMax-Opts.tol*min(Len));
  
  % Chunk size (observed and unobserved)
  Opts.TCh = TMax;
  Opts.TxCh = GetTx(Opts.TCh,max(Len)*Opts.tol);

  % Evaluate the objective at the prespecified points
  L = length(Len);
  
  Objs = zeros(1,L);
  ParamsCur = Params;

  Opts.DSLength = min(Len); % ensures the same inference
                            % initialisation is used for
                            % all of the length-scales
  
  % Loop over length-scales
  for l = 1:L
    ParamsCur.len = Len(l);
    Objs(l) = GetLaplaceObjGPPAD(y,ParamsCur,Opts);
  end
  
  % Information
  Info.TCh = Opts.TCh;
  Info.TxCh = Opts.TxCh;
  [val,pos] = max(Objs);
  
  if pos==1
    % optimum is at a lower length-scale than those in Len
    Info.MaxLoc = -1;
  elseif pos == L
    % optimum is at a larger length-scale than those in Len
    Info.MaxLoc = 1;
  else
    % optimum is within the specified length-scales
    Info.MaxLoc = 0;
  end
  