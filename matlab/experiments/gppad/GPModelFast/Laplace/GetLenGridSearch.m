function [len,Info] = GetLenGridSearch(LenRng,y,Params,Opts)

  % function [len,Info] = GetLenGridSearch(LenRng,y,Params,Opts)
  % 
  % Finds a time-scale by performing a grid serach over
  % the time-scales between the values specified in
  % LenRng;
    
  if isfield(Opts,'NumLens')
    NumLens = Opts.NumLens;
  else
    NumLens = [];
  end

  disp('Learning time-scale by a grid search')  

  [Objs,Lens] = GetLaplaceObjGPPADStitched(LenRng, ...
                                           NumLens,y,Params,Opts);
  
  Info.Objs = Objs;
  Info.Lens = Lens;
  
  [val,pos] = max(Objs);
  
  len = Lens(pos);
    
  K = length(Objs);
  
  if pos==K||pos==1
    disp(['Warning: optimal time-scale is at the edge of the range'])
  end
  
  if Objs(1)>Objs(2)
    LenMax = Lens(1);
  else
    LenMax = [];
  end
  
  for k=2:K-1
    if Objs(k)>Objs(k-1)&Objs(k)>Objs(k+1)
      LenMax = [LenMax;Lens(k)];
    end
  end

  if Objs(K)>Objs(K-1)
        LenMax = [LenMax;Lens(K)];
  end

  if length(LenMax)>1
    disp('Warning: There are multiple local optima')
    disp(['Optima found at: ',num2str(LenMax')])
  end
  
  Info.LenMax = LenMax;