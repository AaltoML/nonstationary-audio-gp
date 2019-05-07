function Obj = GetLaplaceObjGPPAD(y,Params,Opts)
  
  % function Obj = GetLaplaceObjGPPAD(y,Params,Opts)
  %
  % Laplace's approximation for the marginal likelihood of
  % the GPPAD model by breaking the data into multiple
  % chunks and treating each chunk as an independent
  % data-set
  %
  % INPUTS
  % y = signal [T,1]
  % Params = Structure of parameters including
  %     varc = carrier variance
  %     varx = transformed envelope variance
  %     mu = transformed envelope mean
  %     vary = observation noise (scalar or [T,1])
  %     len = length-scale of the modulator  
  % Opts = structure of options including
  %     NumIts = number of iterations in the inference step  
  %     tol = amount of padding (missing) data (in time-scales)  
  %     nev = number of eigenvalues to find (around 100 is a
  %           good number)  
  %     TruncPoint = Number of standard deviations to retain 
  %     MinLength = small length-scale to allow (controls
  %                 downsampling rates)
  %
  % OUTPUT
  % Obj = Objective function
    
  Opts.Long = 1; % To get the inference program to return
                 % the modulators with missing data    
  
  T = length(y);  % length of signal
     
  TCh =  Opts.TCh; % length of one chunk
  TxCh = Opts.TxCh;
    
  NumChunks = ceil(T/TCh); % number of chunks
  
  Obj = 0; % Objective

  fprintf('\n')
  
  for ch=1:NumChunks % loop over chunks
    
    % progress
    %disp(['Progress: ',num2str(ch),'/',num2str(NumChunks)])
    %str = ['\rLaplace Objective Progress: ',num2str(ch),'/',num2str(NumChunks)];
    %fprintf(str)
    
    first = (ch-1)*TCh+1;  % first sample of chunk
    last = min([T,ch*TCh]);% second sample of chunk
    
    DimsCur.T = last-first+1;
    DimsCur.Tx = TxCh;
    Opts.Tx = TxCh;
    
    yCur = y(first:last,1); % current portion of signal

    ParamsCur = Params;
    
    if length(Params.vary)>1 % for non-stationary
                             % observation noise
      ParamsCur.vary = Params.vary(first:last,1);
    end
    
    % estimate the modulator for the chunk
    [aInf,xInf, Info] = InferAmplitudeGPMAP(yCur,ParamsCur,Opts);  
    
    % compute laplace's approximation for that chunk
    ObjCur = GetLaplaceObjGPPADOneChunk(xInf,yCur,ParamsCur,DimsCur,Opts);
    
    % add this to the objective
    Obj = Obj+ObjCur;
    
  end
  
