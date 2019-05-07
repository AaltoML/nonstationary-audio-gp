function [A,X,C,Params,Info] = FBGPMAP(Y,len,Opts);
  
  % function [A,X,Params,Info] = FBGPMAP(Y,len,Opts);
  %
  % Performs Probabilistic Amplitude Demodulation down each
  % of the columns of a multidimensional signal Y. For each
  % column, the method infers the marginal parameters and
  % the envelopes by MAP inference. If the user wishes to
  % learn the time-scales of the modulators, or to have
  % errorbars on the MAP inferences, please use
  % FBGPMAPLearnLen.m
  %
  % INPUTS
  % Y = data [T,D]
  % len = time-scale of the modulation to pick out (in
  %        samples) either [D,1] or [1] 
  % Opts = various options controlling the inference and
  %        learning  
  %
  % OUTPUTS 
  % A = modulators [T,D] (or [Tx,D] if user requests in
  %                       the options)
  % X = transformed modulators [T,D] (or [Tx,D] if user requests in
  %                       the options) 
  % C = carriers [T,D]
  % Params = Structure of parameters 
  %     varc = carrier variance
  %     varx = transformed envelope variance
  %     mu = transformed envelope mean
  %     len = length scale of the GP  
  % Info = various pieces of information about the
  % algorithm 

    if isstruct(len)
      Params = len;
      len = Params.len;
      LrnPar = 0;
    else
      LrnPar = 1;
    end  
    
    [T,D] = size(Y);
 
    if exist('Params') 
      if isfield(Params,'vary')
	if D>size(Params.vary,2)
	  Params.vary = Params.vary*ones(1,D);
	end
      end
    else
      Params.vary = zeros(T,D);
    end
    

    C = zeros(T,D);
    
    % If the user requests the inferences for the entire
    % circularised data-set to be returned, then
    % preallocated accordingly. Also tests to see whether
    % the user specifies the length of the augmented
    % dataset (Tx).
    if isfield(Opts,'Long')

      if isfield(Opts,'Tx')
        Tx = Opts.Tx;
        A = zeros(Tx,D);
        X = zeros(Tx,D);
      else
        Tx = GetTx(T,len*Opts.tol);
        A = zeros(Tx,D);
        X = zeros(Tx,D);
      end
    else
      A = zeros(T,D);
      X = zeros(T,D);
    end

 if length(len)==1
    len = len*ones(D,1);
  end

  Params.len = len;
  
  if LrnPar==1
    Params.varx = zeros(D,1);
    Params.varc = zeros(D,1);
    Params.mux = zeros(D,1);
  else
    if length(Params.varx)~=D
      Params.varx = Params.varx*ones(D,1);
      Params.varc = Params.varc*ones(D,1);
      Params.mux = Params.mux*ones(D,1);
    end
  end
    
  tic; % times the simulation duration
  for d=1:D  
    fprintf(['\nCHANNEL ',num2str(d),'/',num2str(D),'\n'])
    
    yCur = Y(:,d);
    
    if LrnPar==1
      % Learn the marginal parameters 
      [varcCur,varxCur,muxCur,InfoCur1] = LearnMarginalParams(yCur,Opts);    
    else
      varcCur = Params.varc(d);
      varxCur = Params.varx(d);
      muxCur = Params.mux(d);
      InfoCur1 = [];
    end
    
    varyCur = Params.vary(:,d);
    ParamsCur = PackParamsGP(varxCur,len(d),muxCur,varcCur,varyCur); 
    
    % Infer the MAP amplitudes
    [a,x, InfoCur2] = InferAmplitudeGPMAP(yCur,ParamsCur,Opts);  
    
    % Output the results
    A(:,d) = sqrt(varcCur)*a;
    X(:,d) = x;
    C(:,d) = Y(:,d)./A(1:T,d);
    
    Params.varc(d) = varcCur;
    Params.varx(d) = varxCur;
    Params.mux(d) = muxCur;
    
    Info{d,1} = InfoCur1;
    Info{d,2} = InfoCur2;
  end
  
  % Duration of simulation
  time = toc;

  if time>60
    timestr = [num2str(time/60),' minutes'];
  else
    timestr = [num2str(time),' seconds'];
  end
  disp(['Duration of estimation ',timestr])