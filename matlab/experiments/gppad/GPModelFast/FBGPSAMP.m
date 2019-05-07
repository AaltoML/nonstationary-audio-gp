function [A,X,XVar,C,Params,Info] = FBGPSAMP(Y,lenx,Opts);
  
  % function [A,X,XVar,Params,Info] = FBGPSAMP(Y,lenx,Opts);
  %
  % Performs Probabilistic Amplitude Demodulation down each
  % of the columsn of a multidimensional signal Y. For
  % each column, the method infers the parameters and the
  % envelopes by sampling based inference.
  %
  % INPUTS
  % Y = data [T,D]
  % lenx = time-scale of the modulation to pick out (in
  %        samples) either [D,1] or [1] 
  % Opts = various options controlling the inference and
  %        learning   

    [T,D] = size(Y);
    C = zeros(T,D);
      
    if isfield(Opts,'Long')
      Tx = GetTx(T,lenx*Opts.tol);
      A = zeros(Tx,D);
      X = zeros(Tx,D);
      XVar = zeros(Tx,D);
    else
      A = zeros(T,D);
      X = zeros(T,D);
      XVar = zeros(T,D);
    end

 if length(lenx)==1
    lenx = lenx*ones(D,1);
  end

  Params.lenx = lenx;
  Params.varx = zeros(D,1);
  Params.varc = zeros(D,1);
  Params.mux = zeros(D,1);
  
  tic; % times the simulation duration
  for d=1:D  
    fprintf(['\nCHANNEL ',num2str(d),'/',num2str(D),'\n'])
    
    yCur = Y(:,d);
    
    % Learn the marginal parameters 
    [varcCur,varxCur,muxCur,InfoCur1] = LearnMarginalParams(yCur,Opts);    

    ParamsCur = PackParamsGP(varxCur,lenx(d),muxCur,varcCur,0); 

    % Infer the MAP amplitudes
    [aMn,xMn,xVar,InfoCur2] = InferAmplitudeGPSamp(yCur,ParamsCur,Opts);  
    
    % Output the results
    A(:,d) = sqrt(varcCur)*aMn;
    X(:,d) = xMn;
    XVar(:,d) = xVar;
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