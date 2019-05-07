function [A,X,C,Params,Info,varargout] = FBGPMAPLearnLen(Y,len,Opts);
  
  % function [A,X,Params,Info,Varx] = FBGPMAPLearnLen(Y,len,Opts);
  %
  % Performs Probabilistic Amplitude Demodulation down each
  % of the columns of a multidimensional signal Y. For each
  % column, the method infers all of the parameters and the
  % envelopes.
  %
  % INPUTS
  % Y = data [T,D]
  % lenx = time-scale of the modulation to pick out (in
  %        samples) either [D,1] or [1] 
  % Opts = various options controlling the inference and
  %        learning   
  %
  % OUTPUTS 
  % A = modulators [T,D] 
  % X = transformed modulators [T,D] 
  % C = carriers [T,D]
  % Params = Structure of parameters 
  %     varc = carrier variance
  %     varx = transformed envelope variance
  %     mu = transformed envelope mean
  %     len = length scale of the GP  
  % Info = various pieces of information about the
  %     algorithm 
  % Varx = (optional output) if requested this variable will
  %     contain the variance of the transformed modulators
  %     estimated using Laplace's approximation [T,D]
  % Aupper = (optional output) if requested this variable
  %     will contain the upper error-bars on the modulators
  %     estimated using Laplace's approximation [T,D]
  % Alower = (optional output) if requested this variable
  %     will contain the lower error-bars on the modulators
  %     estimated using Laplace's approximation [T,D]

    if nargout==8
      Opts.Errorbars=1;
    else
      Opts.Errorbars=0;
    end

    
    if isstruct(len)
      Params = len;
      len = Params.len;
      LrnPar = 0;
    else
      LrnPar = 1;
    end  

    
    [T,D] = size(Y);
    C = zeros(T,D);
    
    A = zeros(T,D);
    X = zeros(T,D);

    if Opts.Errorbars==1
      Varx = zeros(T,D);
      Aupper = zeros(T,D);
      Alower = zeros(T,D);
    end
    
    % Preallocate the parameters
    
    % If only one time-scale is specified use the same
    % one for each channel
    if length(len)==1
      len = len*ones(D,1);
    elseif length(len)==0&Opts.LearnLengths==2
      len = 1;
    else
    end

    
      if exist('Params') 
      if isfield(Params,'vary')
	if D>size(Params.vary,2)
	  Params.vary = Params.vary*ones(1,D);
	end
      end
    else
      Params.vary = 1e-6*ones(T,D);
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

    % Learn the time-scales if the user requests
    if Opts.LearnLengths ==1
%      [len,InfoCur2] = LearnLengthsGPPAD(len(d),yCur,ParamsCur,Opts);
      %[len,InfoCur2] = LearnLengthsGPPADFast(len(d,:),yCur,ParamsCur,Opts);
       [len_,InfoCur2] = GetLenGradAscent(len(d,:),yCur,ParamsCur,Opts);
       ParamsCur.len = len_;
    elseif Opts.LearnLengths == 2
      if size(Opts.LenRng,1)>1
        LenRng = Opts.LenRng(d,:);
      else
        LenRng = Opts.LenRng;
      end
       [len_,InfoCur2] = GetLenGridSearch(LenRng,yCur,ParamsCur,Opts);
       ParamsCur.len = len_;
    else
      InfoCur2 = [];
    end

    % Infers the modulators and the error-bars if requested
    if Opts.Errorbars==1
      % Infer the MAP amplitudes and Error-bars
      [a,x,varx,InfoCur3] = InferAmplitudeErrorBarsGPMAP(yCur,ParamsCur,Opts);
    else
      % Infer the MAP amplitudes only
      [a,x,InfoCur3] = InferAmplitudeGPMAP(yCur,ParamsCur,Opts);  
      InfoCur3 = [];
    end
        
    % Output the results
    A(:,d) = sqrt(varcCur)*a;
    X(:,d) = x;
    C(:,d) = Y(:,d)./A(1:T,d);
        
    if Opts.Errorbars==1
      Varx(:,d) = varx;
      Aupper(:,d) = sqrt(varcCur)*log(1+exp(x+3*sqrt(varx)));
      Alower(:,d) = sqrt(varcCur)*log(1+exp(x-3*sqrt(varx)));
    end
    
    Params.varc(d) = varcCur;
    Params.varx(d) = varxCur;
    Params.mux(d) = muxCur;
    Params.len(d) = ParamsCur.len;
    
    Info{d,1} = InfoCur1;
    Info{d,2} = InfoCur2;
    Info{d,3} = InfoCur3;
  end
  
  % Duration of simulation
  time = toc;

  if time>60
    timestr = [num2str(time/60),' minutes'];
  else
    timestr = [num2str(time),' seconds'];
  end
  disp(['Duration of estimation ',timestr])
  
  if Opts.Errorbars==1
    varargout(1) = {Varx};
    varargout(2) = {Aupper};
    varargout(3) = {Alower};
  end
