function [A,C,varargout] = GPPAD(Y,len,varargin)
  
  % Gaussain Process Probabilistic Amplitude Demodulation
  % (GPPAD). As described in Turner and Sahani 2010,
  % 'Demodulation as an Inference Problem', TASLP and Turner
  % 2010 'Statistical models for natural sounds' UCL PhD
  % thesis.
  %
  % Various options are described below, beginning by those
  % where you just need to specify the signal and a
  % time-scale to demodulate with, moving to versions where
  % the time-scale is learned from the signal, and then more
  % complicated options like error-bar estimation or
  % cascade-demodulation with multiple modulators. For a
  % flavour of the functionality read README.m and look at
  % the various demos it refers to in the DEMO directory.
  %
  % [A,C] = GPPAD(Y,len) demodulates the columns of Y using
  % GPPAD with a time-scale given by len. The modulators and
  % carriers are returned in the columns of A and C
  % respectively.  The remaining parameters are learned from
  % the signal. len can be a scalar, in which case this
  % time-scale is used for all of the columns of Y, or a
  % vector with as many entries as there are columns in
  % Y. Various other variables can also be returned by
  % the function:
  %
  % e.g.  [A,C,Params] = GPPAD(Y,len) also returns the
  % learned parameters in the structure Params.  
  %  
  % e.g.  [A,C,Params,Info] = GPPAD(Y,len) also returns
  % algorithmic information in the structure Info like the
  % convergence of the gradient based optimisation
  %
  % e.g.  [A,C,Params,Info,X] = GPPAD(Y,len) also returns
  % the transformed modulators
  %  
  % e.g.  [A,C,Params,Info,X,VarX,Aupper,Alower] =
  % GPPAD(Y,len) also returns the error-bars on the
  % transformed envelopes (VarX) and the 3 sigma error bars
  % on the modulators (Aupper and Alower). This will take
  % more time than just vanilla inference.
  %  
  % [A,C,Params,...] = GPPAD(Y,len,LrnLen) as above if LrnLen=0,
  % but LrnLen=1 indicates that the time-scale is to be
  % learned using len as an intial starting point. All
  % the various output options above stil apply.
  %
  % [A,C,Params,...] = GPPAD(Y,Params) by specifying the
  % parameters mux,varx,varc,len in the structure Params,
  % the user indicates that these are not to be learned from
  % the signal. This is useful when they have previously
  % been learned from a training sample. It is possible to
  % fix mux, varx and varc and learn len via [A,C,Params] =
  % GPPAD(Y,Params,1).
  %   
  % [A,C,Params,...] = GPPAD(Y,len,LrnLen,Opts) is as
  % above, but the options contained in Opts modify the
  % standard options contained in the following files:
  %
  % [A,C,Params,...] = GPPAD(Y,len,2,Opts) is as above, but
  % the time-scale is learned via a grid-search where the
  % objective is evaluated at multiple time-scales and can
  % be found in the Info structure. This is useful if the
  % objective has multiple local optima and will
  % automatically detect this and give a warning. It is very
  % slow though. Opts should contain Opts.LenRng =
  % [minLen,maxLen] over which to search, and Opts.NumLens =
  % number of time-scales over which to search.
  %  
  % [A,C,Params,Info] = GPPAD(Y,[len1;len2;...;lenM])
  % (i.e. the time-scale input is a column vector of length M)
  % estimates a cascade decomposition of each of the columns
  % of Y: Y(:,d) = C(:,d)*A(:,d,1)*A(:,d,2)...*A(:,d,M) This
  % has only been implemented for the case where the
  % time-scales are known, and not when they must be
  % learned. Errorbars have not been implemented either.
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Handles the various inference/learning possibilities

  [T,D] = size(Y);
  
  % Figuring out whether the time-scale has to be learned  
  if nargin==2 
    LrnLen=0;
  else
    LrnLen = varargin{1};
  end
  
  % Figures out whether the marginal parameters are to be learned  
  if isstruct(len)  
    LrnPar = 0;
    Params=len;        
  else
    LrnPar = 1;
  end
  
  % Figures out whether the user is demanding a cascade
  if isstruct(len)  
    if size(len.len,1)>1 & D==1
      Cascade=1;
    else
      Cascade=0;
    end    
  else
    if size(len,1)>1
      Cascade=1;
    else
      Cascade=0;
    end
  end
  
  
  % Figures out whether we need to return error-bars  
  if nargout<6
    ErrBar = 0;
  else  
    ErrBar = 1;
  end 

  % Figures out whether the user has specified any
  % modification to the standard options
  if nargin>3
    OptsMod = varargin{2};
  else
    OptsMod = [];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Call the functions based on the possibilities above:
  
  % Loading and updating all the possible options
  LoadLearnMargGPPADOpts;
  LoadInferenceGPMAPOpts;
  LoadLearnLenGradAscentGPPADOpts;
  LoadErrorBarsLapOpts;
  CASGPMAPOpts;
  
  OptsNew = ModifyOpts(Opts,OptsMod);

  if Cascade==0
    % No cacade
    if LrnLen==0 & ErrBar==0
      [A,X,C,Params,Info] = FBGPMAP(Y,len,OptsNew);
    elseif LrnLen==0 & ErrBar==1
      OptsNew.LearnLengths = 0;
      [A,X,C,Params,Info,VarX,Aupper,Alower] = ...
          FBGPMAPLearnLen(Y,len,OptsNew);
    elseif LrnLen==1 & ErrBar==0
      OptsNew.LearnLengths = 1;
      [A,X,C,Params,Info] = FBGPMAPLearnLen(Y,len,OptsNew);
    elseif LrnLen==1 & ErrBar==1
      OptsNew.LearnLengths = 1;
      [A,X,C,Params,Info,VarX,Aupper,Alower] = ...
          FBGPMAPLearnLen(Y,len,OptsNew);
    elseif LrnLen==2 & ErrBar==0
      OptsNew.LearnLengths = 2;
      [A,X,C,Params,Info] = FBGPMAPLearnLen(Y,len,OptsNew);
    elseif LrnLen==2 & ErrBar==1
      OptsNew.LearnLengths = 2;
      [A,X,C,Params,Info,VarX,Aupper,Alower] = ...
          FBGPMAPLearnLen(Y,len,OptsNew);
    else
      disp('ERROR in GPPAD: Unknown method required')
    end
    
  else 
    % Cascade
    [A,X,C,Params,Info] = FBCasGPMAP(Y,len,Opts);
  end  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Handles the various output possibilities
    
  if nargout>2
    varargout(1) = {Params};
  end
  
  if nargout>3
    varargout(2) = {Info};
  end

  if nargout>4
    varargout(3) = {X};
  end

  if nargout>5
    varargout(4) = {VarX};
  end

  if nargout>6
    varargout(5) = {Aupper};
    varargout(6) = {Alower};
  end  
  
