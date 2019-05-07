function [a,x, Info] = InferAmplitudeGPMAP(y,Params,Opts)  
  
  % function [a,x, Info] = InferAmplitudeGPMAP(y,Params,Opts)  
  % 
  % Carries out FAST MAP Inference for the amplitudes for
  % GPPAD as described in Turner, 2010, 'Statistical Models
  % for Natural Sounds' UCL PhD Thesis. Available from
  % http://www.gatsby.ucl.ac.uk/~turner
  %
  % Copy write, Richard Turner 2010
  %   
  % INPUTS
  % y = Signal [T,1]
  % Params = structure of parameters containing
  %   len = length scale of the GP
  %   mux = shift parameter of transformed envelopes 
  %   varc = variance of carrier noise 
  %   varx = variance of the amplitude noise
  %   vary = observations noise (usually set to zero)     
  %
  % Opts = structure of options including
  % Opts.Long = [OPTIONAL] whether to return the inferences over the
  %             observed section (default) or the full set of
  %             inferences over the circularised
  %             variables too (Opts.Long=1).
  % Opts.NumIts = Total number of iterations of conjugate
  %             gradient update
  %   
  % Opts.Tx = [OPTIONAL] length of the observed+missing
  %           data can be specified by the user. Should
  %           be a power of 2 
  % Opts.tol = number of time-scales to be added to the
  %           signal to be added to avoid artifacts from
  %           circularising the signal (Tx = T+tol*len)
  % Opts.MinLength = minimum length-scale to down-sample
  %           the data to. 
  %
  % OUTPUTS
  % a = inferred modulators [1,T] (or [1,Tx] if Opts.Long==1)
  % x = inferred transforemed modulators [1,T] (or [1,Tx] if Opts.Long==1) 
  % Info = structure containing various pieces of
  %        information about the estimation including,
  % Info.Obj = Objective function 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Figure out the length of the missing latent variables
% so that they are a power of two in length (this speeds
% up the fft)

[varx,len,mux,varc,vary] = UnpackParamsGP(Params);   
T = length(y);

if isfield(Opts,'Tx')
  Tx = Opts.Tx;
else
  Tx = GetTx(T,len*Opts.tol);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Down Sampling Rates - the gradient based optimisation
% is initialised by downsampling the data and demodulating
% this downsampled version, before upsampling to initialise
% the higher sampling rate version

if isfield(Opts,'DSLength') 
  % if the user specifies how they want the downsampled
  % version to be handled e.g. when using Laplace
  DSs = GetDSs(Opts.DSLength,Opts.MinLength);
else
  DSs = GetDSs(len,Opts.MinLength);
end

Lens = len./DSs;
NumDSs = length(DSs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over down sampling rates

ParamsCur = Params;  
NumIts = ones(NumDSs,1)*ceil(Opts.NumIts/NumDSs);
Info.Obj = [];
Info.NumIts = NumIts;

fprintf('\n')
  
for n=1:NumDSs

%  disp(['Resolution ',num2str(n),'/',num2str(NumDSs)])
  %fprintf(['\rInference Resolution ',num2str(n),'/',num2str(NumDSs)]) % WW
  yDS = y(1:DSs(n):T);
  TCur = length(yDS);
  TxCur =  Tx/DSs(n);

  DimsCur = PackDimsGPFast(TCur,TxCur); 
  ParamsCur.len = Lens(n);
  
  % If there is non-stationary observations noise...
  if length(Params.vary)>1
    ParamsCur.vary = Params.vary(1:DSs(n):T);
  end

  % The first sampling rate is initialised by square and
  % low pass (or a version like this)
  if n==1
    % If there is missing data - use a longer window
    if max(ParamsCur.vary)>5*var(yDS)
      WinSz = max([floor(Lens(1)*2),4]);
    else
      WinSz = max([floor(Lens(1)/2),2]);
    end
    [xDS,aDS,cDS] = InitPADSampleNonStat(yDS,ParamsCur,WinSz);
    xDS = [xDS;Params.mux*ones(TxCur-TCur,1)];
    
%    figure
%    hold on;
%    plot(abs(yDS),'-k')
%    plot(aDS,'-r','linewidth',2)
%    keyboard
  
  end

  % The gradient based optimisation carried out using the
  % conjugate gradient algorithm
 
  [xDS,Obj] = MAPGPFast(xDS,yDS,ParamsCur,DimsCur,NumIts(n));
 
  % Upsample to initialise the next round
  if n<NumDSs
    xDS = resample(xDS,2,1);
  end
  
  % Retain the objectives per number of time-points for reference
  Info.Obj = [Info.Obj;Obj/Tx];
end 

% Return either the short version or the long version
% (with the circularised data points) depending on the
% user requests.

if isfield(Opts,'Long')
  a = log(1+exp(xDS));
  x = xDS;
else
  a = log(1+exp(xDS(1:T)));
  x = xDS(1:T);
end

fprintf('\n')