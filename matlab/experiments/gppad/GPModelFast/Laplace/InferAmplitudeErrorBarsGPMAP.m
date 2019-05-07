function [a,x,Varx,Info] = InferAmplitudeErrorBarsGPMAP(y,Params,Opts)

  % function Varx = GetErrorBarsGPPADCh(y,Params,Opts)
  %
  % Gets the MAP estimate error-bars for the transformed
  % envelopes using the Laplace approximation for LONG DATA
  % by chunking up the data into sections. The overlapping
  % sections are windowed to ensure smoothness across the
  % different sections. The window means the algorithm has
  % linear cost in time.
  %  
  % INPUTS
  % x = MAP transformed envelopes [Tx,1]
  % y = data [T,1]
  % Params = Structure of parameters 
  %     varc = carrier variance
  %     varx = transformed envelope variance
  %     mu = transformed envelope mean
  %     len = length scale of the GP  
  %     vary = observation noise (scalar or [T,1])  
  %
  % Opts = structure of options
  %     tol = amount of padding (missing) data (in time-scales)  
  %     nev = number of eigenvalues to find (around 100 is a
  %           good number)  
  %     TruncPoint = Number of standard deviations to retain 
  %     Overlap = amount to overlap segments (in time-scales)
  %     MinLength = small length-scale to allow (controls
  %                 downsampling rates)
  %  
  % OUTPUTS
  % a = modulator [T,1]
  % x = transformed modulator [T,1]  
  % Varx = estimated posterior marginal variances of
  %        transformed envelopes [T,1]
  % Info = structure containing information about the
  %        inference and error-bar estimation  

T = length(y);    

% Set the chunk length according to the number of
% eigenvalues we're retaining and make it a power of 2
TxCh = floor(2*pi*Opts.nev*Params.len/Opts.TruncPoint);
TxCh = 2^floor(log2(TxCh));  
    
% Get the observed data chunk length by subtracting off
% the length of padded data 
TCh = TxCh-ceil(Params.len*Opts.tol);

% Compute the overlap of the segments
Overlap = ceil(Params.len*Opts.Overlap);

% Number of chunks 
NumChunks = ceil(T/(TCh-Overlap));

% Display some information
fprintf(['\n','Chunk length: ',num2str(TCh)])
fprintf(['\n','Number of Chunks: ',num2str(NumChunks),'\n'])

% Check for problems
if Overlap>TCh
  disp('Error: Overlap is longer than the chunk length')
  disp('Reduce Overlap or TruncPoint, or increase nev')
end

if TCh<0
  disp('Error: Observed Chunk length negative')
  disp('Reduce tol or TruncPoint, or increase nev')
end

OptsCur.Tx = TxCh;
OptsCur.nev = Opts.nev;
OptsCur.MinLength = Opts.MinLength;
OptsCur.NumIts = Opts.NumIts;
OptsCur.Long = 1;

DimsCur.Tx = TxCh;

Taper = cos(2*pi*[1:Overlap]'/(4*Overlap)).^2;
WindowStart = [ones(TCh-Overlap,1);Taper];
Window = [Taper(end:-1:1);ones(TCh-2*Overlap,1);Taper];
WindowEnd = [Taper(end:-1:1);ones(TCh-Overlap,1)];

x = zeros(T,1);
Varx = zeros(T,1);
Borders = zeros(NumChunks,1);

% Loop over chunks
for Ch = 1:NumChunks

  fprintf(['Chunk: ',num2str(Ch),'/',num2str(NumChunks),'\n'])
  
  % times in y
  tstart = 1+(Ch-1)*(TCh-Overlap);
  tend = min(tstart+TCh-1,T);
  
  % To deal with the last chunk being of different length
  TChCur = tend-tstart+1;
  DimsCur.T = TChCur;
  
  % Data for the current chunk
  yCur = y(tstart:tend);

  % Parameters for the current chunk
  ParamsCur = Params;
  
  if length(Params.vary)>1
    ParamsCur.vary = Params.vary(tstart:tend);
  end
  
  % Infer transformed envelopes for chunk
  [aCur,xCur,InfoCur] = InferAmplitudeGPMAP(yCur,ParamsCur,OptsCur);  

  % Compute error-bars for chunk
  VarxCur = GetErrorBarGPPAD(xCur,yCur,ParamsCur,DimsCur,OptsCur);
  
  % Read the current chunk into the transformed modulators and
  % error-bars. These are windowed to ensure smoothness
  % across the overlapping region

  if NumChunks==1 % if there is just one chunk
    Varx(tstart:tend) = VarxCur(1:TChCur);
    x(tstart:tend) = xCur(1:TChCur);
  elseif Ch>1 & Ch<NumChunks;
    Varx(tstart:tend) = Varx(tstart:tend)+VarxCur(1:TChCur).*Window(1:TChCur);
    x(tstart:tend) = x(tstart:tend)+xCur(1:TChCur).*Window(1:TChCur);
  elseif Ch==1 
    Varx(tstart:tend) = Varx(tstart:tend)...
        +VarxCur(1:TChCur).*WindowStart(1:TChCur);
    x(tstart:tend) = x(tstart:tend)...
        +xCur(1:TChCur).*WindowStart(1:TChCur);
  elseif Ch==NumChunks
    Varx(tstart:tend) = Varx(tstart:tend)...
        +VarxCur(1:TChCur).*WindowEnd(1:TChCur);
    x(tstart:tend) = x(tstart:tend)...
        +xCur(1:TChCur).*WindowEnd(1:TChCur);
  else
  end

  Borders(Ch) = tend-Overlap/2; 
end

Info.Borders = Borders;

a = log(1+exp(x));