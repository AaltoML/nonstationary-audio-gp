function [H,info] = tnmf_inf(A,W,H,lenx,mux,varx,vary,varargin)

% function [H,info] = tnmf_inf(A,W,H,lenx,mux,varx,vary,varargin)
%
% Performs inference for temporal NMF on the data matrix A. That is,
% it approximates A via A \approx H*W. The NMF variant is trained
% using conjugate gradients and the temporal and spectral basis
% function are learned in a joint optimisation.
%
% The cost function used in the likelihood is:
% \sum_{d,t} A_{d,t}./Ahat_{t,d} + \sum_{d,t} \log(Ahat_{t,d}) 
%
% Each of the temporal basis functions in NMF is given by: 
% h(t) = exp(x(t))  x(t) ~ mux + GP(lenx,varx)
%
%
% INPUTS
% Model parameters:
%   A = data [T,D]
%   W = initial spectral weight vector [K,D]
%   lenx = squared exponential length-scale [1,K]
%   mux = steady state mean of the log temporal priors [1,K]
%   varx = steady state variance of the log temporal priors [1,K]
%   H = initial temporal weight vector [K,T]
%   vary = observation noise [T,D] (set to zero if signal is noiseless)
%
% Optional algorithmic options:
% numIts = number of iterations per optimisation batch
% progress_chunk = number of optimisation batches 
%
% OUTPUTS
%   H = inferred temporal weight vector [K,T]
%   info = structure containing information about the optimisation including
%     Obj = objective function
%     it = number of completed iteration steps per batch

if nargin>7 & isfield(varargin{1},'numIts')
  numIts = varargin{1}.numIts;
else
  numIts = 100;
end

if nargin>7 & isfield(varargin{1},'progress_chunk')
  progress_chunk = varargin{1}.progress_chunk;
else
  progress_chunk = 25;
end

if nargin>7 & isfield(varargin{1},'tol')
  tol = varargin{1}.tol;
else
  tol = 11;
end

[K,D] = size(W);
[T,K] = size(H);
X = zeros(T,K);

disp('tnmf: converting to basis function representation')
for k=1:K
  % convert from 
  logH = log(H(:,k))-mux(k);
  [X(:,k),infoTrans] = fitSEGP_BF(logH,lenx(k),varx(k));

  % logHfit = basis2SEGP(X(:,k),lenx(k),varx(k),9);
  % figure
  % hold on
  % plot(logH,'-k')
  % plot(logHfit,'-b')
  % keyboard
end

L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];

x = X(:);

for l=1:L
  
  % save old parameters
  Hold = H;
  
  % run conjugate gradient update
  tic;
  
  [x, ObjCur, itCur] = minimize(x,'getObj_nmf_log_norm_inf', ...
				      numIts(l),A,W,vary,lenx,mux,varx,tol);
  timCur = toc;

  % Pull out H 
  X = reshape(x,[T,K]);

  H = zeros([T,K]);

  for k=1:K
    tau = ceil(lenx(k)/sqrt(2)*tol);
    bh = sqrt(varx(k)*sqrt(2)/(lenx(k)*sqrt(pi))); % basis function height to
						   % get correct variance
    bas = bh*exp(-1/(lenx(k)^2)*([-tau:1:tau]').^2);
    
    H(:,k) =  exp(conv(X(:,k),bas,'same')+mux(k));
  end

  % Store objective and iteration information
  Obj = [Obj;ObjCur];
  it = [it;itCur];
  tim = [tim;timCur];
  dH = sqrt(sum((H(:)-Hold(:)).^2)/sum(H(:).^2));

  Ahat = H*W;
  snrChan = 10*log10(mean(A.^2,1))-10*log10(mean((A-Ahat).^2));
  
  % Display some information to the user
  str1 = ['Progress ',num2str(l),'/',num2str(L),];
  str2 = ['Obj ',num2str(ObjCur(end))];

  if length(ObjCur)>1
    str3 = ['dObj ',num2str(diff(ObjCur(end-1:end)))];
  else
    str3 = ['dObj ','only one value'];
  end
  str4 = ['time ',num2str(ceil(timCur/6)/10),'mins'];
  str5 = ['total time ',num2str(ceil(sum(tim)/6)/10),'mins'];
  str6 = ['dH ',num2str(round(dH*1000)/10),'%%'];
  str7 = ['A snr ',num2str(round(mean(snrChan)*1000)/1000)];

  str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6,str_space,str7])

end

fprintf('\n')

% Save information about the optimisation
info.Obj = Obj;
info.it = it;
