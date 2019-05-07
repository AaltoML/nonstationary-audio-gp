function [W,H,info] = nmf(A,W,H,muinf,varinf,lam,vary,varargin)

% function [W,H,info] = nmf(A,W,H,muinf,varinf,lam,vary,varargin)
%
% Performs NMF on the data matrix A. That is, it approximates A
% via A \approx H*W. The NMF variant is trained using conjugate gradients
% and the temporal and spectral basis function are learned in a joint
% optimisation.
%
% The cost function used in the likelihood is:
% \sum_{d,t} A_{d,t}./Ahat_{t,d} + \sum_{d,t} \log(Ahat_{t,d}) 
%
% If you want to run normal NMF with no temporal priors set 
% muinf = [], varinf = [], lam = []
%
% The NMF variant used has (optional) temporal constraints (K = number
% of basis functions can be under/over/complete).
%
% If you want to run tNMF with truncated normal temporal priors set 
% muinf = [], varinf = x, lam = y
%
% If you want to run tNMF with inverse gamma temporal priors set 
% muinf = x, varinf = y, lam = z, i.e. no empty matrices
% 
% Inverse Gamma Markov chain priors on the temporal basis functions
% are as follows:
% h_{k,1} ~ IG(alp1,bet1)
% h_{k,t}|h_{k,t-1} ~ IG(alp_{k},bet_{k,t})
% alp_k = 2+lam_k^2+muinf_k^2/varinf_k
% bet_{t,k} = (alp_{t,k}-1)(lam_k*(h_{k,t-1}-muinf_k)+muinf_k)
%
% The steady state mean of the IG chain is muinf(k), the variance
% varinf(k), and the correlation coefficient between successive
% variables is lam(k).
%
% The first variable in the chain is drawn from an inverse Gamma
% distribution with mean and variance equal to the steady state.
%
%
%
% INPUTS
% Model parameters:
%   A = data [T,D]
%   W = spectral weight vector [K,D]
%   H = temporal weight vector [K,T]
%   muinf = steady state mean of the IGMC priors [1,K]
%   varinf = steady state variance of the IGMC priors [1,K]
%   lam = steady state correlation-coefficient for successive
%       variables in the IGMC priors [1,K]
%   vary = observation noise [T,D] (set to zero if signal is noiseless)
%
% Optional algorithmic options:
% numIts = number of iterations per optimisation batch
% progress_chunk = number of optimisation batches 
%
% OUTPUTS
%   W = learned spectral weight vector [K,D]
%   H = learned temporal weight vector [K,T]
%   info = structure containing information about the optimisation including
%     Obj = objective function
%     it = number of completed iteration steps per batch

if nargin>7 & isfield(varargin{1},'numIts')
  numIts = varargin{1}.numIts;
else
  numIts = 1000;
end

if nargin>7 & isfield(varargin{1},'progress_chunk')
  progress_chunk = varargin{1}.progress_chunk;
else
  progress_chunk = 1000;
end

[K,D] = size(W);
[T,K] = size(H);

if nargin>7 & isfield(varargin{1},'restarts')
  % initialises with restarts if desired (first initialisation is
  % the one provided by the user)
  restarts = varargin{1}.restarts;
  Wtest = W;
  Htest = H;
  ObjBest = inf;
  disp(['---- using ',num2str(restarts),' restarts ----'])
  Opts.numIts = 1000;
  
  for r=1:restarts
%    this should be updated to use the fixed point algorithm for initialisation
     [Htest,info] = nmf_inf(A,Wtest,Htest,muinf,varinf,lam,vary,Opts);

   if info.Obj(end)<ObjBest
      W = Wtest;
      H = Htest;
      ObjBest = info.Obj(end);
    end
    ks = ceil(T*rand(K,1));
    Wtest = A(ks,:);%exp(randn(K,D))/1;
    Htest = exp(randn(T,K))/1;
  end
 
end

L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];

logHW = [log(H(:));log(W(:))];

for l=1:L
  
  % save old parameters
  Hold = H;
  Wold = W;
  
  % run conjugate gradient update
  tic;
  if isempty(varinf)
    % normal NMF
    [logHW, ObjCur, itCur] = minimize(logHW,'getObj_nmf_temp', ...
				      numIts(l),A,vary);
  elseif isempty(muinf)
    % temporal NMF truncated Gaussian
    [logHW, ObjCur, itCur] = minimize(logHW,'getObj_nmf_temp', ...
				      numIts(l),A,vary,varinf,lam);
  else
    % temporal NMF inverse Gamma
    [logHW, ObjCur, itCur] = minimize(logHW,'getObj_nmf_temp', ...
				      numIts(l),A, vary,muinf,varinf,lam);
  end
  timCur = toc;

  % Pull out H and W from logHW
  logH = reshape(logHW(1:K*T),[T,K]);
  H = exp(logH);
  W = reshape(exp(logHW(K*T+1:end)),[K,D]);
  W = diag(1./sum(W,2))*W;
  Ahat = H*W;
  snrChan = 10*log10(mean(A.^2,1))-10*log10(mean((A-Ahat).^2));
  
  % Store objective and iteration information
  Obj = [Obj;ObjCur];
  it = [it;itCur];
  tim = [tim;timCur];
  dW = sqrt(sum((W(:)-Wold(:)).^2)/sum(W(:).^2));
  dH = sqrt(sum((H(:)-Hold(:)).^2)/sum(H(:).^2));
  
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
  str7 = ['dW ',num2str(round(dW*1000)/10),'%%'];
  str8 = ['A snr ',num2str(round(mean(snrChan)*1000)/1000)];


  str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6,str_space,str7,str_space,str8,str_space])

end

fprintf('\n')

% Save information about the optimisation
info.Obj = Obj;
info.it = it;
