function [W,H,info] = nmf_fp(A,W,H,vary,varargin)

% function [W,H,info] = nmf_fp(A,W,H,vary,varargin)
%
% Performs fixed point updates for the standard NMF model w/o any
% temporal constraints
%
% INPUTS
% Model parameters:
%   A = data [T,D]
%   W = spectral weight vector [K,D]
%   H = temporal weight vector [K,T]
%   vary = observation noise [T,D] (set to zero if signal is noiseless)
%
% Optional algorithmic options:
% numIts = number of iterations per optimisation batch
%
% OUTPUTS
%   W = learned spectral weight vector [K,D]
%   H = learned temporal weight vector [K,T]
%   info = structure containing information about the optimisation including
%     Obj = objective function
%     it = number of completed iteration steps per batch

if nargin>4 & isfield(varargin{1},'numIts')
  numIts = varargin{1}.numIts;
else
  numIts = 1000;
end

[K,D] = size(W);
[T,K] = size(H);

if nargin>4 & isfield(varargin{1},'restarts')
  % initialises with restarts if desired (first initialisation is
  % the one provided by the user)
  restarts = varargin{1}.restarts;
  Wtest = W;
  Htest = H;
  ObjBest = inf;
  disp(['---- using ',num2str(restarts),' restarts ----'])
  Opts.numIts = 10;

  for r=1:restarts
    Wtest = diag(1./sum(Wtest,2))*Wtest;
     [Htest,info] = nmf_inf_fp(A,Wtest,Htest,vary,Opts);

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


Obj = []; tim = [];


W = diag(1./sum(W,2))*W;
    
for l=1:numIts
  
  % save old parameters
  Hold = H;
  Wold = W;
  
  % run fixed point update
  tic;

  AHat = H*W+vary;
  H = ((A.*AHat.^(-2))*W')./(AHat.^(-1)*W').*H;
  
  logHW = [log(H(:));log(W(:))];
  [ObjCur,dObj] = getObj_nmf_temp(logHW,A,vary);
  Obj = [Obj,ObjCur];

  AHat = H*W;
  W = (H'*(A.*AHat.^(-2)))./(H'*AHat.^(-1)).*W;
  W = diag(1./sum(W,2))*W;
  
  logHW = [log(H(:));log(W(:))];
  [ObjCur,dObj] = getObj_nmf_temp(logHW,A,vary);
  Obj = [Obj,ObjCur];

% I derived the following updates for the temporal variant of the
% model with truncated Gaussian priors, but it did not behave very well:
%     dHpos = zeros(T,K);
%     dHneg = zeros(T,K);
    
%     dHpos(1,:) = (lam.^2./varinf).*H(1,:);
%     dHneg(1,:) = (lam./varinf).*H(2,:);
    
%     dHpos(2:T-1,:) = ones(T-2,1)*((1+lam.^2)./varinf).*H(2:T-1,:);
%     dHneg(2:T-1,:) = ones(T-2,1)*(lam./varinf).*(H(3:T,:)+H(1:T-2,:));
		
%     dHpos(T,:) = 1./varinf.*H(T,:);
%     dHneg(T,:) = (lam./varinf).*H(T-1,:);
    
%     AHat = H*W;
%     H = ((A.*AHat.^(-2))*W' + dHneg)./(AHat.^(-1)*W'+dHpos).*H;
%     %keyboard
%     logHW = [log(H(:));log(W(:))];
%     [ObjCur,dObj] = getObj_nmf_temp(logHW,A,vary,varinf,lam);
%     Obj = [Obj,ObjCur];

%     AHat = H*W;
%     W = (H'*(A.*AHat.^(-2)))./(H'*AHat.^(-1)).*W;
%     W = diag(1./sum(W,2))*W;
  
%     logHW = [log(H(:));log(W(:))];
%     [ObjCur,dObj] = getObj_nmf_temp(logHW,A,vary,varinf,lam);
%     Obj = [Obj,ObjCur];

% %    keyboard
    
%   else
%     % temporal NMF inverse Gamma
%     [logHW, ObjCur, itCur] = minimize(logHW,'getObj_nmf_temp', ...
% 				      numIts(l),A, vary,muinf,varinf,lam);
%   end
  timCur = toc;
  tim = [tim;timCur];
  
  if mod(l,100)==0
    Ahat = H*W;
    snrChan = 10*log10(mean(A.^2,1))-10*log10(mean((A-Ahat).^2));
  
    % % Store objective and iteration information
    % Obj = [Obj;ObjCur];
    % it = [it;itCur];
    
    dW = sqrt(sum((W(:)-Wold(:)).^2)/sum(W(:).^2));
    dH = sqrt(sum((H(:)-Hold(:)).^2)/sum(H(:).^2));
    
    % % Display some information to the user
    str1 = ['Progress ',num2str(l),'/',num2str(numIts),];
    str2 = ['Obj ',num2str(ObjCur(end))];
    str3 = ['dObj ',num2str(diff(Obj(end-1:end)))];
    str4 = ['time ',num2str(ceil(timCur/6)/10),'mins'];
    str5 = ['total time ',num2str(ceil(sum(tim)/6)/10),'mins'];
    str6 = ['dH ',num2str(round(dH*1000)/10),'%%'];
    str7 = ['dW ',num2str(round(dW*1000)/10),'%%'];
    str8 = ['A snr ',num2str(round(mean(snrChan)*1000)/1000)];
  
    str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6,str_space,str7,str_space,str8,str_space])

  end


end

fprintf('\n')

% Save information about the optimisation
info.Obj = Obj;
%info.it = it;
