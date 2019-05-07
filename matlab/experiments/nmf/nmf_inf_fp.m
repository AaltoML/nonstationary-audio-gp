function [H,info] = nmf_inf_fp(A,W,H,vary,varargin)

% function [W,info] = nmf_inf_fp(A,W,H,vary,varargin)
%
% Performs fixed point inference updates for the standard NMF model
% w/o any temporal constraints
% Ws must be normalised for this to work properly
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
  numIts = 100;
end

[K,D] = size(W);
[T,K] = size(H);

Obj = []; tim = [];

if sum(W,2)~=ones(K,1)
  disp('Error in nmf_inf_fp: W unnormalised -- normalising')
  W = diag(1./sum(W,2))*W;
end

for l=1:numIts
  
  % save old parameters
  Hold = H;
  
  % run fixed point update
  tic;

  AHat = H*W+vary;
  H = ((A.*AHat.^(-2))*W')./(AHat.^(-1)*W').*H;
  
  logHW = [log(H(:));log(W(:))];
  [ObjCur,dObj] = getObj_nmf_temp(logHW,A,vary);
  Obj = [Obj,ObjCur];


  timCur = toc;
  tim = [tim;timCur];
  
  if mod(l,ceil(numIts/10))==0&l>1
    Ahat = H*W;
    snrChan = 10*log10(mean(A.^2,1))-10*log10(mean((A-Ahat).^2));
  
    % % Store objective and iteration information
    
    dH = sqrt(sum((H(:)-Hold(:)).^2)/sum(H(:).^2));
    
    % % Display some information to the user
    str1 = ['Progress ',num2str(l),'/',num2str(numIts),];
    str2 = ['Obj ',num2str(ObjCur(end))];
    str3 = ['dObj ',num2str(diff(Obj(end-1:end)))];
    str4 = ['time ',num2str(ceil(timCur/6)/10),'mins'];
    str5 = ['total time ',num2str(ceil(sum(tim)/6)/10),'mins'];
    str6 = ['dH ',num2str(round(dH*1000)/10),'%%'];
    str7 = ['A snr ',num2str(round(mean(snrChan)*1000)/1000)];
  
    str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6,str_space,str7,str_space])

  end


end

fprintf('\n')

% Save information about the optimisation
info.Obj = Obj;
%info.it = it;
