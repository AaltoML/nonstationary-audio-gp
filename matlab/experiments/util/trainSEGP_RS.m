function [lenx,varx,info] = trainSEGP_RS(x,varargin)

% function [lenx,varx,info] = trainSEGP_RS(x,varargin)
%
% Trains a regularly sampled squared exponential GP using conjugate
% gradient optimisation of the likelihood which is expressed in the
% spectral domain to make things fast.
%
% Seems to work nicely for lenx<100 but above this maximum-likelihood
% under-estimates the length-scale
%
% INPUT 
% x = signal
% optional
% opts = structure containing algorithmic options incl.
%   numIts = number of iterations   
%   progress_chunk = regularity of feedback
% lenx = length scale of GP
% varx = variance of GP
% vary = noise variance
%
% OUPUT
% lenx = length scale of GP
% varx = variance of GP
% info = structure containing information about the optimisation including
%     Obj = objective function
%     it = number of completed iteration steps per batch

if nargin>1 & isfield(varargin{1},'numIts')
  numIts = varargin{1}.numIts;
else
  numIts = 100;
end

if nargin>1 & isfield(varargin{1},'progress_chunk')
  progress_chunk = varargin{1}.progress_chunk;
else
  progress_chunk = 100;
end

specx = abs(fft(x)).^2;
T = length(x);

% initialisation if user doesn't provide starting values
if nargin<2
  varx = var(x);

  % initialisise x from the spectrum using the average square
  % frequency, where the avarage is taken w.r.t. the normalised spectrum
  
  if mod(T,2)==0
    om = linspace(0,1/2,T/2+1);
    om = [om,om(end-1:-1:2)]'; 
  else
    om = linspace(0,1/2,T/2+1);
    om = [om,om(end:-1:2)]';
  end
  
  ave_om_sq = (om'.^2*specx)/sum(specx);
  lenx = 1/(2*pi*sqrt(ave_om_sq));

  vary = 0.1;

else
  % if user has provided initial values then use them
  lenx = varargin{2};
  varx = varargin{3};
  vary = varargin{4};
end

% if nargin>1 & isfield(varargin{1},'opt_y')
% else
%   optsInit.opt_y = 1;
%   [lenx,varx,vary,infoInit] = trainSEGP_RS(x,optsInit,lenx,varx,vary);
% end

L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];

params = [log(lenx)];

for l=1:L
  
  % save old parameters
  lenxOld = lenx;
  varxOld = varx;
  varyOld = vary; 
    
  % this seemed to improve things due to overfitting -- but it's a
  % little ad hoc so I've commented it out
  % if nargin<2
  %   if lenx<100
  %     vary = 0.1;
  %   elseif lenx<300
  %     vary = 0.3;
  %   else 
  %     vary = 1;
  %   end
  % end

  % run conjugate gradient update
  tic;

  [params, ObjCur, itCur] = minimize(params,'getObjSEGP', ...
				     numIts(l),specx,vary,varx);
  lenx = exp(params(1));

  % GRID SEARCH FOR COMPARISON
  % M = 100;
  % obj = zeros(M,1);
  % lens = logspace(log10(10),log10(1000),M);
  % for m=1:M
  %   [obj(m),dobj] = getObjSEGP(log(lens(m)),specx,vary,varx);
  % end
  
  % [val,loc] = min(obj);
  % [lenx,lens(loc)]

  
  % ALSO IMPLEMENTED OPTIMISATION OVER THE OTHER PARAMETERS, BUT
  % THIS LED TO WORSE RESULTS
  % if nargin>1 & isfield(varargin{1},'opt_y')
%     % optimising over lenx only
%  %   disp(['trainSEGP_RS: optimising over lenx and vary only'])
  
%% params = [log(lenx);log(varx);log(vary)];

 
%   else
%     % optimising over all parameters
% %    disp(['trainSEGP_RS: optimising over all parameters'])

%     [params, ObjCur, itCur] = minimize(params,'getObjSEGP', ...
% 				      numIts(l),specx);
%     lenx = exp(params(1));
%     varx = exp(params(2));
%     vary = exp(params(3));

%   end
  timCur = toc;

  % Store objective and iteration information
  Obj = [Obj;ObjCur];
  it = [it;itCur];
  tim = [tim;timCur];
  dvary = vary-varyOld;
  dlenx = lenx-lenxOld;
  dvarx = varx-varxOld;
  
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
  str6 = ['vary ',num2str(vary)];
  str7 = ['varx ',num2str(varx)];
  str8 = ['lenx ',num2str(lenx)];

  str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6,str_space,str7,str_space,str8])

end

fprintf('\n')

% Save information about the optimisation
info.Obj = Obj;
info.it = it;

if lenx>100
  disp(['trainSEGP_RS  WARNING lenx>100, might be over-fitting and underestimating lenx'])
end