function [z,info] = fitSEGP_BF(y,lenx,varx,varargin)

% function [z,info] = fitSEGP_BF(y,lenx,varx,varargin)
%
% Fits a squared exponential GP to the signal using conjugate
% gradient optimisation of the likelihood. Uses the basis function
% representation of the SE GP (see Mackay 2003) and returns the
% coefficients of the basis functions.
%
%
% INPUT 
% y = signal
% lenx = length scale of GP
% varx = variance of GP
%
% optional inputs
% opts = structure containing algorithmic options incl.
%   numIts = number of iterations   
%   progress_chunk = regularity of feedback
% tol = tolerence which sets the extent of the temporal basis functions
%
% OUPUT
% lenx = length scale of GP
% varx = variance of GP
% info = structure containing information about the optimisation including
%     Obj = objective function
%     it = number of completed iteration steps per batch

if nargin>3 & isfield(varargin{1},'numIts')
  numIts = varargin{1}.numIts;
else
  numIts = 100;
end

if nargin>3 & isfield(varargin{1},'progress_chunk')
  progress_chunk = varargin{1}.progress_chunk;
else
  progress_chunk = 50;
end

if nargin>3 & isfield(varargin{1},'progress_chunk')
  tol = varargin{1}.tol;
else
  tol = 9;
end

T = length(y);

% initialisation if user doesn't provide starting values
z = randn(T,1)/1000;

L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];

for l=1:L
  
  % save old parameters
  zOld = z;

  % run conjugate gradient update
  tic;

  [z, ObjCur, itCur] = minimize(z,'getObj_fitSE_basis_func', ...
				     numIts(l),y,lenx,varx,tol);
  
  timCur = toc;

  % Store objective and iteration information
  Obj = [Obj;ObjCur];
  it = [it;itCur];
  tim = [tim;timCur];
  dz = mean((z-zOld).^2);
  
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
  str6 = ['dz ',num2str(dz)];

  str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6])

end

fprintf('\n')

% Save information about the optimisation
info.Obj = Obj;
info.it = it;

