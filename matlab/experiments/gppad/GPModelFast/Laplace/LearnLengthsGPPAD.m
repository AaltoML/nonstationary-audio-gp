function [len,Info] = LearnLengthsGPPAD(len0,y,Params,Opts);
  
% function [len,Info] = LearnLengthsGPPAD(len0,y,Params,Opts);
%
% Learns the length-scale of the modulation by optimising
% Laplace's approximation to the marginal likelihood of
% the model.  
% 
% The optimisation occurs in two stages. The first stage,
% called the step-out algorithm, finds three length scales,
% l=[l1,l2,l3] where l1<l2:l3 such that Obj(l2)>Obj(l1) and
% Obj(l2) > Obj(l3) i.e. it returns the range in which
% the optimum lives. The second stage is called the
% bisection algorithm because it hones down this range by
% bisecting it repeatedly.
%
% Very roughly speaking, the step out algorithm gets within
% a factor of two of the length-scale i.e. it gets to within
% a length-scale away from the truth. The bisection
% algorithm then reduces the error by
% (1/2)^(NumItsBisect/2). So, for NumItsBisect = 10, this
% amounts to an error of about 3%.  
  
% With regard to the computational cost, each iteration of
% the bisection alogrithm involves the computation of a
% single objective. Each iteration of the the step out
% algorithm involves 3 objective evaluations.
  
if size(len0)==1;
  Len0 = [len0/2,len0,len0*2];
elseif size(len0)==3
  Len0 = len0;
else
  fprintf('\n\aError in LearnLengthsGPPAD.m - len0 incorrect size')
  return
end

% Run the stepping out algorithm to estimate a range in
% which the optimal length-scale lies
[Len1,Obj1,Info1] = StepOutLengthFast(Len0,y,Params,Opts);

if Info1.DS==1
  Opts.TCh = Info1.TCh;
  Opts.TxCh = Info1.TxCh;
  Len2 = Len1;
  Obj2 = Obj1;
  Info2 = 'Step Out was not required';
else
  [Len2,Obj2,Info2] = StepOutLength(Len1,y,Params,Opts);
  Opts.TCh = Info2.TCh;
  Opts.TxCh = Info2.TxCh;
end

% Run the bisection algorithm to reduce the range in
% which the optimal length-scale lies
[Len3,Obj3,Info3] = BisectionLaplace(Len2,Obj2,y,Params,Opts);

% Output the best guess at the length-scale as well as
% algorithmic information
len = Len3(2);

Info.LenRng = Len3;
Info.LenObj = Obj3;

Info.InfoStepOutFast = Info1;
Info.InfoStepOut = Info2;
Info.InfoBisect = Info3;
