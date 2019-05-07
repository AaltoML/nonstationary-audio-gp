function [len,Info] = LearnLengthsGPPADFast(len0,y,Params,Opts);
  
% function [len,Info] = LearnLengthsGPPADFast(len0,y,Params,Opts);
%
% Learns the length-scale of the modulation by optimising
% Laplace's approximation to the marginal likelihood of
% the model.  
% 
  
% if size(len0)==1;
%   Len0 = [len0/1.5,len0,len0*1.5];
% elseif size(len0)==3
%   Len0 = len0;
% else
%   fprintf('\n\aError in LearnLengthsGPPAD.m - len0 incorrect size')
%   return
% end

%[Len,Obj,Info] = StepOutBisectDownSampleV2(Len0,y,Params,Opts);

  [Len,Obj,Info] = StepOutBisectDownSampleV2(len0,y,Params,Opts);

  % Output the best guess at the length-scale as well as
  % algorithmic information
  len = Len(2);

  Info.LenRng = Len;
  Info.LenObj = Obj;
