function [y,a,x,Params,Info] = FMGPFastSampParams(varcRng, ...
                                             varxRng, ...
                                             muRng, varyRng, ...
                                             lenRng,Dims,MinModDepth)

% function [y,a,x,Params,Info] = FMGPFastSampParams(varcRng, ...
%                                              varxRng, ...
%                                              muRng,
%                                              varyRng, ...
%                                              lenRng,Dims,MinModDepth)
% 
% Samples parameters and data from the forward model. The
% data is sampled by first sampling the parameters uniformly
% from the ranges given as inputs, and then sampling the
% modulator and signal. Rejection sampling is used to reject
% any signal which has too small a modulation depth, where
% the threshold is determined by MinModDepth. 

  
numrejects = 0;  
r=1;
  
while r==1
  
  varc = varcRng(1)+rand*diff(varcRng);
  varx = varxRng(1)+rand*diff(varxRng);
  mux = muRng(1)+rand*diff(muRng);
  len = lenRng(1)+rand*diff(lenRng);
  vary = varyRng(1)+rand*diff(varyRng);
  
  Params = PackParamsGP(varx,len,mux,varc,vary); 
  
  [y, a, x] = FMGPFast(Params,Dims);
  
  if sqrt(var(a))/mean(a)>=MinModDepth
    r=0;
  else
    numrejects = numrejects+1;
  end
  
end

Info.NumRejects = numrejects;