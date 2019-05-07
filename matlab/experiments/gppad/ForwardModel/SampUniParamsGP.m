function Params = SampUniParamsGP(LenRng,MuRng,VarxRng,VarcRng,VaryRng)
  
% function Params = SampUniParamsGP(LenRng,MuRng,VarxRng,VarcRng,VaryRng)
%
% Samples parameters from uniform distributions between
% the limits given in the inputs
  
Params.len = LenRng(1)+rand*diff(LenRng);  
Params.mux = MuRng(1)+rand*diff(MuRng);  
Params.varx = VarxRng(1)+rand*diff(VarxRng);  
Params.varc = VarcRng(1)+rand*diff(VarcRng);  
Params.vary = VaryRng(1)+rand*diff(VaryRng);  
