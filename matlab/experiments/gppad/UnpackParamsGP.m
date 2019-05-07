function [varx,len,mux,varc,vary] = UnpackParamsGP(Params); 

  % function [varx,len,mux,varc,vary] = UnpackParamsGP(Params); 
  %
  % Unpacks the parameter structure into individual
  % variables
  %  
  % INPUTS  
  % Params = Structure of parameters including
  %    vary = observation noise (optional - set to zero
  %           if absent)
  %    varx = scale parameter of stundent-t
  %    len = SE covariance function time-scale
  %    mux = mean parameter of stundent-t
  %
  % Outputs - as above, but as the component variables
    
  varx = Params.varx;
  len = Params.len;
  mux = Params.mux;
  varc = Params.varc;

  if isfield(Params,'vary')
    vary = Params.vary;
  else
    vary = 0;
  end    