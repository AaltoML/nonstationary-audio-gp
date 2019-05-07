function Params = PackParamsGP(varx,len,mux,varc,vary); 

  % function Params = PackParamsGP(varx,mux,varc,vary); 
  %
  % Packs the individual variables into a parameter structure 
  %  
  % INPUTS  
  %    vary = observation noise
  %    varx = scale parameter of SE
  %    varc = scale parameter of carrier
  %    mux = mean parameter of stundent-t
  %    len = length-scale SE covariance
  %
  % Outputs - as above, but as fields of the structure Params
    
  Params.varx = varx;
  Params.len = len;
  Params.mux = mux;
  Params.vary = vary;
  Params.varc = varc;

    