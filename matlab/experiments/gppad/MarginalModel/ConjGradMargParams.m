function [varc,varx,mux,Info] = ConjGradMargParams(y,varc,varx,mux,Opts);

% function [varc,varx,mux,Info] = ConjGradMargParams(y,varc,varx,mux,Opts);  
% 
% Performs a conjugate gradients to estimate parameters from
% the marginal data

  theta = [log(varc);log(varx);mux];

  NumIts = Opts.NumItsMM;
  
  [theta,Obj,i] = minimize(theta,'GetObjNumericalInt2',NumIts,y,Opts);
  
  varc = exp(theta(1));
  varx = exp(theta(2));
  mux = theta(3);

  Info.Obj = Obj;