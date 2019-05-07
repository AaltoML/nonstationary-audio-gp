function test_suite =  test_fitSEGP_BF
  initTestSuite;

% test 
% function [z,info] = fitSEGP_BF(y,lenx,varx,varargin)

function test_inf

dispFigs=1;
randn('state',1)

lenx = 10*exp(randn);
varx = exp(randn);

T = 10000;
y = sampleGPSE(varx,lenx,T);
[z,info] = fitSEGP_BF(y,lenx,varx);

tol = 9;
yEst = basis2SEGP(z,lenx,varx,tol);

if dispFigs==1
  figure
  hold on
  plot(y,'-r');
  plot(yEst,'-k');
end

tol = 1e-1;
assertVectorsAlmostEqual(mean(abs(y-yEst)),0,'absolute',tol,0)

