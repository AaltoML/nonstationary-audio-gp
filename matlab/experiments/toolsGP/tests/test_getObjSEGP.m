function test_suite = test_getObjSEGP
  initTestSuite;

% Tests: 
%
% [obj,dob] = getObjSEGP(param,specy,varargin)
%

function test_check_gradient_full

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
y = randn(T,1);
specy = abs(y).^2;

param = randn(3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObjSEGP',param,delta,specy);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


function test_check_gradient_partial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
y = randn(T,1);
specy = abs(y).^2;

param = randn(2,1);
vary = exp(randn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObjSEGP',param,delta,specy,vary);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)



function test_check_gradient_partial_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
y = randn(T,1);
specy = abs(y).^2;

param = randn(1,1);
vary = exp(randn);
varx = exp(randn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObjSEGP',param,delta,specy,vary,varx);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


