function test_suite = test_getObj_nmf_temp
  initTestSuite;

% Tests: 
%
% function [obj,dobj] = getObj_nmf_temp(logHW,A,vary,muinf,varinf,lam)
%

function test_check_gradient_temporal_truncated_gaussian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
D = 3;
K = 2;

logHW = randn(T*K+K*D,1);
A = exp(randn(T,D));

varinf = [1/16,1/32];
lam = [0.1,0.3];

%muinf = [1/4,1/8];
%varinf = [1/16,1/32];
%lam = [0.1,0.2];
vary = rand(T,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObj_nmf_temp',logHW,delta,A,vary,varinf,lam);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


function test_check_gradient_temporal_inverse_gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
D = 3;
K = 2;

logHW = randn(T*K+K*D,1);
A = exp(randn(T,D));

muinf = [1/4,1/4];
varinf = [1/16,1/16];
lam = [0.1,0.1];

%muinf = [1/4,1/8];
%varinf = [1/16,1/32];
%lam = [0.1,0.2];
vary = rand(T,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObj_nmf_temp',logHW,delta,A,vary,muinf,varinf,lam);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)



function test_check_gradient_non_temporal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
D = 3;
K = 2;

logHW = randn(T*K+K*D,1);
A = exp(randn(T,D));

vary = rand(T,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObj_nmf_temp',logHW,delta,A,vary);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)
