function test_suite = test_getObj_nmf_log_norm
  initTestSuite;

% Tests: 
%
% function [obj,dobj] = getObj_nmf_log_norm(xlogW,A,vary,lenx,mux,varx,tol)

function test_check_gradient_temporal_truncated_gaussian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
D = 3;
K = 2;

xlogW = randn(T*K+K*D,1);
A = exp(randn(T,D));

mux = [1,2];
lenx = [3,45];
varx = [0.1,0.3];

%muinf = [1/4,1/8];
%varinf = [1/16,1/32];
%lam = [0.1,0.2];
vary = rand(T,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST
tol = 9;

delta = 1e-5;

d=checkgrad('getObj_nmf_log_norm',xlogW,delta,A,vary,lenx,mux,varx,tol);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)

