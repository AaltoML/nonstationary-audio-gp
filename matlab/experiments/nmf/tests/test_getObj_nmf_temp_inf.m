function test_suite = test_getObj_nmf_temp_inf
  initTestSuite;

% Tests: 
%
% function [obj,dobj] = getObj_nmf_temp_inf(logH,A,W,vary,varargin)
%

function test_check_gradient_temporal_TG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
D = 3;
K = 2;

logH = randn(T*K,1);
A = exp(randn(T,D));
W = exp(randn(K,D));

varinf = [1/16,1/16];
lam = [0.1,0.1];
vary = rand(T,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObj_nmf_temp_inf',logH,delta,A,W,vary,varinf,lam);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)

function test_check_gradient_temporal_IG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
D = 3;
K = 2;

logH = randn(T*K,1);
A = exp(randn(T,D));
W = exp(randn(K,D));

muinf = [1/4,1/4];
varinf = [1/16,1/16];
lam = [0.1,0.1];
vary = rand(T,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObj_nmf_temp_inf',logH,delta,A,W,vary,muinf,varinf,lam);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


function test_check_gradient_non_temporal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20;
D = 3;
K = 2;

logH = randn(T*K,1);
A = exp(randn(T,D));
W = exp(randn(K,D));

vary = rand(T,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObj_nmf_temp_inf',logH,delta,A,W,vary);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)
