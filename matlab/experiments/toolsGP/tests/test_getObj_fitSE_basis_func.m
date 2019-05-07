function test_suite = test_getObj_fitSE_basis_func
  initTestSuite;

% Tests: 
%
% function [obj,dobj] = getObj_fitSE_basis_func(z,x,lenx,varx,tol)
%

function test_check_gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 100;
z = randn(T,1);
x = randn(T,1);
lenx = 20;
varx = 1.3;
tol  = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

delta = 1e-5;

d=checkgrad('getObj_fitSE_basis_func',z,delta,x,lenx,varx,tol);

tol = 1e-6;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)

