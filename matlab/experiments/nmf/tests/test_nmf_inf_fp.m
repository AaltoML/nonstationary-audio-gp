function test_suite = test_nmf_inf_fp
  initTestSuite;

% Tests: 
%
%  function [W,H,info] = nmf_inf_fp(A,W,H,vary,varargin)
%

function test_simple_example

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 200;
D = 5;
K = 3;


muinf = [1,1,1];
varinf = [3,3,3];
lam = [0.2,0.9,0.95];
vary = rand([T,D])*0.02;

W = [1/3,0,1/3,0,1/3;
     0,1/3,0,2/3,0;
     1/3,1/3,1/3,0,0];

[A,H] = randnmf(W,muinf,varinf,lam,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

HInit = exp(randn(T,K))/1;

[HEst,info] = nmf_inf_fp(A,W,HInit,vary);

figure
plot(info.Obj)

figure
for k=1:K
  subplot(K,1,k)
  hold on
  plot(H(:,k),'-k')
  plot(HEst(:,k),'-r')
set(gca,'yscale','log')
end

assertTrue(max(diff(info.Obj))<0)

keyboard

