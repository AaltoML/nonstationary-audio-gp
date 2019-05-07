function test_suite = test_tnmf_inf
  initTestSuite;

% Tests: 
%
% function [H,info] = tnmf_inf(A,W,H,lenx,mux,varx,vary,varargin)
%

function test_simple_example_temp_TN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 200;
D = 5;
K = 3;


mux = [1,1,1];
varx = [3,3,3];
lenx = [10,20,30];
vary = zeros(T,D);

W = [1/3,0,1/3,0,1/3;
     0,1/3,0,2/3,0;
     1/3,1/3,1/3,0,0];

[A,H] = randtnmf(W,lenx,mux,varx,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

HInit = exp(randn(T,K))/1;

[HEst,info] = tnmf_inf(A,W,HInit,lenx,mux,varx,vary);

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

keyboard

function test_simple_example_temp_IG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 200;
D = 5;
K = 3;


muinf = [1,1,1];
varinf = [3,3,3];
lam = [0.2,0.9,0.95];
vary = zeros(T,D);

W = [1/3,0,1/3,0,1/3;
     0,1/3,0,2/3,0;
     1/3,1/3,1/3,0,0];

[A,H] = randnmf(W,muinf,varinf,lam,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

HInit = exp(randn(T,K))/1;

[HEst,info] = nmf_inf(A,W,HInit,muinf,varinf,lam,vary);

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

keyboard
