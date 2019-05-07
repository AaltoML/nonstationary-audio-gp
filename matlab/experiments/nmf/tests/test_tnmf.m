function test_suite = test_tnmf
  initTestSuite;

% Tests: 
%
% function [W,H,info] = tnmf(A,W,H,lenx,mux,varx,vary,varargin)
%

function test_simple_example_truncated_normal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20000;
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

WInit = exp(randn(K,D))/1;
HInit = exp(randn(T,K))/1;


[WEst,HEst,info] = tnmf(A,WInit,HInit,lenx,mux,varx,vary);

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

figure
subplot(2,1,2)
imagesc(WEst)

subplot(2,1,1)
imagesc(W)

keyboard

