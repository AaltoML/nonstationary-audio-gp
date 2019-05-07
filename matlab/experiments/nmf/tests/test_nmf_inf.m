function test_suite = test_nmf_inf
  initTestSuite;

% Tests: 
%
% [H,info] = nmf_inf(A,W,H,muinf,varinf,lam,vary,varargin)
%

function test_simple_example_temp_TN

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

[HEst,info] = nmf_inf(A,W,HInit,[],varinf,lam,vary);

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
