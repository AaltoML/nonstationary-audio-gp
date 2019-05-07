function test_suite = test_randnmf
  initTestSuite;

% Tests: 
%
% function [A,H] = randnmf(W,muinf,varinf,lam,T)
%

function test_check_gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 100000;
D = 4;
K = 3;

W = exp(randn(K,D));

% mode =b/(a+1) = h_{t-1} 
% variance = b^2/((a-1)^2(a-2)) = h_{t-1}^2 (a-1)^2 / ((a+1)^2 (a-2))
muinf = [1/2,3/4,1/3];
varinf = [1/8,1/16,1/4];
lam = [0,0.15,0.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST
[A,H] = randnmf(W,muinf,varinf,lam,T);

mnH = mean(H,1);
varH = var(H,1);
meanlagH = (mean(H(2:end,:).*H(1:end-1,:),1)-mnH.^2)./varH;

[muinf;mnH]
[varinf;varH]
[lam;meanlagH]

tol = 1e-1;
assertVectorsAlmostEqual(muinf,mnH,'absolute',tol,0)
assertVectorsAlmostEqual(varinf,varH,'absolute',tol,0)
assertVectorsAlmostEqual(lam,meanlagH,'absolute',tol,0)


