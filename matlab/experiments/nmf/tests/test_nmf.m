function test_suite = test_nmf
  initTestSuite;

% Tests: 
%
% [W,H,info] = nmf(A,W,H,muinf,varinf,lam,vary,varargin)
%

function test_simple_example_truncated_normal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20000;
D = 5;
K = 3;


muinf = [1,1,1];
varinf = [3,3,3];
lam = [0.2,0.9,0.99];
vary = zeros(T,D);

W = [1/3,0,1/3,0,1/3;
     0,1/3,0,2/3,0;
     1/3,1/3,1/3,0,0];
%W = exp(randn(K,D))
[A,H] = randnmf(W,muinf,varinf,lam,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST

WInit = exp(randn(K,D))/1;
HInit = exp(randn(T,K))/1;

%WInit = W;
%HInit = H;
varTG = varinf.*(1-lam.^2);

Opts.restarts = 5;
[WEst,HEst,info] = nmf(A,WInit,HInit,[],varTG,lam,vary,Opts);

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

% function test_simple_example_temp_restarts

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Forward model Settings

% T = 20000;
% D = 5;
% K = 3;


% muinf = [1,1,1];
% varinf = [3,3,3];
% lam = [0.2,0.9,0.99];
% vary = zeros(T,D);

% W = [1/3,0,1/3,0,1/3;
%      0,1/3,0,2/3,0;
%      1/3,1/3,1/3,0,0];
% %W = exp(randn(K,D))
% [A,H] = randnmf(W,muinf,varinf,lam,T);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TEST

% WInit = exp(randn(K,D))/1;
% HInit = exp(randn(T,K))/1;

% %WInit = W;
% %HInit = H;

% Opts.restarts = 5;
% [WEst,HEst,info] = nmf(A,WInit,HInit,muinf,varinf,lam,vary,Opts);

% figure
% plot(info.Obj)

% figure
% for k=1:K
%   subplot(K,1,k)
%   hold on
%   plot(H(:,k),'-k')
%   plot(HEst(:,k),'-r')
% set(gca,'yscale','log')
% end

% figure
% subplot(2,1,2)
% imagesc(WEst)

% subplot(2,1,1)
% imagesc(W)

% keyboard
%tol = 1e-6;
%assertVectorsAlmostEqual(d,0,'absolute',tol,0)

% function test_simple_example_temp

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Forward model Settings

% T = 20000;
% D = 5;
% K = 3;


% muinf = [4,4,4];
% varinf = [3,3,3];
% lam = [0.2,0.9,0.99];
% vary = zeros(T,D);

% W = [1/3,0,1/3,0,1/3;
%      0,1/3,0,2/3,0;
%      1/3,1/3,1/3,0,0];
% %W = exp(randn(K,D))
% [A,H] = randnmf(W,muinf,varinf,lam,T);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TEST

% WInit = exp(randn(K,D))/1;
% HInit = exp(randn(T,K))/1;

% %WInit = W;
% %HInit = H;

% [WEst,HEst,info] = nmf(A,WInit,HInit,muinf,varinf,lam,vary);

% figure
% plot(info.Obj)

% figure
% for k=1:K
%   subplot(K,1,k)
%   hold on
%   plot(H(:,k),'-k')
%   plot(HEst(:,k),'-r')
% set(gca,'yscale','log')
% end

% figure
% subplot(2,1,2)
% imagesc(WEst)

% subplot(2,1,1)
% imagesc(W)

% keyboard
% %tol = 1e-6;
% %assertVectorsAlmostEqual(d,0,'absolute',tol,0)


% function test_simple_example

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Forward model Settings

% T = 20000;
% D = 5;
% K = 3;


% muinf = [4,4,4];
% varinf = [3,3,3];
% lam = [0.2,0.9,0.99];
% vary = zeros(T,D);

% W = [1/3,0,1/3,0,1/3;
%      0,1/3,0,2/3,0;
%      1/3,1/3,1/3,0,0];
% %W = exp(randn(K,D))
% [A,H] = randnmf(W,muinf,varinf,lam,T);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TEST

% WInit = exp(randn(K,D))/1;
% HInit = exp(randn(T,K))/1;

% %WInit = W;
% %HInit = H;

% [WEst,HEst,info] = nmf(A,WInit,HInit,[],[],[],vary);

% figure
% plot(info.Obj)

% figure
% for k=1:K
%   subplot(K,1,k)
%   hold on
%   plot(H(:,k),'-k')
%   plot(HEst(:,k),'-r')
% set(gca,'yscale','log')
% end

% figure
% subplot(2,1,2)
% imagesc(WEst)

% subplot(2,1,1)
% imagesc(W)

% keyboard
% %tol = 1e-6;
% %assertVectorsAlmostEqual(d,0,'absolute',tol,0)


