function test_suite = test_nmf_fp
  initTestSuite;

% Tests: 
%
% [W,H,info] = nmf_fp(A,W,H,muinf,varinf,lam,vary,varargin)
%

% COULDN'T GET THIS TO WORK RELIABLY
% function test_simple_example_truncated_normal_tnmf

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
% varTG = varinf.*(1-lam.^2);

% disp('fp')
% Opts.restarts = 0;
% Opts2.numIts = 100;
% [WEst,HEst,info] = nmf_fp(A,WInit,HInit,[],varTG,lam,vary,Opts);

% disp('conj-grad')
% Opts2.restarts = 0;
% Opts2.numIts = 100;
% [WEst2,HEst2,info2] = nmf(A,WInit,HInit,[],varTG,lam,vary,Opts);

% figure
% hold on
% plot(info.Obj,'-r')
% plot(info2.Obj,'-k')

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

% assertTrue(max(diff(info.Obj))<0)

function test_simple_example_normal_nmf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 20000;
D = 5;
K = 3;


muinf = [1,1,1];
varinf = [3,3,3];
lam = [0.2,0.9,0.99];
vary = rand([T,D])*0.03;

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

Opts.restarts = 10;
Opts.numIts = 1000;

disp('fp')
[WEst,HEst,info] = nmf_fp(A,WInit,HInit,vary,Opts);

disp('conj-grad')
[WEst2,HEst2,info2] = nmf(A,WInit,HInit,[],[],[],vary,Opts);


figure
hold on
plot(info.Obj,'-r')
plot(info2.Obj,'-k')

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
assertTrue(max(diff(info.Obj))<0)

