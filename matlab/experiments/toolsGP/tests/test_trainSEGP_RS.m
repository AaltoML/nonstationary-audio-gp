function test_suite =  test_trainSEGP_RS
  initTestSuite;

% test 
% function [lenx,varx,vary,info] = trainSEGP_RS(x,varargin)

function test_inf_true

dispFigs=1;
randn('state',1)

L = 20;
Ntrials = 10;
lenxs = logspace(log10(1),log10(500),L);
lenxsEst = zeros(Ntrials,L);
varxsEst = zeros(Ntrials,L);
varysEst = zeros(Ntrials,L);
T = 10000;


for l=1:L
  lenx = lenxs(l);
  for n=1:Ntrials
    varx = exp(randn);
 
    x = sampleGPSE(varx,lenx,T);
  
    [lenxsEst(n,l),varxsEst(n,l),info] = trainSEGP_RS(x);
    
   %  x = x -mean(x);
   %  spec = fft(x);
    
   %  if mod(T,2)==0
   %    om = linspace(0,1/2,T/2+1);
   %    om = [om,om(end-1:-1:2)]'; 
   %  else
   %    om = linspace(0,1/2,T/2+1);
   %    om = [om,om(end:-1:2)]';
   %  end
    
   % ave_om_sq = (om'.^2*abs(spec).^2)/sum(abs(spec).^2);
   % lenxsEst(n,l) = 1./(2*pi*sqrt(ave_om_sq));
    
  end
end

figure
hold on
plot(lenxs,lenxs,'-r');
for l=1:L
  plot(lenxs(l),lenxsEst(:,l),'.k');
end
plot(lenxs,mean(lenxsEst,1),'-b');
plot(lenxs,mean(lenxsEst,1)+sqrt(var(lenxsEst,1)),'-b');
plot(lenxs,mean(lenxsEst,1)-sqrt(var(lenxsEst,1)),'-b');

%tol = 1e-1;
%assertVectorsAlmostEqual(autoCor,autoCorEmp','absolute',tol,0)

dl = abs(ones(Ntrials,1)*lenxs-lenxsEst)./(ones(Ntrials,1)*lenxs);
disp(['average percentage error ',num2str(100*mean(dl(:))),'%']);
keyboard