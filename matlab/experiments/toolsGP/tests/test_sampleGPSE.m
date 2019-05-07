function test_suite =  test_sampleGPSE
  initTestSuite;

% test 
% function x = sampleGPSE(var,len,T)

function test_correctAutoCorrelation

dispFigs=1;
randn('state',1)
var = 1.5;
len = 5.32;
tau = 20;
T = 40000;

x = sampleGPSE(var,len,T);

autoCor = var*exp(-1/(2*len^2)*([0:tau-1]).^2);
 
autoCorEmp = zeros(tau,1);

for t=1:tau
  autoCorEmp(t) = mean(x(1:T-t+1).*x(t:T));
end

if dispFigs==1
  figure;
  hold on
  plot(autoCor,'-k','linewidth',2)
  plot(autoCorEmp,'-r','linewidth',1)
end

tol = 1e-1;
assertVectorsAlmostEqual(autoCor,autoCorEmp','absolute',tol,0)
