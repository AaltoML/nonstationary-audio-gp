function test_suite =  test_getGPSEVarSpec
  initTestSuite;

% test 
% function freqSq = getGPSEVarSpec(len);

function test_compareToSample

dispFigs=1;
randn('state',1)
var = 1.5;
len = 3.32;
tau = 20;
T = 40000;

freqSq1 = getGPSEVarSpec(len);

x = sampleGPSE(var,len,T);

freqSq2=getVarSpec(x);

[freqSq1,freqSq2]

tol = 1e-1;
assertVectorsAlmostEqual(freqSq1,freqSq2,'absolute',tol,0)
