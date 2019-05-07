

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

T = 50000;
D = 100;
K = 100;

logHW = randn(T*K+K*D,1);
A = exp(randn(T,D));
alp = 10;
alp1 = 10;
bet1 = 10;

profile on;

[obj,dobj] = getObj_nmf_temp(logHW,A,alp,alp1,bet1);

profile off;
profile viewer;


