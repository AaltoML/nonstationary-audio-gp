function [xv,lmb] = GetEigSigDSig(d2px,d2py,Opts);
  
%function [xv,lmb] = GetEigDSigD(d2px,d2py,Opts);
%
% Finds the first K eigenvalues and eigenvectors of;
% A = Cov^{1/2} D Cov^{1/2} 
% For the GPPAD model  
%  
% INPUTS
% d2px = hessian of the log-prior [Tx,1]
% d2py = second derivative of the log-likelihoods for
%        each time  [T,1]
% Opts = structure of options containing
%        nev = Number of eigenvalues to find  
% 
% OUTPUTS
% xv = eigenvectors [Tx,Opts.nev]
% lmb = eigenvalues [Opts.nev,1] 
  
n = size(d2px,1);
d2pxhalf = d2px.^(1/2);

epert = 1e-9;
nev = Opts.nev;  % number of eigenvalues to find
jmax = 2*nev; % number of iterations
vstart = randn(n,1); % randn(n,1);
tresh = 1e-12;%1e-8 %tresh = 1e-10; % threshold for orthogonalising
tolconv = 1e-8;

complete = 0;
while complete==0
  [xv,lmb] = LanczosGPPAD(d2pxhalf,d2py,vstart,epert,nev,jmax,tresh,tolconv);
  if sum(isnan(xv))<1
    complete =1;
%    disp('Completed Lanczos')
  else
    disp('Error in Lanczos')
    vstart = randn(n,1); % randn(n,1);
  end
end

lmb = real(lmb);