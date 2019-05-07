function [A,H] = randnmf(W,muinf,varinf,lam,T)

% [A,H] = randnmf(W,muinf,varinf,lam,T)
%
% Draws from the tNMF generative model:
%
% Uses Inverse Gamma Markov chain priors on the temporal basis functions:
% h_{k,1} ~ IG(alp1,bet1)
% h_{k,t}|h_{k,t-1} ~ IG(alp_{k},bet_{k,t})
% alp_k = 2+lam_k^2+muinf_k^2/varinf_k
% bet_{t,k} = (alp_{t,k}-1)(lam_k*(h_{k,t-1}-muinf_k)+muinf_k)
%
% The steady state mean of the IG chain is muinf(k), the variance
% varinf(k), and the correlation coefficient between successive
% variables is lam(k).
%
% The first variable in the chain is drawn from an inverse Gamma
% distribution with mean and variance equal to the steady state.
%
% The obervations, A, are then computed using the linear relation:
% A = H*W
%
% INPUTS
% W = NMF spectral basis vectors [K,D]
% muinf = steady state mean of the IGMC priors [1,K]
% varinf = steady state variance of the IGMC priors [1,K]
% lam = steady state correlation-coefficient for successive
%       variables in the IGMC priors [1,K]
% T = number of time-steps to draw
% 
% OUTPUTS
% A = NMF data [T,D]
% H = NMF temporal basis functions [T,K]

% see test_randnmf.m for the tests of this function

% draw the temporal priors first
alp = (2 + lam(:).^2 + muinf(:).^2./varinf(:))';

K = length(alp);
H = zeros(T,K);

% set the initial IG parameters so that the have the mean and
% variance of the steady state
alp1 = muinf.^2./varinf+2;
bet1 = (alp1-1).*muinf;

invHt = gamrnd(alp1(:)',bet1(:)');
H(1,:) = 1./invHt;

for t=2:T
  invbet = (alp-1).*(lam.*(H(t-1,:)-muinf)+muinf);
  invHt = gamrnd(alp,1./invbet);
  H(t,:) = 1./invHt;
end

% combine with the spectral basis functions to form the data
A = H*W;