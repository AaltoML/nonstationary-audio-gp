function [A,H] = randtnmf(W,lenx,mux,varx,T)

% function [A,H] = randtnmf(W,lenx,mux,varx,T)
%
% Draws from the tNMF generative model:
%
% Each of the temporal basis functions in NMF is given by: 
% h(t) = exp(x(t))  x(t) ~ mux + GP(lenx,varx)
%
% The obervations, A, are then computed using the linear relation:
% A = H*W
%
% INPUTS
% W = NMF spectral basis vectors [K,D]
% lenx = squared exponential length-scale [1,K]
% mux = steady state mean of the log temporal priors [1,K]
% varx = steady state variance of the log temporal priors [1,K]
% T = number of time-steps to draw
% 
% OUTPUTS
% A = NMF data [T,D]
% H = NMF temporal basis functions [T,K]
tol = 20;

% draw the temporal priors first
K = length(lenx);
X = randn(T,K);
H = zeros([T,K]);

for k=1:K
  tau = ceil(lenx(k)/sqrt(2)*tol);
  bh = sqrt(varx(k)*sqrt(2)/(lenx(k)*sqrt(pi))); % basis function height to
						 % get correct variance
  bas = bh*exp(-1/(lenx(k)^2)*([-tau:1:tau]').^2);
    
  H(:,k) =  exp(conv(X(:,k),bas,'same')+mux(k));
end

% combine with the spectral basis functions to form the data
A = H*W;