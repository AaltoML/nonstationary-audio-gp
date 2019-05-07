function [obj,dobj] = getObj_nmf_log_norm_inf(x,A,W,vary,lenx,mux,varx,tol)

% function [obj,dobj] = getObj_nmf_log_norm_inf(x,A,W,vary,lenx,mux,varx,tol) 
%
% Returns the objective function for inference in a version of NMF
% with temporal priors which are given by exponentiated squared
% exponential Gaussian Processes
%
% NMF prior (K = number of basis functions can be
% under/over/complete).
%
% A \approx Ahat = H*W
%
% Each of the temporal basis functions in NMF is given by: 
% h(t) = exp(x(t))  x(t) ~ mux + GP(lenx,varx)
%
% The weights are normalised to have sum 1
%
% The cost function used in the likelihood is:
% \sum_{d,t} A_{d,t}./Ahat_{t,d} + \sum_{d,t} \log(Ahat_{t,d}) 
%
%
% INPUTS
% x = vectorised matrices of temporal basis functions [T*K,1]
% A = data [T,D]
% W = spectral basis functions [K,D]
% vary = observation noise [T,D]
% lenx = squared exponential length-scale [1,K]
% mux = steady state mean of the log temporal priors [1,K]
% varx = steady state variance of the log temporal priors [1,K]

[T,D] = size(A);
K = length(x)/(T);

X = reshape(x,[T,K]);
H = zeros([T,K]);

for k=1:K
  tau = ceil(lenx(k)/sqrt(2)*tol);
  bh = sqrt(varx(k)*sqrt(2)/(lenx(k)*sqrt(pi))); % basis function height to
	          				 % get correct variance
  bas = bh*exp(-1/(lenx(k)^2)*([-tau:1:tau]').^2);
					   
  H(:,k) =  exp(conv(X(:,k),bas,'same')+mux(k));
end

Ahat = H*W+vary;

%%%

obj1 = sum(A(:)./Ahat(:)+log(Ahat(:)));

dobj1dA = (Ahat-A)./Ahat.^2;


dobj1dH = dobj1dA*W';
dobj1dX = zeros(T,K);

for k=1:K
  tau = ceil(lenx(k)/sqrt(2)*tol);
  bh = sqrt(varx(k)*sqrt(2)/(lenx(k)*sqrt(pi))); % basis function height to
	          				 % get correct variance
  bas = bh*exp(-1/(lenx(k)^2)*([-tau:1:tau]').^2);
					   
  dobj1dX(:,k) =  conv(dobj1dH(:,k).*H(:,k),bas,'same');
end


obj2 = 1/2*sum(X(:).^2);
dobj2dx = X(:);

% putting the objective function back together
obj = (obj1+obj2)/T;
dobj = [dobj1dX(:)+dobj2dx]/T;


% obj = obj1/T;
% dobj = [dobj1dX(:);dobj1dlogw]/T;





