function Y = basis2SEGP(Z,lenx,varx,tol)

% function Y = basis2SEGP(Z,lenx,varx,tol)
%
% converts the basis function representation of a squared
% exponential Gaussian Process to the GP itself
%
% INPUTS
% Z = basis function coefficients [T,K]
% lenx = horizontal length-scale parameter of the GP [K,1]
% varx = vertical scale parameter of GP (marginal variance) [K,1]
% tol = tolerence that sets size of the basis functions (probably
%       should be > 9)
%
% OUTPUTS
% Y = Gaussian process [T,K]

[T,K] = size(Z);
Y = zeros(T,K);

for k=1:K
  tau = ceil(lenx(k)/sqrt(2)*tol);
  bh = sqrt(varx(k)*sqrt(2)/(lenx(k)*sqrt(pi))); % basis function height to
	          				 % get correct variance
  bas = bh*exp(-1/(lenx(k)^2)*([-tau:1:tau]').^2);
					   
  Y(:,k) =  conv(Z(:,k),bas,'same');
end