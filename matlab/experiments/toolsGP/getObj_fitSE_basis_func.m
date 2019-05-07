function [obj,dobj] = getObj_fitSE_basis_func(z,x,lenx,varx,tol)

% function [obj,dobj] = getObj_fitSE_basis_func(z,x,lenx,varx,tol)
%
% Fits a signal x with a squared exponential GP of length scale
% lenx and variance varx using the basis function representation of
% the GP
%
% INPUTS
% z = basis-function coefficients [T,1]
% x = signal to fit [T,1]
% lenx = length scale of GP
% varx = variance of GP
% tol = tolerance parameter which sets the extent of the Gaussian
%       basis functions
%
% OUTPUTS
% obj = objective
% dob = derivative of the objective by z [T,1] 


%%%%%%
% form envelopes from basis functions

T = length(z);

tau = ceil(lenx/sqrt(2)*tol);
bh = sqrt(varx*sqrt(2)/(lenx*sqrt(pi))); % basis function height to
                                         % get correct variance
bas = bh*exp(-1/(lenx^2)*([-tau:1:tau]').^2);
  
xHat =  conv(z,bas,'same');

% compute squared error and add prior
dx = xHat-x;
obj = 1/2*sum(dx.^2)+1/2*sum(z.^2);

% derivatives
dobj = conv(dx,bas,'same')+z;


obj = obj/T;
dobj = dobj/T;