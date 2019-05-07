function y = sigmoid(x,sig_range,c,a)

if nargin<2
    sig_range = [0,20];
end
if nargin<3
    c = 0; 
else
    assert(isscalar(c)==1,'c must be a scalar.') 
end
if nargin<4
    a = 1; 
else
    assert(isscalar(a)==1,'a must be a scalar.') 
end

lo = sig_range(1);
up = sig_range(end);
y = (up-lo)./(1 + exp(-a.*(x-c)))+lo;

end