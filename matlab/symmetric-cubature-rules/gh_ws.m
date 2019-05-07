% Syntax:
%   [W,SX] = gh_ws(m,P,n)
%
% Description:
%   Generate Gauss-Hermite cubature order n based
%   weights W and sigma points SX such that
%
%     int g(x) N(x | m,P) dx ?~ sum W(i) g(SX(i))

% Copyright (C) Simo Särkkä, 2011

function [W,SX] = gh_ws(m,P,n)

    %
    % Form Probabilists' Hermite polynomials of
    % order n-1 and n
    %
    Hpm = 1;
    Hp  = [1 0];
    for i=1:n-1
        tmp = Hp;
        Hp = [Hp 0] - [0 0 i*Hpm];
        Hpm = tmp;
    end

    %
    % Single dimensional weights and points
    %
    xi1 = roots(Hp)';
    W1 = factorial(n)./(n^2*polyval(Hpm,xi1).^2);
 
    d = size(m,1);
    
    %Generate all n^d collections of indexes by
    %transforming numbers 0...n^d-1) into n-base system
    %and by adding 1 to each digit
    num = 0:(n^d-1);
    ind = zeros(d,n^d);
    for i=1:d
        ind(i,:) = rem(num,n)+1;
        num = floor(num / n);
    end

    %Form the sigma points and weights
    L = chol(P)';
    SX = real(L*xi1(ind)+repmat(m,1,size(ind,2)));
    W = real(prod(W1(ind),1)); % ND weights

